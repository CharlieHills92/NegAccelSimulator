#!/usr/bin/env python3
"""Build a MITICA magnetic-field map from Pimazzoni et al. Fig. 1a.

FORMAT (determined from FieldManager.cpp:98 and meshvectorfield.hpp:159)
    one line per mesh point:   x  y  z  Bx  By  Bz
    coordinates in MILLIMETRES -- FieldManager passes xscale = FIELD_SCALE = 1e-3
    field in TESLA             -- fscale = the config's magneticField.fields[].scale
    ordering: z fastest, then y, then x
    the mesh step h is taken from the FIRST TWO lines, so the mesh must be cubic
    (equal spacing in x, y and z) and strictly coordinate-sorted.

THE APPROXIMATION -- READ THIS
    Fig. 1a is a 1-D plot: B along the axis of ONE beamlet. Reconstructing a 3-D map
    from it requires an assumption about the transverse dependence, and the honest one
    is that there is none: B(x,y,z) = B_axis(z), uniform across the aperture.

    What that preserves: the on-axis field integral, and therefore the total deflection
    an electron accumulates along z. That is what decides WHICH GRID a secondary
    electron lands on, which is the quantity that was wrong with an EG-only field.

    What it discards:
      * the transverse gradient near the embedded magnets (real focusing/defocusing),
      * the ALTERNATING POLARITY between neighbouring apertures. Pimazzoni Sec. 4 and
        Fubiani Sec. III A both model two adjacent beamlets precisely because of this.
        A single-beamlet, transversely-uniform map cannot represent it.

    So this is a substantial improvement over a field that is zero past AG1, but it is
    NOT the NBImag map (paper ref [18]). Treat results as indicative and do not present
    them as a reproduction of the paper's per-grid split.

CURVE READING
    Blue trace is 20*Bx: about -60 mT at z=0 rising smoothly to +10 mT at z=600,
    i.e. Bx from -3.0 mT to +0.5 mT. This is the long-range filter field of the ~4 kA
    flowing in the PG.
    Red trace is By: a strong bipolar feature at the EG reaching about -70 mT, then a
    bipolar feature at each acceleration grid of order +/-40 mT, decaying slowly
    downstream. Amplitudes are read off a printed figure and are the least certain part
    of this reconstruction.
"""
from __future__ import annotations
import argparse, math
from pathlib import Path

# MITICA grid mid-planes [mm], from tools/Templates/MTF_1.json
GRID_Z = {"EG": 21.8, "AG1": 128.0, "AG2": 233.0, "AG3": 338.0, "AG4": 443.0, "GG": 548.0}
# peak |By| of the bipolar feature at each grid [T], read from Fig. 1a
BY_PEAK = {"EG": -0.070, "AG1": 0.040, "AG2": 0.040, "AG3": 0.038, "AG4": 0.036, "GG": 0.034}
BY_WIDTH = 9.0          # mm, half-width of each feature; sets how sharp the spikes are


def by_axis(z_mm: float, shape: str = "unipolar") -> float:
    """By on axis: one feature per grid, summed.

    THE SHAPE IS THE LEAST CERTAIN PART OF THIS RECONSTRUCTION, and it matters, because
    it decides whether each grid imposes a NET deflection or merely a local wiggle:

      unipolar (default)  a Gaussian peak centred on the grid. This is what a magnet
                          embedded in the grid and magnetised transversely actually
                          produces on the aperture axis, and it matches Fig. 1a reading
                          as peaks AT the shaded grid regions. integral By dz != 0, so
                          each grid steers the beam -- which is the point of the EG
                          suppression magnets and of the downstream deflection
                          compensation (paper ref [27]).

      bipolar             the derivative of a Gaussian: antisymmetric, zero at the grid
                          centre, integral ~0. This is the fringe field of a short dipole
                          seen from outside it, and imposes no net deflection.

    Use --shape to switch. Neither is the NBImag map; if you obtain that, use it instead.
    """
    total = 0.0
    for name, z0 in GRID_Z.items():
        u = (z_mm - z0) / BY_WIDTH
        if shape == "bipolar":
            total += BY_PEAK[name] * (-u) * math.exp(0.5 - 0.5 * u * u)
        else:
            total += BY_PEAK[name] * math.exp(-0.5 * u * u)
    return total


def bx_axis(z_mm: float) -> float:
    """Long-range PG filter field: -3.0 mT at z=0 rising smoothly to +0.5 mT at 600 mm."""
    t = max(0.0, min(1.0, z_mm / 600.0))
    return (-3.0 + 3.5 * t ** 0.75) * 1e-3


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=Path("MTFInputFile/PIMAZZONI_fig1a_Bfield.fld"))
    ap.add_argument("--half-xy", type=float, default=40.0, help="half-width in x and y [mm]")
    ap.add_argument("--zmin", type=float, default=0.0)
    ap.add_argument("--zmax", type=float, default=600.0)
    ap.add_argument("--step", type=float, default=2.0,
                    help="cubic mesh step [mm]; must be equal in x, y and z")
    ap.add_argument("--shape", choices=("unipolar", "bipolar"), default="unipolar",
                    help="By feature shape at each grid; see by_axis() (default unipolar)")
    a = ap.parse_args()

    nx = int(round(2 * a.half_xy / a.step)) + 1
    nz = int(round((a.zmax - a.zmin) / a.step)) + 1
    total = nx * nx * nz
    print(f"mesh {nx} x {nx} x {nz} = {total:,} lines, step {a.step} mm, By shape: {a.shape}")

    # B depends only on z, so evaluate once per z-plane instead of per point
    profile = []
    for k in range(nz):
        z = a.zmin + k * a.step
        profile.append((z, bx_axis(z), by_axis(z, a.shape)))

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w") as f:
        for i in range(nx):
            x = -a.half_xy + i * a.step
            for j in range(nx):
                y = -a.half_xy + j * a.step
                for z, bx, by in profile:                    # z fastest
                    f.write(f"{x:.6f} {y:.6f} {z:.6f} {bx:.9E} {by:.9E} 0.000000000E+00\n")
    print(f"wrote {a.out}  ({a.out.stat().st_size/1e6:.1f} MB)")

    print()
    print("on-axis check against Fig. 1a:")
    print("  %-6s %-9s %10s %10s" % ("grid", "z[mm]", "Bx[mT]", "By[mT]"))
    for name, z0 in GRID_Z.items():
        print("  %-6s %-9.1f %10.3f %10.3f" % (name, z0, bx_axis(z0)*1e3, by_axis(z0, a.shape)*1e3))
    pk = max((abs(by_axis(a.zmin + k*0.5, a.shape)), a.zmin + k*0.5)
             for k in range(int((a.zmax-a.zmin)/0.5)))
    print("  peak |By| = %.1f mT at z = %.0f mm   (Fig. 1a: ~70 mT at the EG)" % (pk[0]*1e3, pk[1]))
    integ = sum(math.hypot(bx_axis(z), by_axis(z, a.shape)) * 0.5e-3
                for z in [a.zmin + k*0.5 for k in range(int((a.zmax-a.zmin)/0.5))])
    print("  integral |Bt| dz = %.5f T.m   (EXTfield.fld gives 0.00158)" % integ)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
