#!/usr/bin/env python3
"""Beam-loss budget for a completed run, split the way the literature splits it.

WHY THE SPLIT MATTERS
    Pimazzoni et al., Fus. Eng. Des. 192 (2023) 113621 defines beamlet loss as
        L = 100 (1 - j_acc / j_ext)
    and states (Sec. 3) that "the stripping of the negative ions was considered as a
    reduction of the negative space charge" -- stripped ions are never removed from the
    traced beam there, so L counts ONLY geometric interception of the grids. Their
    stripping loss is quoted separately in the Introduction as ~30%.

    So "beam loss" means two different numbers depending on the source, and comparing a
    total loss against their <5% figure is comparing unlike quantities. This script
    reports them separately.

DEFINITIONS USED
    direct interception  gen-0 primary negative ions arriving on EG or AG1..AGn.
                         This is the quantity comparable to the paper's L.
    stripping loss       negative ions whose trajectory ended in the volume, i.e. were
                         consumed by a gas reaction (row 0, "In_volume").
    transverse escape    negative ions reaching an x/y domain wall. In a single-beamlet
                         domain this is an ARTEFACT: the paper models a 20x64 aperture
                         lattice so such particles enter a neighbouring aperture instead.
                         Enable periodicityBounds to recover it.
    un-extracted         gen-0 negative ions returning to the PG. Never extracted, so
                         excluded from the denominator rather than counted as loss.

USAGE
    python3 tools/beam_loss.py <caseTag> [<caseTag> ...]
    e.g. python3 tools/beam_loss.py heatC_conv
"""
from __future__ import annotations
import sys, math
from pathlib import Path

SUMMARY = Path("Summary")
BEAM_SPECIES = "HM"          # primary negative ion, H- or D- (particle_kind PARTICLE_HM)


def read_breakdown(tag: str):
    p = SUMMARY / f"{tag}_ALL_grid_power_breakdown.txt"
    if not p.is_file():
        raise SystemExit(f"no breakdown for '{tag}': {p}")
    out = []
    for line in p.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 14:
            continue
        out.append(dict(rid=int(f[0]), sp=f[1], gen=int(f[2]), origin=f[3],
                        gross=float(f[4]), cur=abs(float(f[7])), n=int(f[11]),
                        inc=(f[12] == "true"), desc=f[13]))
    return out


def report(tag: str) -> None:
    rows = read_breakdown(tag)
    beam = [r for r in rows if r["sp"] == BEAM_SPECIES]
    if not beam:
        raise SystemExit(f"no {BEAM_SPECIES} rows in '{tag}'")

    grids = {r["desc"] for r in rows if r["inc"]}          # EG + AG1..AGn

    def tally(pred):
        sel = [r for r in beam if pred(r)]
        return sum(r["cur"] for r in sel) * 1e3, sum(r["n"] for r in sel)

    # gen 0 = launched primaries; anything else was made by a reaction
    intercept, n_int = tally(lambda r: r["desc"] in grids and r["gen"] == 0)
    strip,     n_str = tally(lambda r: r["desc"] == "In_volume")
    walls,     n_wal = tally(lambda r: r["desc"] in ("x-min", "x-max", "y-min", "y-max"))
    out,       n_out = tally(lambda r: r["desc"] == "z-max")
    unext,     n_une = tally(lambda r: r["desc"] == "PG")
    # higher-generation negative ions on grids: real load, but not DIRECT interception
    sec_on_grid, _   = tally(lambda r: r["desc"] in grids and r["gen"] != 0)

    extracted = out + strip + walls + intercept + sec_on_grid
    if extracted <= 0:
        raise SystemExit(f"'{tag}': no extracted beam found")

    print(f"\n{'='*74}\n{tag}\n{'='*74}")
    print(f"  un-extracted (returned to PG, excluded)   {unext:9.3f} mA   N={n_une}")
    print(f"  EXTRACTED BEAM                            {extracted:9.3f} mA   = 100 %")
    print(f"  {'-'*70}")
    print(f"  DIRECT INTERCEPTION  (= paper's L)        {intercept:9.3f} mA   {100*intercept/extracted:6.2f} %   N={n_int}")
    print(f"  stripping loss (In_volume)                {strip:9.3f} mA   {100*strip/extracted:6.2f} %   N={n_str}")
    print(f"  transverse escape (ARTEFACT)              {walls:9.3f} mA   {100*walls/extracted:6.2f} %   N={n_wal}")
    if sec_on_grid > 0:
        print(f"  higher-gen neg ions on grids              {sec_on_grid:9.3f} mA   {100*sec_on_grid/extracted:6.2f} %")
    print(f"  transmitted                               {out:9.3f} mA   {100*out/extracted:6.2f} %   N={n_out}")
    print(f"  {'-'*70}")
    phys = intercept + strip
    print(f"  PHYSICAL LOSS (interception + stripping)  {phys:9.3f} mA   {100*phys/extracted:6.2f} %")
    print(f"  reference: paper's L < 1 % at perveance match, <= 5 % at window edges")

    # per-grid split of the direct interception -- tells you WHERE it lands
    per = {}
    for r in beam:
        if r["desc"] in grids and r["gen"] == 0:
            per[r["desc"]] = per.get(r["desc"], 0.0) + r["cur"] * 1e3
    if per:
        print(f"  {'-'*70}\n  direct interception by grid:")
        for k in sorted(per, key=lambda k: -per[k]):
            print(f"     {k:<6} {per[k]:8.4f} mA   {100*per[k]/extracted:6.3f} % of extracted")

    # extracted current density, for placing the case on a perveance scan
    print(f"  {'-'*70}")
    for d_mm in (14.0,):
        a = math.pi * (d_mm/2e3)**2
        print(f"  extracted current density (through {d_mm:.0f} mm aperture) = {extracted*1e-3/a:.0f} A/m2")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    for t in sys.argv[1:]:
        report(t)
