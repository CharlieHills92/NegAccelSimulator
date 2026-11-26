#!/usr/bin/env python3
"""Set minimumEnergyEV / maximumEnergyEV on every reaction from its .cross data range.

WHY
    Cross sections are stored as log10(sigma) = sum c_i (log10 E)^i, fitted over the
    POSITIVE samples of the .cross file only (log10(0) is undefined). Outside that window
    the measured cross section is zero, and a polynomial cannot represent zero -- it
    diverges upward instead. h0_process_2, whose data is non-zero only between 5.34 and
    43.3 keV, reaches 1.2e-7 m^2 at 800 keV: thirteen orders of magnitude above an atomic
    cross section, a mean free path of 4e-13 m, hence a collision on EVERY integration
    step. With projectileFate="survive" and four products, that is an unbounded particle
    factory and the dominant cause of secondary-population explosion.

    Bounding is not a fudge. Outside the window the tabulated cross section is zero, so
    returning zero is the correct physics, and every value inside the window is unchanged.

USAGE
    python3 tools/bound_cross_sections.py [--write] config.json [config.json ...]
    Without --write it reports what it would change and touches nothing.
"""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path


def data_window(path: Path) -> tuple[float, float] | None:
    """Energy range over which the .cross file has a strictly positive cross section."""
    if not path.is_file():
        return None
    lo = hi = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            e, s = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        if s > 0.0:
            lo = e if lo is None else min(lo, e)
            hi = e if hi is None else max(hi, e)
    return (lo, hi) if lo is not None else None


def sigma_at(reaction: dict, energy_ev: float) -> float:
    coef = reaction.get("coefficients") or []
    if not coef or energy_ev <= 0:
        return 0.0
    lg = math.log10(energy_ev)
    acc = 0.0
    for c in coef:
        acc = acc * lg + c
    try:
        return 10.0 ** acc
    except OverflowError:
        return float("inf")


def process(cfg: Path, write: bool, probe_ev: float) -> int:
    spec = json.loads(cfg.read_text())
    reactions = (spec.get("gasInteractions") or {}).get("reactions")
    if reactions is None:                      # runtime-format case
        reactions = (spec.get("physics") or {}).get("reactions")
    if not reactions:
        print(f"{cfg}: no reactions found, skipped")
        return 0

    root = cfg.parent
    changed = 0
    print(f"\n{cfg}")
    print("  %-34s %-22s %-24s %s" % ("processId", "data window [eV]", "bounds before -> after", "sigma@probe"))
    for r in reactions:
        src = r.get("sourcePath")
        if not src:
            print("  %-34s %s" % (r.get("processId", "?"), "no sourcePath, skipped"))
            continue
        p = Path(src)
        win = data_window(p if p.is_absolute() else root / p) or data_window(Path.cwd() / src)
        if win is None:
            print("  %-34s %s" % (r.get("processId", "?"), f"cross file not found: {src}"))
            continue
        lo, hi = win
        before_lo, before_hi = r.get("minimumEnergyEV"), r.get("maximumEnergyEV")
        # Intersect: never widen a bound the author set deliberately.
        after_lo = max(float(before_lo), lo) if before_lo is not None else lo
        after_hi = min(float(before_hi), hi) if before_hi is not None else hi
        s_probe = sigma_at(r, probe_ev)
        touched = (before_lo != after_lo) or (before_hi != after_hi)
        if touched:
            changed += 1
            r["minimumEnergyEV"] = after_lo
            r["maximumEnergyEV"] = after_hi
        flag = "  <-- was DIVERGENT" if (s_probe > 1e-18 and (before_hi is None or before_hi > hi)) else ""
        print("  %-34s %-22s %-24s %.3e%s" % (
            r.get("processId", "?"),
            "%.3g .. %.3g" % (lo, hi),
            "%s/%s -> %.4g/%.4g" % (
                ("%.3g" % before_lo) if before_lo is not None else "-",
                ("%.3g" % before_hi) if before_hi is not None else "-",
                after_lo, after_hi),
            s_probe, flag))

    if write and changed:
        cfg.write_text(json.dumps(spec, indent=1) + "\n")
        print(f"  -> WROTE {changed} bounded reaction(s) to {cfg}")
    elif changed:
        print(f"  -> {changed} reaction(s) would change (dry run; pass --write to apply)")
    else:
        print("  -> already fully bounded, nothing to do")
    return changed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("configs", nargs="+", type=Path)
    ap.add_argument("--write", action="store_true", help="apply changes (default: dry run)")
    ap.add_argument("--probe-ev", type=float, default=8.0e5,
                    help="energy at which to report sigma, to expose divergence (default 8e5)")
    a = ap.parse_args()
    total = sum(process(c, a.write, a.probe_ev) for c in a.configs)
    print(f"\n{total} reaction(s) {'bounded' if a.write else 'would be bounded'} in total")
    return 0


if __name__ == "__main__":
    sys.exit(main())
