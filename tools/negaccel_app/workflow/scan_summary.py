"""One combined table of the headline parameters for every case in a scan.

Reads each case's own Summary/ directory and emits a single whitespace-aligned .txt so a
parameter scan can be compared at a glance. One row per case, with:

    scan parameter values     (whatever the scan varied)
    I_beam, Dx, Dy            negative-ion beam current and 1/e divergence at the exit
    jext, jacc                extracted and accelerated current density
    beam loss breakdown       stripping / direct interception / transverse escape
    per-grid power + current  EG and AG1..AGn

followed by two further blocks, one row per (case, grid) and one row per case:

    per-grid power by process interception / volume / surface, per grid and summed
    final particle count      macroparticles left at the end of the run, per species

DEFINITIONS -- these follow the same conventions used elsewhere in this project, and they
matter because "beam loss" means different things in different papers:

  Dx, Dy          1/e divergence = sqrt(2) x current-weighted RMS, i.e. the convention of
                  Pimazzoni et al., Fus. Eng. Des. 192 (2023) 113621, p.4. Taken from the
                  NEGIONBEAM summary, which is negative-ion only.
  interception    gen-0 primary negative ions arriving on a grid. Comparable to that
                  paper's L = 100(1 - j_acc/j_ext), which EXCLUDES stripping because
                  stripped ions there only reduce the space charge.
  stripping       negative ions consumed by a gas reaction: they end inside the volume,
                  so they land in the "In_volume" row.
  transverse      negative ions reaching an x/y domain wall. In a single-beamlet domain
                  this is an ARTEFACT -- with a real aperture lattice they would enter a
                  neighbouring aperture. Reported separately, never folded into loss.
  un-extracted    gen-0 negative ions returning to the PG. Never extracted, so excluded
                  from the denominator instead of counted as loss.
  interception /  the three deposition processes, taken from the breakdown file's own
  volume /        Origin column rather than re-derived here. The C++ writer defines them
  surface         as: primary = injected gen-0 particle arriving on the surface (direct
                  beam interception), volume = gas-phase secondary (gen 1-99), surface =
                  surface secondary (gen >= 101). They partition the gross grid power, so
                  the three columns sum to P_total by construction.
  jext, jacc      average negative-ion current density through one PG aperture,
                  I / (pi r^2) with r the configured apertureRadiusMeters. jext is taken
                  at the configured transmission plane (nominally z = 9 mm, the PG exit)
                  and jacc at the last diagnostic plane. This is the definition Pimazzoni
                  FED 192 (2023) uses -- "the average current density through the PG
                  apertures" -- so L = 100(1 - jacc/jext) is comparable to that paper only
                  if stripping is excluded, which our jacc does NOT do (see intercept).
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterable

from .post_processing import _resolve_case_context, aggregate_case_diagnostics

# Species tag of the primary negative ion. get_particle_name() prefixes with the ACTIVE
# PARTICLE FAMILY (funct.cpp particle_family_prefix), so a hydrogen run writes "HM" and a
# deuterium run writes "DM". Try both rather than assuming hydrogen.
BEAM_SPECIES_TAGS = ("HM", "DM")
WALL_ROWS = {"x-min", "x-max", "y-min", "y-max"}

# Origin values written by DiagnosticsManager::writeGridPowerBreakdown, in report order.
ORIGIN_KEYS = ("primary", "volume", "surface")
ORIGIN_HEADERS = {"primary": "P_intercept[W]", "volume": "P_volume[W]", "surface": "P_surface[W]"}


def _beam_summary(aggregate: dict) -> dict | None:
    """The grid-power summary for whichever negative-ion species this run used."""
    summaries = aggregate.get("gridPowerSummaries") or {}
    for tag in BEAM_SPECIES_TAGS:
        candidate = summaries.get(tag)
        if isinstance(candidate, dict):
            return candidate
    return None


def _beam_exit_row(aggregate: dict[str, Any]) -> dict[str, Any] | None:
    """Last plane of the negative-ion-only summary, i.e. the accelerator exit."""
    ni = aggregate.get("negativeIonSummary")
    if isinstance(ni, dict) and ni.get("rows"):
        return ni["rows"][-1]
    # fall back to the all-species plane file; Dx/Dy there are still beam-only columns
    ds = aggregate.get("diagnosticSummary")
    if isinstance(ds, dict) and ds.get("rows"):
        return ds["rows"][-1]
    return None


def _loss_budget(aggregate: dict[str, Any]) -> dict[str, float]:
    """Negative-ion loss budget in mA, from the beam-species grid-power summary."""
    beam = _beam_summary(aggregate)
    out = dict(extracted=0.0, transmitted=0.0, stripped=0.0,
               intercepted=0.0, transverse=0.0, unextracted=0.0)
    if beam is None:
        return out
    for row in beam.get("rows", []):
        desc = str(row.get("Description", ""))
        cur = abs(float(row.get("Current[A]") or 0.0)) * 1e3      # mA
        if desc == "PG":
            out["unextracted"] += cur
            continue                                              # not part of the beam
        if desc == "z-max":
            out["transmitted"] += cur
        elif desc == "In_volume":
            out["stripped"] += cur
        elif desc in WALL_ROWS:
            out["transverse"] += cur
        elif row.get("IncludeInTotal"):
            out["intercepted"] += cur
        else:
            continue
        out["extracted"] += cur
    return out


def _grid_rows(aggregate: dict[str, Any]) -> list[tuple[str, float, float]]:
    """(name, gross power [W], current [A]) for each grid, in file order."""
    summaries = aggregate.get("gridPowerSummaries") or {}
    all_species = summaries.get("ALL")
    if not isinstance(all_species, dict):
        return []
    return [
        (str(r.get("Description", "?")),
         float(r.get("Power[W]") or 0.0),
         float(r.get("Current[A]") or 0.0))
        for r in all_species.get("rows", [])
        if r.get("IncludeInTotal")
    ]


def _origin_power(aggregate: dict[str, Any]) -> dict[str, dict[str, float]]:
    """{grid name: {origin: gross power [W]}} for the grids counted in the total.

    Origin is read from the breakdown file's own column, not re-derived from generation
    numbers here -- the encoding is shared with the periodic-wrap counters and with
    SURFACE_GENERATION_OFFSET, so the C++ side is the only place that should decode it.

    Restricted to IncludeInTotal rows so the grids match the P_* columns of the main table.
    """
    out: dict[str, dict[str, float]] = {}
    breakdown = (aggregate.get("gridPowerBreakdowns") or {}).get("ALL")
    if not isinstance(breakdown, dict):
        return out
    for row in breakdown.get("rows", []):
        if not row.get("IncludeInTotal"):
            continue
        name = str(row.get("Description", "?"))
        origin = str(row.get("Origin", "")).strip()
        cell = out.setdefault(name, {key: 0.0 for key in ORIGIN_KEYS})
        if origin in cell:
            cell[origin] += float(row.get("GrossPower[W]") or 0.0)
    return out


def _species_counts(aggregate: dict[str, Any]) -> dict[str, int]:
    """Macroparticles present at the end of the run, per species tag.

    Every macroparticle is classified into exactly one grid-power row, so summing the
    Particles column over ALL rows of a per-species file -- grids, domain boundary planes
    and In_volume alike -- gives that species' final population.

    ExcludedParticles is added back on top. Un-extracted primary negative ions are
    deliberately kept out of the Particles column (they are an artefact of the meniscus
    model, not a beam load) but they do still exist in the database, and this is a
    population count rather than a load.
    """
    out: dict[str, int] = {}
    for tag, summary in (aggregate.get("gridPowerSummaries") or {}).items():
        if tag == "ALL" or not isinstance(summary, dict):
            continue
        total = 0
        for row in summary.get("rows", []):
            total += int(row.get("Particles") or 0)
            total += int(row.get("ExcludedParticles") or 0)
        out[tag] = total
    return out


def _current_densities(aggregate: dict[str, Any]) -> dict[str, Any]:
    """Extracted and accelerated negative-ion current density [A/m^2].

    Both are I / (pi r^2) over one PG aperture. The plane for jext is the case's own
    configured transmissionPlaneZPositionMeters (nominally the 9 mm PG exit) and the
    nearest available diagnostic plane to it is used, with the z actually used reported so
    a mismatch is visible rather than silent. jacc is the last plane.

    I[mA] in the plane file is the beam-species (negative-ion) current, not Itot, which is
    what makes these comparable to a published jext/jacc.
    """
    result: dict[str, Any] = dict(j_ext=None, j_acc=None, z_ext_mm=None, z_acc_mm=None,
                                 I_ext_mA=None, I_acc_mA=None, aperture_r_m=None)

    radius: float | None = None
    z_ext_target_mm = 9.0
    config_path = aggregate.get("configPath")
    if config_path:
        try:
            config = json.loads(Path(config_path).read_text(encoding="utf-8"))
            diagnostics = config.get("diagnostics") or {}
            # The authoring schema nests these under diagnostics.summary; older/flatter
            # hand-written configs put them directly on diagnostics. Accept either rather
            # than silently reporting "-" for a case that does define them.
            scopes = [diagnostics.get("summary") or {}, diagnostics]
            for scope in scopes:
                if radius is None:
                    radius = scope.get("apertureRadiusMeters")
                plane = scope.get("transmissionPlaneZPositionMeters")
                if plane is not None:
                    z_ext_target_mm = 1.0e3 * float(plane)
                    break
        except Exception:                                          # noqa: BLE001
            radius = None                                          # fall through, report "-"

    summary = aggregate.get("diagnosticSummary")
    rows = summary.get("rows") if isinstance(summary, dict) else None
    if not rows:
        return result

    def _plane_current(row: dict[str, Any]) -> float | None:
        try:
            return abs(float(row.get("I[mA]")))
        except (TypeError, ValueError):
            return None

    # nearest plane to the configured transmission z
    best = None
    for row in rows:
        try:
            z_mm = float(row.get("z[mm]"))
        except (TypeError, ValueError):
            continue
        delta = abs(z_mm - z_ext_target_mm)
        if best is None or delta < best[0]:
            best = (delta, z_mm, row)

    area = math.pi * radius * radius if radius else None
    result["aperture_r_m"] = radius

    if best is not None:
        _d, z_mm, row = best
        current = _plane_current(row)
        result["z_ext_mm"] = z_mm
        result["I_ext_mA"] = current
        if current is not None and area:
            result["j_ext"] = 1.0e-3 * current / area

    last = rows[-1]
    try:
        result["z_acc_mm"] = float(last.get("z[mm]"))
    except (TypeError, ValueError):
        pass
    current = _plane_current(last)
    result["I_acc_mA"] = current
    if current is not None and area:
        result["j_acc"] = 1.0e-3 * current / area

    return result


def collect_case_row(runtime_json: Path) -> dict[str, Any]:
    """Headline numbers for one case. Raises if its outputs cannot be read."""
    agg = aggregate_case_diagnostics(runtime_json)
    exit_row = _beam_exit_row(agg) or {}
    loss = _loss_budget(agg)
    ext = loss["extracted"] or float("nan")
    return {
        "caseTag": agg.get("caseTag", runtime_json.stem),
        "z_exit_mm": exit_row.get("z[mm]"),
        "I_beam_mA": exit_row.get("I[mA]"),
        "Dx_mrad": exit_row.get("Dx[mrad]"),
        "Dy_mrad": exit_row.get("Dy[mrad]"),
        "extracted_mA": loss["extracted"],
        "transmitted_pct": 100.0 * loss["transmitted"] / ext,
        "stripping_pct": 100.0 * loss["stripped"] / ext,
        "interception_pct": 100.0 * loss["intercepted"] / ext,
        "transverse_pct": 100.0 * loss["transverse"] / ext,
        "unextracted_mA": loss["unextracted"],
        "grids": _grid_rows(agg),
        "originPower": _origin_power(agg),
        "speciesCounts": _species_counts(agg),
        "currentDensity": _current_densities(agg),
    }


def _fmt(value: Any, width: int, decimals: int = 4) -> str:
    if value is None:
        return "-".rjust(width)
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value).rjust(width)
    if f != f:                                    # NaN
        return "-".rjust(width)
    return f"{f:>{width}.{decimals}g}"


def write_scan_summary(
    cases: Iterable[tuple[str, Path]],
    output_path: Path,
    parameter_names: list[str] | None = None,
    parameter_values: dict[str, list[Any]] | None = None,
    title: str | None = None,
) -> Path:
    """Write the combined table.

    cases            (label, runtime_json_path) per case, in scan order
    parameter_names  scan parameter column headers, if the caller has them
    parameter_values label -> list of values, aligned with parameter_names
    title            first header line; defaults to "Combined scan summary"
    """
    parameter_names = parameter_names or []
    parameter_values = parameter_values or {}

    rows: list[dict[str, Any]] = []
    failures: list[tuple[str, str]] = []
    for label, runtime_json in cases:
        try:
            row = collect_case_row(Path(runtime_json))
        except Exception as exc:                                   # noqa: BLE001
            failures.append((label, f"{type(exc).__name__}: {exc}"))
            continue
        row["label"] = label
        rows.append(row)

    # grid names come from the data, so this adapts to MAMuG (EG..AG5), MITICA (+SG), etc.
    grid_names: list[str] = []
    for row in rows:
        for name, _p, _c in row["grids"]:
            if name not in grid_names:
                grid_names.append(name)

    header = ["case".ljust(24)]
    header += [n.rjust(14) for n in parameter_names]
    header += ["I_beam[mA]".rjust(12), "Dx[mrad]".rjust(10), "Dy[mrad]".rjust(10),
               "jext[A/m2]".rjust(12), "jacc[A/m2]".rjust(12),
               "Iext[mA]".rjust(11), "trans[%]".rjust(9), "strip[%]".rjust(9),
               "intercept[%]".rjust(13), "transv[%]".rjust(10)]
    for n in grid_names:
        header += [f"P_{n}[W]".rjust(12), f"I_{n}[A]".rjust(12)]
    header += ["Pgrids[W]".rjust(12)]

    lines: list[str] = []
    lines.append("# " + (title or "Combined scan summary"))
    lines.append("#")
    # The case column is truncated to 24 characters. Only when that actually loses part of a
    # label is the full identity recorded here, so existing scan output -- whose labels are
    # short -- is unchanged.
    if any(len(str(r["label"])) > 24 for r in rows):
        lines.append("# cases:")
        for r in rows:
            tag = r.get("caseTag") or ""
            lines.append(f"#   {r['label']}"
                         + (f"  =  {tag}" if str(tag) != str(r["label"]) else ""))
        lines.append("#")
    lines.append("# Dx,Dy       1/e divergence (sqrt(2) x current-weighted RMS), negative ions only,")
    lines.append("#             at the last diagnostic plane. Convention of Pimazzoni FED 192 (2023).")
    lines.append("# Iext        extracted negative-ion current (all destinations except the PG).")
    lines.append("# trans       transmitted out of the accelerator, % of Iext.")
    lines.append("# strip       lost to gas stripping (ends in volume), % of Iext.")
    lines.append("# intercept   gen-0 primary negative ions onto a grid, % of Iext.")
    lines.append("#             This is the quantity comparable to the paper's L = 100(1-jacc/jext),")
    lines.append("#             which EXCLUDES stripping. Do not add strip and intercept when")
    lines.append("#             comparing against published L values.")
    lines.append("# transv      reached an x/y domain wall. ARTEFACT of a single-beamlet domain;")
    lines.append("#             with an aperture lattice these enter a neighbouring aperture.")
    lines.append("# jext        extracted negative-ion current density at the transmission plane")
    lines.append("#             (nominally z = 9 mm), I/(pi r^2) over one PG aperture.")
    lines.append("# jacc        accelerated negative-ion current density at the last plane.")
    lines.append("#             NOTE jacc is depleted by stripping in this code, whereas Pimazzoni")
    lines.append("#             keeps stripped ions in the traced beam. 1-jacc/jext is therefore")
    lines.append("#             stripping+interception here, NOT that paper's interception-only L.")
    lines.append("# P_*, I_*    per-grid gross arriving power and current (all species).")
    lines.append("# Pgrids      sum of P_* over the grids.")
    lines.append("#")
    lines.append("# " + " ".join(header))

    for row in rows:
        cells = [str(row["label"])[:24].ljust(24)]
        vals = parameter_values.get(row["label"], [])
        for i in range(len(parameter_names)):
            cells.append(_fmt(vals[i] if i < len(vals) else None, 14, 6))
        cells.append(_fmt(row["I_beam_mA"], 12, 6))
        cells.append(_fmt(row["Dx_mrad"], 10, 4))
        cells.append(_fmt(row["Dy_mrad"], 10, 4))
        cells.append(_fmt(row["currentDensity"].get("j_ext"), 12, 6))
        cells.append(_fmt(row["currentDensity"].get("j_acc"), 12, 6))
        cells.append(_fmt(row["extracted_mA"], 11, 6))
        cells.append(_fmt(row["transmitted_pct"], 9, 4))
        cells.append(_fmt(row["stripping_pct"], 9, 4))
        cells.append(_fmt(row["interception_pct"], 13, 4))
        cells.append(_fmt(row["transverse_pct"], 10, 4))
        by_name = {n: (p, c) for n, p, c in row["grids"]}
        total = 0.0
        for n in grid_names:
            p, c = by_name.get(n, (None, None))
            if p is not None:
                total += p
            cells.append(_fmt(p, 12, 6))
            cells.append(_fmt(c, 12, 5))
        cells.append(_fmt(total, 12, 6))
        lines.append("  " + " ".join(cells))

    # ---------------------------------------------------------------------------------
    # Block 2: per-grid power split by deposition process.
    #
    # The three origins partition the gross power, so P_intercept+P_volume+P_surface is
    # P_total by construction and the TOTAL line matches Pgrids in the main table. Kept as
    # its own block rather than more main-table columns because it is 4 numbers per grid,
    # which would make one row per case unreadably wide.
    # ---------------------------------------------------------------------------------
    lines.append("#")
    lines.append("# ==== Per-grid power by deposition process [W] ====")
    lines.append("# intercept = direct beam interception (injected gen-0 particle onto the grid)")
    lines.append("# volume    = gas-phase secondaries (stripping/ionisation products)")
    lines.append("# surface   = surface secondaries (EAMCC backscatter and secondary electrons)")
    lines.append("# The three partition the gross load, so they sum to P_total exactly.")
    lines.append("# TOTAL is the sum over all grids and matches Pgrids above.")
    lines.append("#")
    origin_header = ["case".ljust(24), "grid".rjust(8)]
    origin_header += [ORIGIN_HEADERS[k].rjust(15) for k in ORIGIN_KEYS]
    origin_header += ["P_total[W]".rjust(14)]
    lines.append("# " + " ".join(origin_header))
    for row in rows:
        per_grid = row.get("originPower") or {}
        running = {key: 0.0 for key in ORIGIN_KEYS}
        for name in grid_names:
            cell = per_grid.get(name)
            if cell is None:
                continue
            for key in ORIGIN_KEYS:
                running[key] += cell.get(key, 0.0)
            values = [cell.get(key, 0.0) for key in ORIGIN_KEYS]
            out = [str(row["label"])[:24].ljust(24), str(name).rjust(8)]
            out += [_fmt(v, 15, 6) for v in values]
            out.append(_fmt(sum(values), 14, 6))
            lines.append("  " + " ".join(out))
        out = [str(row["label"])[:24].ljust(24), "TOTAL".rjust(8)]
        out += [_fmt(running[key], 15, 6) for key in ORIGIN_KEYS]
        out.append(_fmt(sum(running.values()), 14, 6))
        lines.append("  " + " ".join(out))

    # ---------------------------------------------------------------------------------
    # Block 3: final macroparticle population per species.
    # ---------------------------------------------------------------------------------
    species_tags: list[str] = []
    for row in rows:
        for tag in sorted((row.get("speciesCounts") or {}).keys()):
            if tag not in species_tags:
                species_tags.append(tag)

    lines.append("#")
    lines.append("# ==== Macroparticles at end of run, by species ====")
    lines.append("# Every macroparticle lands in exactly one grid-power row, so this counts the")
    lines.append("# whole final database: grids, domain boundary planes and In_volume together.")
    lines.append("# Un-extracted primary negative ions are included here (they exist in the")
    lines.append("# database) even though they are excluded from the power/current columns.")
    lines.append("#")
    count_header = ["case".ljust(24)] + [t.rjust(12) for t in species_tags] + ["total".rjust(12)]
    lines.append("# " + " ".join(count_header))
    for row in rows:
        counts = row.get("speciesCounts") or {}
        out = [str(row["label"])[:24].ljust(24)]
        total_count = 0
        for tag in species_tags:
            value = counts.get(tag)
            total_count += int(value or 0)
            out.append(("-" if value is None else str(int(value))).rjust(12))
        out.append(str(total_count).rjust(12))
        lines.append("  " + " ".join(out))

    if failures:
        lines.append("#")
        lines.append("# cases skipped (outputs missing or unreadable):")
        for label, why in failures:
            lines.append(f"#   {label}: {why}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def default_case_summary_path(runtime_json: Path) -> Path:
    """Where a solo run's summary belongs: its own Summary/ directory, beside the per-species
    files the table is built from."""
    context = _resolve_case_context(Path(runtime_json))
    return Path(context["summaryDirectory"]) / f"{context['caseTag']}_case_summary.txt"


def write_case_summary(
    runtime_json: Path,
    output_path: Path | None = None,
    label: str | None = None,
) -> Path:
    """Write the same table as write_scan_summary, for ONE run that was not part of a scan.

    A solo run is just a scan of one case with nothing varying, so this delegates rather
    than duplicating the table: every column, block and definition stays identical, and a
    solo run can be compared against a scan row directly.

    runtime_json  the case's runtime .json
    output_path   defaults to <run root>/Summary/<caseTag>_case_summary.txt
    label         row label, defaults to the case tag
    """
    runtime_json = Path(runtime_json)
    if output_path is None:
        output_path = default_case_summary_path(runtime_json)
    if label is None:
        label = _resolve_case_context(runtime_json)["caseTag"]
    return write_scan_summary(
        [(label, runtime_json)],
        Path(output_path),
        title="Single-case summary (run outside the scan manager)",
    )
