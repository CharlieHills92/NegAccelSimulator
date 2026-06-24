"""Post-processing helpers for retained NegAccel diagnostics outputs."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from .common import DEFAULT_OUTPUTS, WorkflowError, load_json, require_object, require_string, write_json


OUTPUT_CATEGORY_GENERAL = "General outputs"
OUTPUT_CATEGORY_ITERATION = "Outputs per iteration"
OUTPUT_CATEGORY_FINAL = "Outputs at end of simulation"

ARTIFACT_TYPE_VTK = "VTK"
ARTIFACT_TYPE_PNG = "PNG plots"
ARTIFACT_TYPE_TEXT = "Text outputs"

_GEOMETRY_VTK_SUFFIXES = (
    "_geometry.vtk",
    "_geometry_detailed.vtk",
    "_solids.vtk",
)
_ITERATION_TOKEN_PATTERN = re.compile(r"_(?:it|iter|iteration)(?:_|-|\d)", re.IGNORECASE)


def resolve_runtime_output_directory(
    runtime_case: dict[str, Any],
    runtime_path: Path,
    output_key: str,
    default_directory: str,
) -> Path:
    outputs = runtime_case.get("outputs", {})
    output_config = outputs.get(output_key, {})
    directory_name = default_directory
    if isinstance(output_config, dict):
        configured_name = output_config.get("directory")
        if isinstance(configured_name, str) and configured_name.strip():
            directory_name = configured_name.strip()

    configured_root = outputs.get("rootDirectory")
    root_directory = runtime_path.parent
    if isinstance(configured_root, str) and configured_root.strip():
        root_directory = _resolve_path(configured_root, runtime_path.parent)

    directory = Path(directory_name)
    if not directory.is_absolute():
        directory = root_directory / directory
    return directory.resolve()


def format_runtime_setup_summary(runtime_case: dict[str, Any], runtime_path: Path) -> str:
    metadata = runtime_case.get("metadata", {}) if isinstance(runtime_case, dict) else {}
    geometry = runtime_case.get("geometry", {}) if isinstance(runtime_case, dict) else {}
    particle_sources = runtime_case.get("particleSources")
    simulation = runtime_case.get("simulation", {}) if isinstance(runtime_case, dict) else {}
    diagnostics = runtime_case.get("diagnostics", {}) if isinstance(runtime_case, dict) else {}
    outputs = runtime_case.get("outputs", {}) if isinstance(runtime_case, dict) else {}

    summary_lines: list[str] = []
    summary_lines.append("NegAccel Runtime Setup Summary")
    summary_lines.append("=" * 29)
    summary_lines.append("")

    case_tag = metadata.get("caseTag", "unknown")
    summary_lines.extend(
        [
            f"Case tag: {case_tag}",
            f"Runtime JSON: {runtime_path.resolve()}",
        ]
    )
    title = metadata.get("title")
    if isinstance(title, str) and title.strip():
        summary_lines.append(f"Title: {title.strip()}")
    description = metadata.get("description")
    if isinstance(description, str) and description.strip():
        summary_lines.append(f"Description: {description.strip()}")
    summary_lines.append("")

    domain = geometry.get("domain", {}) if isinstance(geometry, dict) else {}
    mesh = geometry.get("mesh", {}) if isinstance(geometry, dict) else {}
    solids = geometry.get("solids", []) if isinstance(geometry, dict) else []
    summary_lines.append("Geometry")
    summary_lines.append("--------")
    summary_lines.append(f"Source mode: {_coerce_text(geometry.get('source', {}).get('mode'), fallback='unknown')}")
    summary_lines.append(
        "Domain [m]: "
        f"x={_format_number(domain.get('xSizeMeters'))}, "
        f"y={_format_number(domain.get('ySizeMeters'))}, "
        f"z={_format_number(domain.get('zSizeMeters'))}, "
        f"z-start={_format_number(domain.get('zStartMeters'), default='0')}"
    )
    summary_lines.append(f"Mesh size [m]: {_format_number(mesh.get('sizeMeters'))}")
    summary_lines.append(f"Solids: {len(solids) if isinstance(solids, list) else 0}")
    summary_lines.append("")

    summary_lines.append("Particle Sources")
    summary_lines.append("----------------")
    source_lines = _format_particle_sources(particle_sources)
    summary_lines.extend(source_lines or ["No particle sources available in the runtime case."])
    summary_lines.append("")

    summary_lines.append("Simulation")
    summary_lines.append("----------")
    summary_lines.extend(_format_simulation_settings(simulation))
    summary_lines.append("")

    summary_lines.append("Diagnostics")
    summary_lines.append("-----------")
    summary_lines.extend(_format_diagnostics_settings(diagnostics))
    summary_lines.append("")

    summary_lines.append("Outputs")
    summary_lines.append("-------")
    summary_lines.extend(_format_outputs_settings(outputs, runtime_path))
    summary_lines.append("")

    return "\n".join(summary_lines).rstrip() + "\n"


def write_runtime_setup_summary(runtime_case: dict[str, Any], runtime_path: Path) -> Path | None:
    outputs = runtime_case.get("outputs", {}) if isinstance(runtime_case, dict) else {}
    summary_config = outputs.get("summary", {}) if isinstance(outputs, dict) else {}
    if isinstance(summary_config, dict) and not bool(summary_config.get("enabled", True)):
        return None

    summary_directory = resolve_runtime_output_directory(runtime_case, runtime_path, "summary", "Summary")
    summary_directory.mkdir(parents=True, exist_ok=True)
    case_tag = require_string(runtime_case.get("metadata", {}), "caseTag", "runtime_case.metadata")
    summary_path = summary_directory / f"{case_tag}_setup_summary.txt"
    summary_path.write_text(format_runtime_setup_summary(runtime_case, runtime_path), encoding="utf-8")
    return summary_path


def build_output_catalog(runtime_case: dict[str, Any], runtime_path: Path) -> list[dict[str, Any]]:
    catalog: list[dict[str, Any]] = []
    for artifact_type, output_key, default_directory, pattern in (
        (ARTIFACT_TYPE_VTK, "vtk", "VTK", "*.vtk"),
        (ARTIFACT_TYPE_PNG, "plots", "Plots", "*.png"),
        (ARTIFACT_TYPE_TEXT, "summary", "Summary", "*.txt"),
    ):
        directory = resolve_runtime_output_directory(runtime_case, runtime_path, output_key, default_directory)
        if not directory.exists():
            continue
        for path in sorted(directory.glob(pattern)):
            if not path.is_file():
                continue
            catalog.append(_classify_output_artifact(path, artifact_type))

    return sorted(
        catalog,
        key=lambda entry: (
            _category_rank(entry["category"]),
            _artifact_type_rank(entry["artifactType"]),
            int(entry.get("iteration") or -1),
            entry["name"],
        ),
    )


def _classify_output_artifact(path: Path, artifact_type: str) -> dict[str, Any]:
    name = path.name
    lowered = name.lower()
    category = OUTPUT_CATEGORY_FINAL
    kind = "raw"
    iteration = _extract_iteration_number(name)

    if lowered.endswith("_setup_summary.txt"):
        category = OUTPUT_CATEGORY_GENERAL
        kind = "setup-summary"
    elif artifact_type == ARTIFACT_TYPE_VTK and any(lowered.endswith(suffix) for suffix in _GEOMETRY_VTK_SUFFIXES):
        category = OUTPUT_CATEGORY_GENERAL
        kind = "geometry-vtk"
    elif _is_iteration_artifact(name):
        category = OUTPUT_CATEGORY_ITERATION
        kind = "iteration-artifact"
    elif artifact_type == ARTIFACT_TYPE_VTK:
        if lowered.endswith("_trajectories.vtk"):
            kind = "trajectories-vtk"
        elif lowered.endswith("_potential.vtk"):
            kind = "potential-vtk"
        elif lowered.endswith("_scharge.vtk"):
            kind = "space-charge-vtk"
        elif lowered.endswith("_simulation_state.vtk"):
            kind = "simulation-state-vtk"
        else:
            kind = "vtk"
    elif artifact_type == ARTIFACT_TYPE_PNG:
        kind = "plot-png"
    elif lowered.endswith("_diagnostic_summary.txt"):
        kind = "diagnostic-summary"
    elif lowered.endswith("_grid_power_summary.txt"):
        kind = "grid-power-summary"
    elif lowered.endswith("final_map_outside.txt"):
        kind = "emitter-map"

    return {
        "path": path.as_posix(),
        "name": name,
        "artifactType": artifact_type,
        "category": category,
        "kind": kind,
        "iteration": iteration,
    }


def _is_iteration_artifact(name: str) -> bool:
    lowered = name.lower()
    if lowered.endswith("_it_eg.txt") or lowered.endswith("_it_out.txt"):
        return True
    return bool(_ITERATION_TOKEN_PATTERN.search(lowered))


def _extract_iteration_number(name: str) -> int | None:
    match = re.search(r"_(?:it|iter|iteration)[_-]?(\d+)", name, re.IGNORECASE)
    return int(match.group(1)) if match else None


def _category_rank(category: str) -> int:
    if category == OUTPUT_CATEGORY_GENERAL:
        return 0
    if category == OUTPUT_CATEGORY_ITERATION:
        return 1
    return 2


def _artifact_type_rank(artifact_type: str) -> int:
    if artifact_type == ARTIFACT_TYPE_VTK:
        return 0
    if artifact_type == ARTIFACT_TYPE_PNG:
        return 1
    return 2


def _coerce_text(value: Any, *, fallback: str) -> str:
    return value.strip() if isinstance(value, str) and value.strip() else fallback


def _format_number(value: Any, *, default: str = "n/a") -> str:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return default
    return f"{float(value):g}"


def _format_vector(value: Any) -> str:
    if not isinstance(value, list) or not value:
        return "n/a"
    formatted = []
    for component in value:
        formatted.append(_format_number(component))
    return "[" + ", ".join(formatted) + "]"


def _format_particle_sources(particle_sources: Any) -> list[str]:
    if isinstance(particle_sources, dict):
        negative_ion_beam = particle_sources.get("negativeIonBeam")
        if isinstance(negative_ion_beam, dict):
            return [
                "Compact negative-ion beam source:",
                f"  Species: {_coerce_text(negative_ion_beam.get('species'), fallback='unknown')}",
                f"  Charge state: {_format_number(negative_ion_beam.get('chargeState'))}",
                f"  Mass [u]: {_format_number(negative_ion_beam.get('massU'))}",
                f"  Current density [A/m^2]: {_format_number(negative_ion_beam.get('currentDensityAm2'))}",
                f"  Axial energy [eV]: {_format_number(negative_ion_beam.get('axialEnergyEV'))}",
                f"  Source model: {_coerce_text(negative_ion_beam.get('sourceModel'), fallback='unknown')}",
            ]
        return []

    if not isinstance(particle_sources, list) or not particle_sources:
        return []

    lines = [f"Explicit sources: {len(particle_sources)}"]
    for index, source in enumerate(particle_sources, start=1):
        if not isinstance(source, dict):
            continue
        uniform = source.get("uniform", {}) if isinstance(source.get("uniform"), dict) else {}
        lines.extend(
            [
                f"{index}. {source.get('name') or source.get('id') or 'source'}",
                f"   Kind: {_coerce_text(source.get('kind'), fallback='unknown')}",
                f"   Particle count: {source.get('particleCount', 'n/a')}",
                f"   Current density [A/m^2]: {_format_number(source.get('currentDensityAm2'))}",
                f"   Axial energy [eV]: {_format_number(source.get('axialEnergyEV'))}",
                f"   Center [m]: {_format_vector(uniform.get('centerMeters'))}",
                f"   Main direction: {_format_vector(uniform.get('mainDirection'))}",
                f"   Width x height [m]: {_format_number(uniform.get('widthMeters'))} x {_format_number(uniform.get('heightMeters'))}",
            ]
        )
    return lines


def _format_simulation_settings(simulation: Any) -> list[str]:
    if not isinstance(simulation, dict):
        return ["Simulation settings unavailable."]

    lines = [
        f"Particle count: {simulation.get('particleCount', 'n/a')}",
        f"Iterations: {simulation.get('iterations', 'n/a')}",
    ]

    solver = simulation.get("solver") if isinstance(simulation.get("solver"), dict) else {}
    if solver:
        solver_type = _coerce_text(solver.get("type"), fallback="unknown")
        lines.append(f"Solver: {solver_type}")
        if solver_type == "multigrid":
            multigrid = solver.get("multigrid") if isinstance(solver.get("multigrid"), dict) else {}
            if multigrid:
                lines.append(
                    "  Multigrid: "
                    f"levels={multigrid.get('levels', 'n/a')}, "
                    f"tolerance={_format_number(multigrid.get('mgTolerance'))}, "
                    f"maxCycles={multigrid.get('maxCycles', 'n/a')}"
                )
        elif solver_type == "bicgstab":
            bicgstab = solver.get("bicgstab") if isinstance(solver.get("bicgstab"), dict) else {}
            if bicgstab:
                lines.append(
                    "  BiCGSTAB: "
                    f"eps={_format_number(bicgstab.get('eps'))}, "
                    f"maxIterations={bicgstab.get('maxIterations', 'n/a')}"
                )

    convergence = simulation.get("convergence") if isinstance(simulation.get("convergence"), dict) else {}
    if convergence:
        lines.append(
            "Convergence: "
            f"jTolerance={_format_number(convergence.get('jTolerance'))}, "
            f"alphaCoeff={_format_number(convergence.get('alphaCoeff'))}"
        )

    return lines


def _format_diagnostics_settings(diagnostics: Any) -> list[str]:
    if not isinstance(diagnostics, dict):
        return ["Diagnostics settings unavailable."]

    planes = diagnostics.get("planes") if isinstance(diagnostics.get("planes"), dict) else {}
    summary = diagnostics.get("summary") if isinstance(diagnostics.get("summary"), dict) else {}
    species = diagnostics.get("species") if isinstance(diagnostics.get("species"), dict) else {}
    grid_power = diagnostics.get("gridPower") if isinstance(diagnostics.get("gridPower"), dict) else {}

    lines = [
        f"Along-Z sample planes [m]: {_format_vector(planes.get('sampleZPositionsMeters'))}",
        f"Summary plane z [m]: {_format_number(planes.get('summaryZPositionMeters'))}",
        f"Emitter export plane z [m]: {_format_number(planes.get('emitterExportZPositionMeters'))}",
        f"Transmission plane z [m]: {_format_number(summary.get('transmissionPlaneZPositionMeters'))}",
        f"Aperture radius [m]: {_format_number(summary.get('apertureRadiusMeters'))}",
        "Per-species diagnostics: "
        f"diagnostics={bool(species.get('writePerSpeciesDiagnostics', False))}, "
        f"gridPower={bool(species.get('writePerSpeciesGridPower', False))}, "
        f"plots={bool(species.get('writePerSpeciesPlots', False))}, "
        f"negativeIonSummary={bool(species.get('writeNegativeIonSummary', False))}",
        f"Grid-power ranges: {len(grid_power.get('ranges', [])) if isinstance(grid_power.get('ranges'), list) else 0}",
    ]
    return lines


def _format_outputs_settings(outputs: Any, runtime_path: Path) -> list[str]:
    if not isinstance(outputs, dict):
        return ["Output settings unavailable."]

    lines = [f"Root directory: {_coerce_text(outputs.get('rootDirectory'), fallback=str(runtime_path.parent))}"]
    for output_key, default_directory in (
        ("summary", "Summary"),
        ("plots", "Plots"),
        ("data", "Data"),
        ("vtk", "VTK"),
    ):
        output_config = outputs.get(output_key)
        if not isinstance(output_config, dict):
            continue
        lines.append(
            f"{output_key.capitalize()}: enabled={bool(output_config.get('enabled', True))}, "
            f"directory={_coerce_text(output_config.get('directory'), fallback=default_directory)}"
        )
    vtk = outputs.get("vtk") if isinstance(outputs.get("vtk"), dict) else {}
    if vtk:
        lines.append(
            "VTK exports: "
            f"geometry={bool(vtk.get('exportGeometry', False))}, "
            f"simulationState={bool(vtk.get('exportSimulationState', False))}, "
            f"tracedParticles={bool(vtk.get('exportTracedParticles', False))}"
        )
    logging = outputs.get("logging") if isinstance(outputs.get("logging"), dict) else {}
    if logging:
        lines.append(
            "Logging: "
            f"consoleLevel={_coerce_text(logging.get('consoleLevel'), fallback='n/a')}, "
            f"fileLevel={_coerce_text(logging.get('fileLevel'), fallback='n/a')}, "
            f"captureStdout={bool(logging.get('captureStdout', False))}, "
            f"debugArtifacts={bool(logging.get('writeDebugArtifacts', False))}"
        )
    iteration = outputs.get("iteration") if isinstance(outputs.get("iteration"), dict) else {}
    if iteration:
        lines.append(
            "Iteration outputs: "
            f"enabled={bool(iteration.get('enabled', True))}, "
            f"everyNIterations={iteration.get('everyNIterations', 'n/a')}, "
            f"planeDiagnostics={bool(iteration.get('exportPlaneDiagnostics', False))}, "
            f"simulationState={bool(iteration.get('exportSimulationState', False))}, "
            f"tracedParticles={bool(iteration.get('exportTracedParticles', False))}"
        )
        lines.append(
            f"Iteration plane overrides [m]: {_format_vector(iteration.get('planeZPositionsMeters'))}"
        )
    return lines


def _parse_numeric_token(token: str, context: str) -> int | float:
    try:
        value = float(token)
    except ValueError as exc:
        raise WorkflowError(f"Invalid numeric token '{token}' in {context}") from exc

    if "." not in token and "e" not in token.lower():
        return int(value)
    return value


def _resolve_path(value: str, base_dir: Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    else:
        path = path.resolve()
    return path


def _resolve_case_config_path(case_source: Path) -> Path:
    if case_source.is_file():
        return case_source.resolve()

    if not case_source.is_dir():
        raise WorkflowError(f"Case source does not exist: {case_source}")

    expected = case_source / f"{case_source.name}.json"
    if expected.is_file():
        return expected.resolve()

    json_files = sorted(path for path in case_source.glob("*.json") if path.is_file())
    if len(json_files) == 1:
        return json_files[0].resolve()

    raise WorkflowError(
        f"Could not resolve a unique case JSON inside {case_source}; expected {expected.name} or one JSON file"
    )


def _resolve_case_context(case_source: Path) -> dict[str, Any]:
    config_path = _resolve_case_config_path(case_source)
    case_config = load_json(config_path)
    if not isinstance(case_config, dict):
        raise WorkflowError(f"Case configuration must be a JSON object: {config_path}")

    metadata = require_object(case_config, "metadata", str(config_path))
    case_tag = require_string(metadata, "caseTag", f"{config_path}.metadata")

    outputs = require_object(case_config, "outputs", str(config_path))
    configured_root = outputs.get("rootDirectory")
    if isinstance(configured_root, str) and configured_root:
        root_directory = _resolve_path(configured_root, config_path.parent)
    else:
        root_directory = config_path.parent.resolve()

    summary_dir_name = DEFAULT_OUTPUTS["summary"]["directory"]
    summary_block = outputs.get("summary")
    if isinstance(summary_block, dict):
        configured_summary_dir = summary_block.get("directory")
        if isinstance(configured_summary_dir, str) and configured_summary_dir:
            summary_dir_name = configured_summary_dir

    summary_directory = _resolve_path(summary_dir_name, root_directory)
    return {
        "configPath": config_path,
        "config": case_config,
        "caseTag": case_tag,
        "rootDirectory": root_directory,
        "summaryDirectory": summary_directory,
    }


def parse_diagnostic_summary_txt(filepath: Path) -> dict[str, Any]:
    if not filepath.is_file():
        raise WorkflowError(f"Diagnostic summary not found: {filepath}")

    columns: list[str] | None = None
    rows: list[dict[str, int | float]] = []
    for line_number, raw_line in enumerate(filepath.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith("#"):
            header_candidate = line[1:].strip()
            if header_candidate.startswith("it "):
                columns = header_candidate.split()
            continue

        if columns is None:
            raise WorkflowError(f"Missing diagnostic header before data in {filepath} at line {line_number}")

        tokens = line.split()
        if len(tokens) != len(columns):
            raise WorkflowError(
                f"Diagnostic row length mismatch in {filepath} at line {line_number}: "
                f"expected {len(columns)} columns, found {len(tokens)}"
            )

        row: dict[str, int | float] = {}
        for column, token in zip(columns, tokens):
            row[column] = _parse_numeric_token(token, f"{filepath}:{line_number}")
        rows.append(row)

    if columns is None:
        raise WorkflowError(f"No diagnostic header found in {filepath}")

    return {
        "path": filepath.as_posix(),
        "columns": columns,
        "rows": rows,
    }


def parse_grid_power_summary_txt(filepath: Path) -> dict[str, Any]:
    if not filepath.is_file():
        raise WorkflowError(f"Grid power summary not found: {filepath}")

    rows: list[dict[str, Any]] = []
    total_included_beam_power_watts: float | None = None
    columns = ["ID", "Power[W]", "Current[A]", "Particles", "IncludeInTotal", "Description"]

    for line_number, raw_line in enumerate(filepath.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith("#"):
            summary_line = line[1:].strip()
            if summary_line.startswith("Total beam power"):
                numeric_token = summary_line.rsplit(":", 1)[1].strip().split()[0]
                total_included_beam_power_watts = float(numeric_token)
            continue

        parts = [part.strip() for part in raw_line.rstrip("\n").split("\t")]
        if len(parts) == 6:
            rows.append(
                {
                    "ID": int(parts[0]),
                    "Power[W]": float(parts[1]),
                    "Current[A]": float(parts[2]),
                    "Particles": int(parts[3]),
                    "IncludeInTotal": parts[4].lower() == "true",
                    "Description": parts[5],
                }
            )
            continue

        if len(parts) == 5:
            rows.append(
                {
                    "ID": int(parts[0]),
                    "Power[W]": float(parts[1]),
                    "Current[A]": float(parts[2]),
                    "Particles": int(parts[3]),
                    "IncludeInTotal": None,
                    "Description": parts[4],
                }
            )
            continue

        raise WorkflowError(
            f"Grid power row length mismatch in {filepath} at line {line_number}: expected 5 or 6 tab-delimited fields"
        )

    return {
        "path": filepath.as_posix(),
        "columns": columns,
        "rows": rows,
        "totalIncludedBeamPowerWatts": total_included_beam_power_watts,
    }


def aggregate_case_diagnostics(case_source: Path) -> dict[str, Any]:
    context = _resolve_case_context(case_source)
    case_tag = context["caseTag"]
    summary_directory = context["summaryDirectory"]

    diagnostic_summary_path = summary_directory / f"{case_tag}_diagnostic_summary.txt"
    negative_ion_summary_path = summary_directory / f"{case_tag}_NEGIONBEAM_diagnostic_summary.txt"
    grid_power_suffix = "_grid_power_summary.txt"

    diagnostic_summary = parse_diagnostic_summary_txt(diagnostic_summary_path)
    negative_ion_summary = (
        parse_diagnostic_summary_txt(negative_ion_summary_path)
        if negative_ion_summary_path.is_file()
        else None
    )

    grid_power_summaries: dict[str, dict[str, Any]] = {}
    for filepath in sorted(summary_directory.glob(f"{case_tag}_*{grid_power_suffix}")):
        species_tag = filepath.name[len(case_tag) + 1 : -len(grid_power_suffix)]
        grid_power_summaries[species_tag] = parse_grid_power_summary_txt(filepath)

    return {
        "caseTag": case_tag,
        "configPath": context["configPath"].as_posix(),
        "rootDirectory": context["rootDirectory"].as_posix(),
        "summaryDirectory": summary_directory.as_posix(),
        "diagnosticSummary": diagnostic_summary,
        "negativeIonSummary": negative_ion_summary,
        "gridPowerSummaries": grid_power_summaries,
    }


def write_case_diagnostics_aggregation(case_source: Path, output_path: Path | None = None) -> Path:
    payload = aggregate_case_diagnostics(case_source)
    resolved_output = output_path or (
        Path(payload["summaryDirectory"]) / f"{payload['caseTag']}_diagnostics_aggregated.json"
    )
    write_json(resolved_output, payload)
    return resolved_output


def synthesize_scan_results(manifest_path: Path) -> dict[str, Any]:
    manifest = load_json(manifest_path)
    if not isinstance(manifest, dict):
        raise WorkflowError(f"Scan manifest must be a JSON object: {manifest_path}")

    cases = manifest.get("cases")
    if not isinstance(cases, list) or not cases:
        raise WorkflowError(f"Scan manifest requires a non-empty 'cases' array: {manifest_path}")

    case_summaries: list[dict[str, Any]] = []
    diagnostic_rows: list[dict[str, Any]] = []
    negative_ion_rows: list[dict[str, Any]] = []
    grid_power_rows: list[dict[str, Any]] = []

    for case in cases:
        if not isinstance(case, dict):
            raise WorkflowError(f"Each manifest case must be an object: {manifest_path}")

        config_path_value = case.get("configPath")
        if not isinstance(config_path_value, str) or not config_path_value:
            raise WorkflowError(f"Each manifest case requires configPath: {manifest_path}")

        case_tag = require_string(case, "caseTag", f"{manifest_path}.cases[]")
        aggregated_case = aggregate_case_diagnostics(Path(config_path_value))
        case_summary = {
            "caseTag": case_tag,
            "label": case.get("label"),
            "value": case.get("value"),
            "configPath": aggregated_case["configPath"],
            "summaryDirectory": aggregated_case["summaryDirectory"],
        }
        case_summaries.append(case_summary)

        for row in aggregated_case["diagnosticSummary"]["rows"]:
            diagnostic_rows.append(dict(case_summary, **row))

        negative_ion_summary = aggregated_case.get("negativeIonSummary")
        if isinstance(negative_ion_summary, dict):
            for row in negative_ion_summary["rows"]:
                negative_ion_rows.append(dict(case_summary, **row))

        for species_tag, grid_power_summary in aggregated_case["gridPowerSummaries"].items():
            for row in grid_power_summary["rows"]:
                grid_power_rows.append(dict(case_summary, speciesTag=species_tag, **row))

    return {
        "metadata": manifest.get("metadata", {}),
        "manifestPath": manifest_path.as_posix(),
        "cases": case_summaries,
        "diagnosticRows": diagnostic_rows,
        "negativeIonRows": negative_ion_rows,
        "gridPowerRows": grid_power_rows,
    }


def write_scan_synthesis(manifest_path: Path, output_path: Path | None = None) -> Path:
    payload = synthesize_scan_results(manifest_path)
    resolved_output = output_path or (manifest_path.parent / "scan_diagnostics_summary.json")
    write_json(resolved_output, payload)
    return resolved_output