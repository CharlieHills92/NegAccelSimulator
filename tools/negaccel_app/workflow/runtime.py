"""Runtime case materialization for NegAccel workflow tooling."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

from .common import (
    WorkflowError,
    apply_cli_overrides,
    ensure_case_metadata,
    load_json,
    require_object,
    require_string,
    resolve_default_case_output,
    sanitize_case_tag,
    write_json,
)
from .geometry_files import load_geometry_reference
from .domains.boundary import build_boundary_conditions
from .domains.diagnostics import build_diagnostics
from .domains.fields import build_magnetic_field
from .domains.geometry import build_runtime_geometry
from .domains.geometry_validation import validate_uniform_sources_within_domain
from .domains.interactions import build_gas_density, build_surface_interactions
from .domains.outputs import build_outputs
from .domains.setup import build_simulation
from .domains.source import build_particle_configuration, build_particle_sources


def authored_to_runtime_case(authoring_spec: dict[str, Any], case_tag_override: str | None = None) -> dict[str, Any]:
    metadata = require_object(authoring_spec, "metadata", "authoring specification")
    geometry_spec = require_object(authoring_spec, "geometry", "authoring specification")
    boundary_conditions = require_object(authoring_spec, "boundaryConditions", "authoring specification")
    gas_interactions = require_object(authoring_spec, "gasInteractions", "authoring specification")
    diagnostics_spec = require_object(authoring_spec, "diagnostics", "authoring specification")
    run = require_object(authoring_spec, "run", "authoring specification")
    execution_spec = authoring_spec.get("execution")
    authoring_source_path = authoring_spec.get("__source_path") if isinstance(authoring_spec.get("__source_path"), Path) else None

    runtime_execution: dict[str, str] | None = None
    if execution_spec is not None:
        if not isinstance(execution_spec, dict):
            raise WorkflowError("authoring specification execution must be an object when provided")
        if "serverCommand" in execution_spec:
            server_command = require_string(execution_spec, "serverCommand", "execution").strip()
            if not server_command:
                raise WorkflowError("execution requires non-empty string 'serverCommand'")
            runtime_execution = {
                "serverCommand": server_command,
            }

    resolved_case_tag = sanitize_case_tag(case_tag_override or require_string(metadata, "caseTag", "metadata"))

    geometry_reference = geometry_spec.get("path")
    if isinstance(geometry_reference, str) and geometry_reference.strip():
        geometry_path, geometry = load_geometry_reference(
            geometry_reference,
            authoring_source_path,
        )
        geometry["path"] = geometry_path.as_posix()
    else:
        geometry = geometry_spec

    gas_density, stripping, reactions = build_gas_density(gas_interactions)
    outputs = build_outputs(authoring_spec.get("outputs"))
    runtime_geometry = build_runtime_geometry(geometry, outputs)

    runtime_particle_types: list[dict[str, Any]] | None = None
    particles = authoring_spec.get("particles")
    if isinstance(particles, dict):
        particle_configuration = build_particle_configuration(particles)
        runtime_particle_types = particle_configuration["particleTypes"]
        runtime_particle_sources: Any = particle_configuration["particleSources"]
        simulation = build_simulation(
            run,
            particle_configuration["particleCount"],
            particle_configuration["plasma"],
            runtime_geometry,
        )
    else:
        particle_source = require_object(authoring_spec, "particleSource", "authoring specification")
        runtime_particle_sources = build_particle_sources(particle_source)
        simulation = build_simulation(run, runtime_geometry=runtime_geometry)

    validate_uniform_sources_within_domain(runtime_geometry, runtime_particle_sources)
    diagnostics = build_diagnostics(diagnostics_spec, runtime_geometry)
    outputs = build_outputs(authoring_spec.get("outputs"), runtime_geometry)

    runtime_case: dict[str, Any] = {
        "metadata": {
            "schemaVersion": require_string(metadata, "schemaVersion", "metadata"),
            "caseTag": resolved_case_tag,
        },
        "geometry": runtime_geometry,
        "particleSources": runtime_particle_sources,
        "simulation": simulation,
        "boundaryConditions": build_boundary_conditions(boundary_conditions, geometry),
        "externalMagneticField": build_magnetic_field(authoring_spec.get("magneticField", {}), authoring_source_path),
        "gasDensity": gas_density,
        "physics": {
            "reactions": reactions,
            "stripping": stripping,
            "surfaceCollisions": build_surface_interactions(authoring_spec.get("surfaceInteractions")),
        },
        "diagnostics": diagnostics,
        "outputs": outputs,
    }

    if runtime_particle_types is not None:
        runtime_case["particleTypes"] = runtime_particle_types

    if runtime_execution is not None:
        runtime_case["execution"] = runtime_execution

    for key in ("title", "description", "author"):
        if key in metadata:
            runtime_case["metadata"][key] = metadata[key]

    return runtime_case


def materialize_authored_case(
    authoring_path: Path,
    output_path: Path | None,
    case_tag: str | None,
    cli_overrides: Sequence[str],
) -> Path:
    authoring_spec = load_json(authoring_path)
    if not isinstance(authoring_spec, dict):
        raise WorkflowError("Authoring specification must be a JSON object")

    authoring_spec["__source_path"] = authoring_path

    apply_cli_overrides(authoring_spec, cli_overrides)
    runtime_case = authored_to_runtime_case(authoring_spec, case_tag)
    resolved_output = output_path or resolve_default_case_output(runtime_case["metadata"]["caseTag"])
    ensure_case_metadata(runtime_case, runtime_case["metadata"]["caseTag"], resolved_output.parent)
    write_json(resolved_output, runtime_case)
    return resolved_output


def materialize_case(
    template_path: Path,
    output_path: Path | None,
    case_tag: str | None,
    cli_overrides: Sequence[str],
) -> Path:
    case_config = load_json(template_path)
    if not isinstance(case_config, dict):
        raise WorkflowError("Case template must be a JSON object")

    apply_cli_overrides(case_config, cli_overrides)

    resolved_case_tag = case_tag or case_config.get("metadata", {}).get("caseTag")
    if not isinstance(resolved_case_tag, str) or not resolved_case_tag:
        raise WorkflowError("Resolved case tag is missing. Provide metadata.caseTag or --case-tag")
    resolved_case_tag = sanitize_case_tag(resolved_case_tag)

    resolved_output = output_path or resolve_default_case_output(resolved_case_tag)
    ensure_case_metadata(case_config, resolved_case_tag, resolved_output.parent)
    write_json(resolved_output, case_config)
    return resolved_output
