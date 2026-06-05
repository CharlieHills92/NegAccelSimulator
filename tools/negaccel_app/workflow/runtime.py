"""Runtime case materialization for NegAccel workflow tooling."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

from .common import (
    BUILTIN_TEMPLATE_TO_ACCELERATOR,
    WorkflowError,
    apply_cli_overrides,
    ensure_case_metadata,
    load_json,
    normalize_template_name,
    require_object,
    require_string,
    resolve_default_case_output,
    sanitize_case_tag,
    write_json,
)
from .domains.boundary import build_boundary_conditions
from .domains.fields import build_magnetic_field
from .domains.geometry import build_runtime_geometry
from .domains.interactions import build_gas_density, build_surface_interactions
from .domains.outputs import build_outputs
from .domains.setup import build_simulation
from .domains.source import build_particle_sources


def authored_to_runtime_case(authoring_spec: dict[str, Any], case_tag_override: str | None = None) -> dict[str, Any]:
    metadata = require_object(authoring_spec, "metadata", "authoring specification")
    geometry = require_object(authoring_spec, "geometry", "authoring specification")
    boundary_conditions = require_object(authoring_spec, "boundaryConditions", "authoring specification")
    particle_source = require_object(authoring_spec, "particleSource", "authoring specification")
    gas_interactions = require_object(authoring_spec, "gasInteractions", "authoring specification")
    run = require_object(authoring_spec, "run", "authoring specification")

    resolved_case_tag = sanitize_case_tag(case_tag_override or require_string(metadata, "caseTag", "metadata"))
    template = normalize_template_name(require_string(geometry, "template", "geometry"))

    gas_density, stripping = build_gas_density(gas_interactions)
    outputs = build_outputs(authoring_spec.get("outputs"))

    runtime_case: dict[str, Any] = {
        "metadata": {
            "schemaVersion": require_string(metadata, "schemaVersion", "metadata"),
            "caseTag": resolved_case_tag,
        },
        "accelerator": {
            "type": template,
            "legacyIndex": BUILTIN_TEMPLATE_TO_ACCELERATOR[template],
            "useBuiltInDefaults": True,
        },
        "geometry": build_runtime_geometry(geometry, template, outputs),
        "particleSources": build_particle_sources(particle_source),
        "simulation": build_simulation(run),
        "boundaryConditions": build_boundary_conditions(boundary_conditions),
        "externalMagneticField": build_magnetic_field(authoring_spec.get("magneticField", {})),
        "gasDensity": gas_density,
        "physics": {
            "stripping": stripping,
            "surfaceCollisions": build_surface_interactions(authoring_spec.get("surfaceInteractions")),
        },
        "outputs": outputs,
    }

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
