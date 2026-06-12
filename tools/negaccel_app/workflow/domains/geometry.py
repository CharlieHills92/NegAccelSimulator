"""Geometry builders for workflow materialization."""

from __future__ import annotations

import copy
from typing import Any

from ..common import WorkflowError, require_number


_LEGACY_SOLID_KINDS = {"grid", "aperture", "wall", "plasma", "custom"}


def _normalize_solid_kind(kind: Any, context: str) -> str:
    if not isinstance(kind, str) or not kind.strip():
        raise WorkflowError(f"{context}.kind must be a non-empty string")

    normalized_kind = kind.strip()
    if normalized_kind in {"solid", "diagnosticPlane"}:
        return normalized_kind
    if normalized_kind in _LEGACY_SOLID_KINDS:
        return "solid"
    raise WorkflowError(f"{context}.kind must be either solid or diagnosticPlane")


def _normalize_profile_arrays(
    z_profile: Any,
    r_profile: Any,
    rounding_radii: Any,
    context: str,
) -> tuple[list[float], list[float], list[float]]:
    if not isinstance(z_profile, list) or len(z_profile) < 2:
        raise WorkflowError(f"{context}.zProfileMeters must contain at least two entries")
    if not isinstance(r_profile, list) or len(r_profile) < 2:
        raise WorkflowError(f"{context}.rProfileMeters must contain at least two entries")
    if len(z_profile) != len(r_profile):
        raise WorkflowError(f"{context}.zProfileMeters and {context}.rProfileMeters must have the same length")

    normalized_z: list[float] = []
    normalized_r: list[float] = []
    previous_z: float | None = None
    has_positive_span = False
    for point_index, (z_value, r_value) in enumerate(zip(z_profile, r_profile), start=1):
        if not isinstance(z_value, (int, float)) or not isinstance(r_value, (int, float)):
            raise WorkflowError(f"{context}: profile point {point_index} must be numeric")
        z_number = float(z_value)
        r_number = float(r_value)
        if previous_z is not None:
            if z_number < previous_z:
                raise WorkflowError(f"{context}.zProfileMeters must be non-decreasing")
            if z_number > previous_z:
                has_positive_span = True
        if r_number < 0.0:
            raise WorkflowError(f"{context}.rProfileMeters must be non-negative")
        normalized_z.append(z_number)
        normalized_r.append(r_number)
        previous_z = z_number

    if not has_positive_span:
        raise WorkflowError(f"{context}.zProfileMeters must span a non-zero z range")

    if rounding_radii is None:
        return normalized_z, normalized_r, [0.0] * len(normalized_z)
    if not isinstance(rounding_radii, list):
        raise WorkflowError(f"{context}.roundingRadiiMeters must be an array when provided")
    if len(rounding_radii) != len(normalized_z):
        raise WorkflowError(
            f"{context}.roundingRadiiMeters must have the same length as {context}.zProfileMeters"
        )

    normalized_rounding: list[float] = []
    for point_index, radius in enumerate(rounding_radii, start=1):
        if not isinstance(radius, (int, float)):
            raise WorkflowError(f"{context}.roundingRadiiMeters[{point_index - 1}] must be numeric")
        radius_value = float(radius)
        normalized_rounding.append(radius_value)

    return normalized_z, normalized_r, normalized_rounding


def _normalize_generated_solids(authoring_geometry: dict[str, Any]) -> list[dict[str, Any]]:
    solids = authoring_geometry.get("solids")
    if solids is None:
        return []
    if not isinstance(solids, list) or not solids:
        raise WorkflowError("geometry.solids must be a non-empty array when provided")

    normalized_solids: list[dict[str, Any]] = []
    for index, solid in enumerate(solids):
        context = f"geometry.solids[{index}]"
        if not isinstance(solid, dict):
            raise WorkflowError(f"{context} must be an object")
        if not isinstance(solid.get("name"), str) or not solid["name"].strip():
            raise WorkflowError(f"{context}.name must be a non-empty string")

        normalized_solid = copy.deepcopy(solid)
        normalized_solid["kind"] = _normalize_solid_kind(normalized_solid.get("kind"), context)

        z_values, r_values, rounding_radii = _normalize_profile_arrays(
            normalized_solid.get("zProfileMeters"),
            normalized_solid.get("rProfileMeters"),
            normalized_solid.get("roundingRadiiMeters"),
            context,
        )
        normalized_solid["zProfileMeters"] = z_values
        normalized_solid["rProfileMeters"] = r_values
        if any(radius != 0.0 for radius in rounding_radii):
            normalized_solid["roundingRadiiMeters"] = rounding_radii
        else:
            normalized_solid.pop("roundingRadiiMeters", None)

        aperture_pattern = normalized_solid.get("aperturePattern")
        if aperture_pattern is not None and not isinstance(aperture_pattern, dict):
            raise WorkflowError(f"{context}.aperturePattern must be an object when provided")

        normalized_solids.append(normalized_solid)

    return normalized_solids


def build_runtime_geometry(authoring_geometry: dict[str, Any], outputs: dict[str, Any]) -> dict[str, Any]:
    mesh_size = require_number(authoring_geometry, "meshSizeMeters", "geometry")

    domain = authoring_geometry.get("domain")
    if not isinstance(domain, dict):
        raise WorkflowError("geometry.domain must be an object")
    runtime_domain = {
        "xSizeMeters": require_number(domain, "xSizeMeters", "geometry.domain"),
        "ySizeMeters": require_number(domain, "ySizeMeters", "geometry.domain"),
        "zSizeMeters": require_number(domain, "zSizeMeters", "geometry.domain"),
    }
    if "zStartMeters" in domain:
        runtime_domain["zStartMeters"] = require_number(domain, "zStartMeters", "geometry.domain")

    geometry_export_vtk = bool(
        authoring_geometry.get("exportGeometryVtk", outputs.get("vtk", {}).get("exportGeometry", True))
    )

    structured_solids = _normalize_generated_solids(authoring_geometry)
    if not structured_solids:
        raise WorkflowError("geometry.solids must be a non-empty array")

    source: dict[str, Any] = {
        "mode": "generated-data",
        "path": "embedded://geometry.solids",
        "format": "json",
    }

    runtime_geometry = {
        "source": source,
        "mesh": {
            "sizeMeters": mesh_size,
            "exportGeometryVtk": geometry_export_vtk,
        },
        "domain": runtime_domain,
    }

    if "gaps" in authoring_geometry:
        if not isinstance(authoring_geometry["gaps"], dict):
            raise WorkflowError("geometry.gaps must be an object when provided")
        if "accelerationGapMeters" in authoring_geometry["gaps"]:
            raise WorkflowError(
                "geometry.gaps.accelerationGapMeters is no longer supported in generated geometry"
            )
        if "extractionGapMeters" in authoring_geometry["gaps"]:
            runtime_geometry["gaps"] = {
                "extractionGapMeters": copy.deepcopy(authoring_geometry["gaps"]["extractionGapMeters"])
            }

    runtime_geometry["solids"] = structured_solids

    return runtime_geometry
