"""Output builders for workflow materialization."""

from __future__ import annotations

import copy
from typing import Any

from ..common import DEFAULT_OUTPUTS, WorkflowError, merge_objects
from .geometry_validation import get_domain_z_bounds


def _require_boolean(value: Any, context: str) -> bool:
    if not isinstance(value, bool):
        raise WorkflowError(f"{context} must be a boolean")
    return value


def _require_integer(value: Any, context: str, *, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise WorkflowError(f"{context} must be an integer")
    if value < minimum:
        raise WorkflowError(f"{context} must be >= {minimum}")
    return value


def _require_plane_list(value: Any, context: str) -> list[float]:
    if not isinstance(value, list):
        raise WorkflowError(f"{context} must be an array")

    planes: list[float] = []
    for index, item in enumerate(value):
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            raise WorkflowError(f"{context}[{index}] must be numeric")
        plane = float(item)
        if plane not in planes:
            planes.append(plane)
    return planes


def _validate_iteration_outputs(outputs: dict[str, Any], geometry: dict[str, Any] | None) -> None:
    iteration = outputs.get("iteration")
    if not isinstance(iteration, dict):
        raise WorkflowError("outputs.iteration must be an object")

    iteration["enabled"] = _require_boolean(iteration.get("enabled"), "outputs.iteration.enabled")
    iteration["everyNIterations"] = _require_integer(
        iteration.get("everyNIterations"),
        "outputs.iteration.everyNIterations",
        minimum=1,
    )
    iteration["exportPlaneDiagnostics"] = _require_boolean(
        iteration.get("exportPlaneDiagnostics"),
        "outputs.iteration.exportPlaneDiagnostics",
    )
    iteration["exportSimulationState"] = _require_boolean(
        iteration.get("exportSimulationState"),
        "outputs.iteration.exportSimulationState",
    )
    iteration["exportTracedParticles"] = _require_boolean(
        iteration.get("exportTracedParticles"),
        "outputs.iteration.exportTracedParticles",
    )

    planes = _require_plane_list(
        iteration.get("planeZPositionsMeters", []),
        "outputs.iteration.planeZPositionsMeters",
    )
    if geometry is not None and planes:
        domain_z_min, domain_z_max = get_domain_z_bounds(geometry)
        for index, plane in enumerate(planes):
            if plane < domain_z_min or plane > domain_z_max:
                raise WorkflowError(
                    f"outputs.iteration.planeZPositionsMeters[{index}]={plane:g} m "
                    f"is outside geometry.domain z range [{domain_z_min:g}, {domain_z_max:g}] m"
                )
    iteration["planeZPositionsMeters"] = planes


def build_outputs(
    authoring_outputs: dict[str, Any] | None,
    geometry: dict[str, Any] | None = None,
) -> dict[str, Any]:
    outputs = copy.deepcopy(DEFAULT_OUTPUTS)
    if authoring_outputs is None:
        _validate_iteration_outputs(outputs, geometry)
        return outputs
    if not isinstance(authoring_outputs, dict):
        raise WorkflowError("outputs must be an object when provided")
    merge_objects(outputs, authoring_outputs)
    _validate_iteration_outputs(outputs, geometry)
    return outputs
