"""Boundary-condition builders for workflow materialization."""

from __future__ import annotations

import ast
import copy
from typing import Any

from ..common import WorkflowError, evaluate_numeric_expression, normalize_symbol_name


_DOMAIN_BOUNDARIES = [
    {"boundaryId": 1, "name": "x-min", "conditionType": "neumann", "value": 0.0},
    {"boundaryId": 2, "name": "x-max", "conditionType": "neumann", "value": 0.0},
    {"boundaryId": 3, "name": "y-min", "conditionType": "neumann", "value": 0.0},
    {"boundaryId": 4, "name": "y-max", "conditionType": "neumann", "value": 0.0},
    {"boundaryId": 5, "name": "z-min", "conditionType": "dirichlet", "value": 0.0},
    {"boundaryId": 6, "name": "z-max", "conditionType": "neumann", "value": 0.0},
]

_LEGACY_STAGE_NAMES = [
    {"GG", "G1", "AG1"},
    {"REP", "G2", "AG2"},
    {"G3", "AG3"},
    {"G4", "AG4"},
    {"G5", "AG5"},
]


def _inventory_from_geometry(authoring_geometry: dict[str, Any]) -> list[dict[str, Any]]:
    inventory = [copy.deepcopy(boundary) for boundary in _DOMAIN_BOUNDARIES]
    solids = authoring_geometry.get("solids")
    if isinstance(solids, list) and solids:
        used_ids: set[int] = set()
        for index, solid in enumerate(solids):
            if not isinstance(solid, dict):
                continue
            boundary_id = solid.get("boundaryId", 7 + index)
            if not isinstance(boundary_id, int) or boundary_id < 7 or boundary_id in used_ids:
                continue
            used_ids.add(boundary_id)
            default_value = solid.get("voltageVolts", 0.0)
            inventory.append(
                {
                    "boundaryId": boundary_id,
                    "name": str(solid.get("name") or f"Solid {boundary_id}"),
                    "conditionType": "dirichlet",
                    "value": float(default_value) if isinstance(default_value, (int, float)) else 0.0,
                }
            )
    return inventory


def _validate_boundaries(boundaries: Any) -> list[dict[str, Any]]:
    if not isinstance(boundaries, list) or not boundaries:
        raise WorkflowError("boundaryConditions.boundaries must be a non-empty array")

    validated: list[dict[str, Any]] = []
    seen_ids: set[int] = set()
    for index, boundary in enumerate(boundaries):
        context = f"boundaryConditions.boundaries[{index}]"
        if not isinstance(boundary, dict):
            raise WorkflowError(f"{context} must be an object")

        boundary_id = boundary.get("boundaryId")
        if not isinstance(boundary_id, int) or boundary_id < 1:
            raise WorkflowError(f"{context}.boundaryId must be an integer >= 1")
        if boundary_id in seen_ids:
            raise WorkflowError(f"Duplicate boundaryConditions.boundaries boundaryId: {boundary_id}")
        seen_ids.add(boundary_id)

        condition_type = str(boundary.get("conditionType", "dirichlet")).strip().lower()
        if condition_type not in {"dirichlet", "neumann"}:
            raise WorkflowError(f"{context}.conditionType must be dirichlet or neumann")

        raw_value = boundary.get("value")
        if isinstance(raw_value, bool) or not isinstance(raw_value, (int, float, str)):
            raise WorkflowError(f"{context}.value must be numeric or an expression string")
        if isinstance(raw_value, str) and not raw_value.strip():
            raise WorkflowError(f"{context}.value must not be empty")

        validated_boundary = {
            "boundaryId": boundary_id,
            "conditionType": condition_type,
            "value": float(raw_value) if isinstance(raw_value, (int, float)) else raw_value.strip(),
        }
        name = boundary.get("name")
        if isinstance(name, str) and name.strip():
            validated_boundary["name"] = name.strip()
        validated.append(validated_boundary)

    return resolve_boundary_value_expressions(validated)


def _expression_symbol_names(expression: str, context: str) -> set[str]:
    try:
        parsed = ast.parse(expression, mode="eval")
    except SyntaxError as exc:
        raise WorkflowError(f"Invalid {context}: {expression}") from exc

    symbols: set[str] = set()
    for node in ast.walk(parsed):
        if isinstance(node, ast.Name):
            normalized_name = normalize_symbol_name(node.id)
            if normalized_name:
                symbols.add(normalized_name)
    return symbols


def resolve_boundary_value_expressions(
    boundaries: list[dict[str, Any]],
    context_prefix: str = "boundaryConditions.boundaries",
) -> list[dict[str, Any]]:
    named_rows: dict[str, list[int]] = {}
    for index, boundary in enumerate(boundaries):
        name = boundary.get("name")
        if not isinstance(name, str) or not name.strip():
            continue
        normalized_name = normalize_symbol_name(name)
        if not normalized_name:
            continue
        named_rows.setdefault(normalized_name, []).append(index)

    resolved: dict[int, float] = {}
    resolving: set[int] = set()

    def _resolve_index(index: int) -> float:
        if index in resolved:
            return resolved[index]
        if index in resolving:
            boundary_id = boundaries[index]["boundaryId"]
            raise WorkflowError(
                f"{context_prefix}[{index}].value has a cyclic expression dependency at boundaryId {boundary_id}"
            )

        boundary = boundaries[index]
        context = f"{context_prefix}[{index}].value"
        raw_value = boundary["value"]
        if isinstance(raw_value, (int, float)):
            value = float(raw_value)
            resolved[index] = value
            boundary["value"] = value
            return value

        resolving.add(index)
        symbols: dict[str, float] = {}
        for symbol_name in _expression_symbol_names(raw_value, context):
            matches = named_rows.get(symbol_name, [])
            if not matches:
                raise WorkflowError(f"Unknown boundary name '{symbol_name}' in {context}")
            if len(matches) != 1:
                raise WorkflowError(f"Ambiguous boundary name '{symbol_name}' in {context}")
            symbols[symbol_name] = _resolve_index(matches[0])

        value = evaluate_numeric_expression(raw_value, symbols, context)
        resolving.remove(index)
        resolved[index] = value
        boundary["value"] = value
        return value

    for index in range(len(boundaries)):
        _resolve_index(index)

    return boundaries


def _legacy_boundaries(
    boundary_conditions: dict[str, Any],
    authoring_geometry: dict[str, Any],
) -> list[dict[str, Any]]:
    boundaries = {
        item["boundaryId"]: copy.deepcopy(item)
        for item in _inventory_from_geometry(authoring_geometry)
    }

    voltages = boundary_conditions.get("gridVoltagesVolts")
    if not isinstance(voltages, dict):
        return [boundaries[boundary_id] for boundary_id in sorted(boundaries)]

    extraction_grid = voltages.get("extractionGrid")
    if isinstance(extraction_grid, (int, float)):
        for boundary in boundaries.values():
            if str(boundary.get("name", "")).upper() == "EG":
                boundary["value"] = float(extraction_grid)
                break

    accelerator_stages = voltages.get("acceleratorStages")
    if accelerator_stages is not None and (not isinstance(accelerator_stages, list) or not accelerator_stages):
        raise WorkflowError("boundaryConditions.gridVoltagesVolts.acceleratorStages must be a non-empty array")

    if isinstance(accelerator_stages, list):
        for stage_index, voltage in enumerate(accelerator_stages[: len(_LEGACY_STAGE_NAMES)], start=1):
            if not isinstance(voltage, (int, float)):
                raise WorkflowError("boundaryConditions.gridVoltagesVolts.acceleratorStages must contain only numbers")
            valid_names = _LEGACY_STAGE_NAMES[stage_index - 1]
            for boundary in boundaries.values():
                if str(boundary.get("name", "")).upper() in valid_names:
                    boundary["value"] = float(voltage)
                    break

    return [boundaries[boundary_id] for boundary_id in sorted(boundaries)]


def build_boundary_conditions(
    boundary_conditions: dict[str, Any],
    authoring_geometry: dict[str, Any],
) -> dict[str, Any]:
    runtime_boundary_conditions = {
        "boundaries": _validate_boundaries(boundary_conditions.get("boundaries"))
        if "boundaries" in boundary_conditions
        else _legacy_boundaries(boundary_conditions, authoring_geometry),
    }

    periodic_boundaries = boundary_conditions.get("periodicBoundaries")
    if periodic_boundaries is not None:
        if not isinstance(periodic_boundaries, dict):
            raise WorkflowError("boundaryConditions.periodicBoundaries must be an object when provided")
        runtime_boundary_conditions["periodicBoundaries"] = copy.deepcopy(periodic_boundaries)

    return runtime_boundary_conditions
