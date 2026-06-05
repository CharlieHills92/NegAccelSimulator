"""Boundary-condition builders for workflow materialization."""

from __future__ import annotations

import copy
from typing import Any

from ..common import WorkflowError, require_number, require_object


def build_electrodes(boundary_conditions: dict[str, Any]) -> list[dict[str, Any]]:
    voltages = require_object(boundary_conditions, "gridVoltagesVolts", "boundaryConditions")
    extraction_grid = require_number(voltages, "extractionGrid", "boundaryConditions.gridVoltagesVolts")
    accelerator_stages = voltages.get("acceleratorStages")
    if not isinstance(accelerator_stages, list) or not accelerator_stages:
        raise WorkflowError("boundaryConditions.gridVoltagesVolts requires a non-empty acceleratorStages array")

    electrodes = [
        {
            "name": "EG",
            "role": "extraction_grid",
            "stage": 0,
            "voltageVolts": extraction_grid,
        }
    ]
    for stage_index, voltage in enumerate(accelerator_stages, start=1):
        if not isinstance(voltage, (int, float)):
            raise WorkflowError("boundaryConditions.gridVoltagesVolts.acceleratorStages must contain only numbers")
        electrodes.append(
            {
                "name": f"AG{stage_index}",
                "role": "accelerator_grid",
                "stage": stage_index,
                "voltageVolts": float(voltage),
            }
        )
    return electrodes


def build_boundary_conditions(boundary_conditions: dict[str, Any]) -> dict[str, Any]:
    runtime_boundary_conditions = {
        "electrodes": build_electrodes(boundary_conditions),
    }

    periodic_boundaries = boundary_conditions.get("periodicBoundaries")
    if periodic_boundaries is not None:
        if not isinstance(periodic_boundaries, dict):
            raise WorkflowError("boundaryConditions.periodicBoundaries must be an object when provided")
        runtime_boundary_conditions["periodicBoundaries"] = copy.deepcopy(periodic_boundaries)

    return runtime_boundary_conditions
