"""Diagnostics builders for workflow materialization."""

from __future__ import annotations

import copy
from typing import Any

from ..common import WorkflowError
from .geometry_validation import derive_diagnostic_planes, validate_diagnostic_planes


def _require_object(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise WorkflowError(f"{context} must be an object")
    return value


def _require_boolean(value: Any, context: str) -> bool:
    if not isinstance(value, bool):
        raise WorkflowError(f"{context} must be a boolean")
    return value


def _require_number(value: Any, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise WorkflowError(f"{context} must be numeric")
    return float(value)


def _require_positive_number(value: Any, context: str) -> float:
    number = _require_number(value, context)
    if number <= 0.0:
        raise WorkflowError(f"{context} must be greater than zero")
    return number


def _require_nonnegative_integer(value: Any, context: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise WorkflowError(f"{context} must be an integer")
    if value < 0:
        raise WorkflowError(f"{context} must be >= 0")
    return value


def _require_string(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise WorkflowError(f"{context} must be a non-empty string")
    return value.strip()


def _require_number_list(value: Any, context: str, *, min_items: int) -> list[float]:
    if not isinstance(value, list):
        raise WorkflowError(f"{context} must be an array")
    values = [_require_number(item, f"{context}[{index}]") for index, item in enumerate(value)]
    if len(values) < min_items:
        raise WorkflowError(f"{context} must contain at least {min_items} values")
    return values


def build_diagnostics(authoring_diagnostics: dict[str, Any], geometry: dict[str, Any]) -> dict[str, Any]:
    diagnostics = copy.deepcopy(_require_object(authoring_diagnostics, "diagnostics"))
    derived_planes = derive_diagnostic_planes(geometry)

    planes = _require_object(diagnostics.get("planes"), "diagnostics.planes")
    if "sampleZPositionsMeters" in planes:
        planes["sampleZPositionsMeters"] = _require_number_list(
            planes.get("sampleZPositionsMeters"),
            "diagnostics.planes.sampleZPositionsMeters",
            min_items=1,
        )
    else:
        planes["sampleZPositionsMeters"] = list(derived_planes["sampleZPositionsMeters"])
    if "summaryZPositionMeters" in planes:
        planes["summaryZPositionMeters"] = _require_number(
            planes.get("summaryZPositionMeters"),
            "diagnostics.planes.summaryZPositionMeters",
        )
    else:
        planes["summaryZPositionMeters"] = float(derived_planes["summaryZPositionMeters"])
    if "emitterExportZPositionMeters" in planes:
        planes["emitterExportZPositionMeters"] = _require_number(
            planes.get("emitterExportZPositionMeters"),
            "diagnostics.planes.emitterExportZPositionMeters",
        )
    else:
        planes["emitterExportZPositionMeters"] = float(derived_planes["emitterExportZPositionMeters"])

    species = _require_object(diagnostics.get("species"), "diagnostics.species")
    species["writePerSpeciesDiagnostics"] = _require_boolean(
        species.get("writePerSpeciesDiagnostics"),
        "diagnostics.species.writePerSpeciesDiagnostics",
    )
    species["writePerSpeciesGridPower"] = _require_boolean(
        species.get("writePerSpeciesGridPower"),
        "diagnostics.species.writePerSpeciesGridPower",
    )
    species["writePerSpeciesPlots"] = _require_boolean(
        species.get("writePerSpeciesPlots"),
        "diagnostics.species.writePerSpeciesPlots",
    )
    species["writeNegativeIonSummary"] = _require_boolean(
        species.get("writeNegativeIonSummary"),
        "diagnostics.species.writeNegativeIonSummary",
    )

    grid_power = _require_object(diagnostics.get("gridPower"), "diagnostics.gridPower")
    raw_ranges = grid_power.get("ranges")
    if not isinstance(raw_ranges, list) or not raw_ranges:
        raise WorkflowError("diagnostics.gridPower.ranges must be a non-empty array")

    normalized_ranges: list[dict[str, Any]] = []
    seen_ids: set[int] = set()
    for index, raw_range in enumerate(raw_ranges):
        context = f"diagnostics.gridPower.ranges[{index}]"
        range_definition = _require_object(raw_range, context)
        range_id = _require_nonnegative_integer(range_definition.get("id"), context + ".id")
        if range_id in seen_ids:
            raise WorkflowError(f"{context}.id duplicates '{range_id}'")
        seen_ids.add(range_id)

        normalized_ranges.append(
            {
                "id": range_id,
                "includeInTotal": _require_boolean(
                    range_definition.get("includeInTotal"),
                    context + ".includeInTotal",
                ),
            }
        )
    grid_power["ranges"] = normalized_ranges

    summary = _require_object(diagnostics.get("summary"), "diagnostics.summary")
    summary["apertureRadiusMeters"] = _require_positive_number(
        summary.get("apertureRadiusMeters"),
        "diagnostics.summary.apertureRadiusMeters",
    )
    if "transmissionPlaneZPositionMeters" in summary:
        summary["transmissionPlaneZPositionMeters"] = _require_number(
            summary.get("transmissionPlaneZPositionMeters"),
            "diagnostics.summary.transmissionPlaneZPositionMeters",
        )
    else:
        summary["transmissionPlaneZPositionMeters"] = float(derived_planes["transmissionPlaneZPositionMeters"])

    plots = _require_object(diagnostics.get("plots"), "diagnostics.plots")
    meniscus = _require_object(plots.get("meniscus"), "diagnostics.plots.meniscus")
    meniscus["enabled"] = _require_boolean(
        meniscus.get("enabled"),
        "diagnostics.plots.meniscus.enabled",
    )
    meniscus["zMinMeters"] = _require_number(
        meniscus.get("zMinMeters"),
        "diagnostics.plots.meniscus.zMinMeters",
    )
    meniscus["zMaxMeters"] = _require_number(
        meniscus.get("zMaxMeters"),
        "diagnostics.plots.meniscus.zMaxMeters",
    )
    meniscus["transverseMinMeters"] = _require_number(
        meniscus.get("transverseMinMeters"),
        "diagnostics.plots.meniscus.transverseMinMeters",
    )
    meniscus["transverseMaxMeters"] = _require_number(
        meniscus.get("transverseMaxMeters"),
        "diagnostics.plots.meniscus.transverseMaxMeters",
    )
    if meniscus["zMaxMeters"] <= meniscus["zMinMeters"]:
        raise WorkflowError("diagnostics.plots.meniscus.zMaxMeters must be greater than zMinMeters")
    if meniscus["transverseMaxMeters"] <= meniscus["transverseMinMeters"]:
        raise WorkflowError(
            "diagnostics.plots.meniscus.transverseMaxMeters must be greater than transverseMinMeters"
        )

    validate_diagnostic_planes(geometry, diagnostics)

    return diagnostics