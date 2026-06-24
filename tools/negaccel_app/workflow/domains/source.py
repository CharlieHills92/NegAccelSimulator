"""Particle builders for workflow materialization."""

from __future__ import annotations

from typing import Any

from ...particles import BUILTIN_PARTICLE_KINDS, get_particle_kind_family
from ..common import WorkflowError, require_number, require_object, require_string


SUPPORTED_PLASMA_MODELS = {"nsimp", "shield"}
SUPPORTED_SOURCE_MODELS = {"uniform"}


def build_particle_sources(particle_source: dict[str, Any]) -> dict[str, Any]:
    negative_ion_beam: dict[str, Any] = {
        "species": particle_source.get("species", "H-"),
        "chargeState": float(particle_source.get("chargeState", -1)),
        "massU": require_number(particle_source, "massU", "particleSource"),
        "currentDensityAm2": require_number(particle_source, "currentDensityAm2", "particleSource"),
        "perpendicularTemperatureEV": float(particle_source.get("perpendicularTemperatureEV", 0.0)),
        "axialEnergyEV": require_number(particle_source, "axialEnergyEV", "particleSource"),
        "sourceModel": particle_source.get("sourceModel", "uniform"),
    }

    if "parallelTemperatureEV" in particle_source:
        negative_ion_beam["parallelTemperatureEV"] = float(particle_source["parallelTemperatureEV"])
    if "plasmaPotentialVolts" in particle_source:
        negative_ion_beam["plasmaPotentialVolts"] = float(particle_source["plasmaPotentialVolts"])
    if "electronsModelWeight" in particle_source:
        negative_ion_beam["electronsModelWeight"] = float(particle_source["electronsModelWeight"])

    return {
        "negativeIonBeam": negative_ion_beam,
    }


def build_particle_configuration(particles: dict[str, Any]) -> dict[str, Any]:
    runtime_particle_types, particle_type_lookup = build_particle_types(particles)
    runtime_particle_sources = build_runtime_particle_sources(particles, particle_type_lookup)

    return {
        "particleTypes": runtime_particle_types,
        "particleSources": runtime_particle_sources,
        "plasma": build_particle_plasma(particles),
        "particleCount": sum(source["particleCount"] for source in runtime_particle_sources),
    }


def build_particle_types(particles: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    raw_types = particles.get("types")
    if not isinstance(raw_types, list) or not raw_types:
        raise WorkflowError("particles requires a non-empty array 'types'")

    runtime_particle_types: list[dict[str, Any]] = []
    particle_type_lookup: dict[str, dict[str, Any]] = {}

    for index, raw_type in enumerate(raw_types):
        context = f"particles.types[{index}]"
        if not isinstance(raw_type, dict):
            raise WorkflowError(f"{context} must be an object")

        type_id = require_string(raw_type, "id", context)
        if type_id in particle_type_lookup:
            raise WorkflowError(f"Duplicate particle type id '{type_id}'")

        kind = require_string(raw_type, "kind", context)
        if kind not in BUILTIN_PARTICLE_KINDS:
            supported = ", ".join(sorted(BUILTIN_PARTICLE_KINDS))
            raise WorkflowError(f"Unsupported {context}.kind '{kind}'. Supported kinds: {supported}")

        kind_definition = BUILTIN_PARTICLE_KINDS[kind]
        runtime_particle_type = {
            "id": type_id,
            "kind": kind_definition.kind,
            "name": str(raw_type.get("name") or kind_definition.label),
            "chargeState": kind_definition.charge_state,
            "massU": kind_definition.mass_u,
            "sourceable": kind_definition.sourceable,
        }

        runtime_particle_types.append(runtime_particle_type)
        particle_type_lookup[type_id] = runtime_particle_type

    return runtime_particle_types, particle_type_lookup


def build_runtime_particle_sources(
    particles: dict[str, Any],
    particle_type_lookup: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    raw_sources = particles.get("sources")
    if not isinstance(raw_sources, list) or not raw_sources:
        raise WorkflowError("particles requires a non-empty array 'sources'")

    runtime_particle_sources: list[dict[str, Any]] = []
    source_ids: set[str] = set()
    heavy_source_families: set[str] = set()

    for index, raw_source in enumerate(raw_sources):
        context = f"particles.sources[{index}]"
        if not isinstance(raw_source, dict):
            raise WorkflowError(f"{context} must be an object")

        source_id = require_string(raw_source, "id", context)
        if source_id in source_ids:
            raise WorkflowError(f"Duplicate particle source id '{source_id}'")
        source_ids.add(source_id)

        particle_type_id = require_string(raw_source, "particleTypeId", context)
        if particle_type_id not in particle_type_lookup:
            raise WorkflowError(f"{context}.particleTypeId references unknown particle type '{particle_type_id}'")

        particle_type = particle_type_lookup[particle_type_id]
        if not bool(particle_type.get("sourceable")):
            raise WorkflowError(
                f"{context}.particleTypeId must reference a sourceable kind; got '{particle_type['kind']}'"
            )

        source_model = str(raw_source.get("sourceModel", "uniform") or "uniform")
        if source_model not in SUPPORTED_SOURCE_MODELS:
            supported = ", ".join(sorted(SUPPORTED_SOURCE_MODELS))
            raise WorkflowError(f"Unsupported {context}.sourceModel '{source_model}'. Supported models: {supported}")

        particle_count = _require_positive_integer(raw_source, "particleCount", context)
        current_density = require_number(raw_source, "currentDensityAm2", context)
        axial_energy = require_number(raw_source, "axialEnergyEV", context)
        uniform = require_object(raw_source, "uniform", context)

        width = require_number(uniform, "widthMeters", f"{context}.uniform")
        height = require_number(uniform, "heightMeters", f"{context}.uniform")
        if width <= 0.0:
            raise WorkflowError(f"{context}.uniform.widthMeters must be greater than zero")
        if height <= 0.0:
            raise WorkflowError(f"{context}.uniform.heightMeters must be greater than zero")

        center = _require_number_vector(uniform, "centerMeters", 3, f"{context}.uniform")
        main_direction = _require_number_vector(uniform, "mainDirection", 3, f"{context}.uniform")
        in_plane_reference = _require_number_vector(
            uniform,
            "inPlaneReferenceDirection",
            3,
            f"{context}.uniform",
        )
        _require_non_zero_vector(main_direction, f"{context}.uniform.mainDirection")
        _require_non_zero_vector(in_plane_reference, f"{context}.uniform.inPlaneReferenceDirection")

        runtime_particle_sources.append(
            {
                "id": source_id,
                "name": str(raw_source.get("name") or source_id),
                "particleTypeId": particle_type_id,
                "kind": particle_type["kind"],
                "chargeState": float(particle_type["chargeState"]),
                "massU": float(particle_type["massU"]),
                "particleCount": particle_count,
                "currentDensityAm2": float(current_density),
                "perpendicularTemperatureEV": float(raw_source.get("perpendicularTemperatureEV", 0.0)),
                "parallelTemperatureEV": float(raw_source.get("parallelTemperatureEV", 0.0)),
                "axialEnergyEV": float(axial_energy),
                "sourceModel": source_model,
                "uniform": {
                    "centerMeters": center,
                    "mainDirection": main_direction,
                    "inPlaneReferenceDirection": in_plane_reference,
                    "widthMeters": float(width),
                    "heightMeters": float(height),
                },
            }
        )

        source_family = get_particle_kind_family(str(particle_type["kind"]))
        if source_family is not None:
            heavy_source_families.add(source_family)

    if len(heavy_source_families) > 1:
        supported = ", ".join(sorted(heavy_source_families))
        raise WorkflowError(
            "Mixed heavy-ion source families are not yet supported on the active runtime path: " + supported
        )

    return runtime_particle_sources


def build_particle_plasma(particles: dict[str, Any]) -> dict[str, Any]:
    plasma = require_object(particles, "plasma", "particles")
    model = str(plasma.get("model", "nsimp") or "nsimp")
    if model not in SUPPORTED_PLASMA_MODELS:
        supported = ", ".join(sorted(SUPPORTED_PLASMA_MODELS))
        raise WorkflowError(f"Unsupported particles.plasma.model '{model}'. Supported models: {supported}")

    runtime_plasma = {
        "plasmaModel": model,
        "initialPlasmaMaxZMeters": 7.0e-3,
    }
    if "initialPlasmaMaxZMeters" in plasma:
        runtime_plasma["initialPlasmaMaxZMeters"] = float(
            require_number(plasma, "initialPlasmaMaxZMeters", "particles.plasma")
        )
    if model == "nsimp":
        runtime_plasma["positiveIonTemperatureEV"] = float(
            require_number(plasma, "positiveIonTemperatureEV", "particles.plasma")
        )
        runtime_plasma["plasmaPotentialVolts"] = float(
            require_number(plasma, "plasmaPotentialVolts", "particles.plasma")
        )
    else:
        runtime_plasma["tanhWidthEV"] = float(require_number(plasma, "tanhWidthEV", "particles.plasma"))
        runtime_plasma["meniscusVoltageVolts"] = float(
            require_number(plasma, "meniscusVoltageVolts", "particles.plasma")
        )

    return runtime_plasma


def _require_positive_integer(document: dict[str, Any], key: str, context: str) -> int:
    value = document.get(key)
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise WorkflowError(f"{context} requires positive integer '{key}'")
    return int(value)


def _require_number_vector(document: dict[str, Any], key: str, length: int, context: str) -> list[float]:
    raw_value = document.get(key)
    if not isinstance(raw_value, list) or len(raw_value) != length:
        raise WorkflowError(f"{context}.{key} must be an array of {length} numbers")

    vector: list[float] = []
    for index, item in enumerate(raw_value):
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            raise WorkflowError(f"{context}.{key}[{index}] must be numeric")
        vector.append(float(item))

    return vector


def _require_non_zero_vector(vector: list[float], context: str) -> None:
    if not any(component != 0.0 for component in vector):
        raise WorkflowError(f"{context} must not be the zero vector")
