"""Particle source builders for workflow materialization."""

from __future__ import annotations

from typing import Any

from ..common import require_number


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
