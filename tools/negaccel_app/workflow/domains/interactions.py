"""Gas and surface-interaction builders for workflow materialization."""

from __future__ import annotations

from typing import Any

from ..common import WorkflowError, require_object, require_string


def build_gas_density(authoring: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    density_profile = require_object(authoring, "densityProfile", "gasInteractions")
    profile_name = require_string(density_profile, "name", "gasInteractions.densityProfile")
    use_for = density_profile.get("useFor", ["stripping", "secondaryEmission", "diagnostics"])
    if not isinstance(use_for, list) or not use_for:
        raise WorkflowError("gasInteractions.densityProfile.useFor must be a non-empty array when provided")

    profile: dict[str, Any] = {
        "name": profile_name,
        "species": require_string(density_profile, "species", "gasInteractions.densityProfile"),
        "source": {
            "mode": "file",
            "path": require_string(density_profile, "path", "gasInteractions.densityProfile"),
            "format": density_profile.get("format", "z-density-columns"),
        },
        "useFor": use_for,
    }
    if "pressureScale" in density_profile:
        profile["pressureScale"] = float(density_profile["pressureScale"])

    stripping = authoring.get("stripping", {})
    if not isinstance(stripping, dict):
        raise WorkflowError("gasInteractions.stripping must be an object when provided")

    enabled = bool(stripping.get("enabled", False))
    if not enabled:
        stripping_mode = "disabled"
    elif bool(stripping.get("generateSecondaries", False)):
        stripping_mode = "withSecondaries"
    else:
        stripping_mode = "primaryOnly"

    stripping_runtime: dict[str, Any] = {
        "mode": stripping_mode,
    }
    if "minimumZMeters" in stripping:
        stripping_runtime["minimumZMeters"] = float(stripping["minimumZMeters"])
    if enabled:
        stripping_runtime["densityProfile"] = profile_name

    return {
        "defaultProfile": profile_name,
        "profiles": [profile],
    }, stripping_runtime


def build_surface_interactions(authoring: dict[str, Any] | None) -> dict[str, Any]:
    if not authoring:
        return {
            "enabled": False,
        }

    runtime_surface: dict[str, Any] = {
        "enabled": bool(authoring.get("enabled", False)),
    }
    if "model" in authoring:
        runtime_surface["model"] = authoring["model"]
    if "minimumImpactZMeters" in authoring:
        runtime_surface["minimumImpactZMeters"] = float(authoring["minimumImpactZMeters"])
    if "maximumSecondaryElectrons" in authoring:
        runtime_surface["maximumSecondaryElectrons"] = int(authoring["maximumSecondaryElectrons"])
    if "secondaryElectronEnergyEV" in authoring:
        runtime_surface["secondaryElectronEnergyEV"] = float(authoring["secondaryElectronEnergyEV"])
    if "debug" in authoring:
        runtime_surface["debug"] = bool(authoring["debug"])
    return runtime_surface
