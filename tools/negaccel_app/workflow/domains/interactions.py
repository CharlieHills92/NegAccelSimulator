"""Gas and surface-interaction builders for workflow materialization."""

from __future__ import annotations

from typing import Any

from ...particles import BUILTIN_PARTICLE_KINDS, get_particle_kind_definition, map_particle_kind_to_family
from ..common import WorkflowError, require_object, require_string


_LEGACY_REACTION_DEFINITIONS = {
    "negative_ion_single_stripping": {
        "name": "Negative ion single stripping",
        "projectileTemplates": ["H-"],
        "projectileFate": "consume",
        "products": [
            {"particleKind": "H0", "count": 1, "speedClass": "fast"},
            {"particleKind": "e-", "count": 1, "speedClass": "fast"},
        ],
    },
    "negative_ion_double_stripping": {
        "name": "Negative ion double stripping",
        "projectileTemplates": ["H-"],
        "projectileFate": "consume",
        "products": [
            {"particleKind": "H+", "count": 1, "speedClass": "fast"},
            {"particleKind": "e-", "count": 2, "speedClass": "fast"},
        ],
    },
    "background_gas_ionization": {
        "name": "Background gas ionization",
        "projectileTemplates": ["H-", "H0"],
        "projectileFate": "survive",
        "products": [
            {"particleKind": "H2+", "count": 1, "speedClass": "slow"},
            {"particleKind": "e-", "count": 1, "speedClass": "slow"},
        ],
    },
    "neutral_projectile_ionization": {
        "name": "Neutral projectile ionization",
        "projectileTemplates": ["H0"],
        "projectileFate": "consume",
        "products": [
            {"particleKind": "H+", "count": 1, "speedClass": "fast"},
            {"particleKind": "e-", "count": 1, "speedClass": "fast"},
        ],
    },
    "positive_ion_charge_exchange": {
        "name": "Positive ion charge exchange",
        "projectileTemplates": ["H+"],
        "projectileFate": "consume",
        "products": [
            {"particleKind": "H2+", "count": 1, "speedClass": "slow"},
            {"particleKind": "H0", "count": 1, "speedClass": "fast"},
        ],
    },
}
_SUPPORTED_PROJECTILE_FATES = {"consume", "survive"}
_SUPPORTED_PRODUCT_SPEED_CLASSES = {"fast", "slow"}


def _infer_particle_family_from_gas_species(species: str) -> str:
    return "D" if species.strip().upper().startswith("D") else "H"


def _normalize_process_products(products: Any, context: str) -> list[dict[str, Any]]:
    if products is None:
        return []
    if not isinstance(products, list):
        raise WorkflowError(f"{context}.products must be an array when provided")

    normalized_products: list[dict[str, Any]] = []
    for product_index, product in enumerate(products):
        product_context = f"{context}.products[{product_index}]"
        if not isinstance(product, dict):
            raise WorkflowError(f"{product_context} must be an object")
        particle_kind = require_string(product, "particleKind", product_context)
        if particle_kind not in BUILTIN_PARTICLE_KINDS:
            raise WorkflowError(f"{product_context}.particleKind is not supported: {particle_kind}")

        count = product.get("count", 1)
        if not isinstance(count, int) or isinstance(count, bool) or count < 1:
            raise WorkflowError(f"{product_context}.count must be an integer >= 1")

        speed_class = require_string(product, "speedClass", product_context).strip().lower()
        if speed_class not in _SUPPORTED_PRODUCT_SPEED_CLASSES:
            raise WorkflowError(f"{product_context}.speedClass must be one of: fast, slow")

        definition = get_particle_kind_definition(particle_kind)
        normalized_products.append(
            {
                "particleKind": particle_kind,
                "count": count,
                "speedClass": speed_class,
                "chargeState": float(definition.charge_state),
                "massU": float(definition.mass_u),
            }
        )

    return normalized_products


def _process_id_suffix(particle_kind: str) -> str:
    return (
        particle_kind.lower()
        .replace("+", "plus")
        .replace("-", "minus")
    )


def _normalize_legacy_reaction_processes(
    reaction: dict[str, Any],
    *,
    particle_family: str,
    context: str,
) -> list[dict[str, Any]]:
    reaction_id = require_string(reaction, "reactionId", context)
    legacy = _LEGACY_REACTION_DEFINITIONS.get(reaction_id)
    if legacy is None:
        raise WorkflowError(
            f"{context}.reactionId must be migrated to processId/projectileKind before use: {reaction_id}"
        )

    normalized_processes: list[dict[str, Any]] = []
    projectile_templates = list(legacy["projectileTemplates"])
    multiple_projectiles = len(projectile_templates) > 1
    base_name = require_string(reaction, "name", context) if "name" in reaction else str(legacy["name"])

    for projectile_template in projectile_templates:
        projectile_kind = map_particle_kind_to_family(str(projectile_template), particle_family)
        mapped_products = []
        for product in legacy["products"]:
            mapped_kind = map_particle_kind_to_family(str(product["particleKind"]), particle_family)
            mapped_products.append(
                {
                    "particleKind": mapped_kind,
                    "count": int(product["count"]),
                    "speedClass": str(product["speedClass"]),
                }
            )

        process_id = reaction_id
        process_name = base_name
        if multiple_projectiles:
            process_id = f"{reaction_id}_{_process_id_suffix(projectile_kind)}"
            process_name = f"{base_name} ({projectile_kind})"

        normalized = {
            "processId": process_id,
            "name": process_name,
            "projectileKind": projectile_kind,
            "sourcePath": require_string(reaction, "sourcePath", context),
            "fitDegree": reaction.get("fitDegree"),
            "coefficients": reaction.get("coefficients"),
            "scaleEnergyByIonMass": bool(reaction.get("scaleEnergyByIonMass", True)),
            "projectileFate": str(legacy["projectileFate"]),
            "products": _normalize_process_products(mapped_products, context),
        }
        if "minimumEnergyEV" in reaction:
            normalized["minimumEnergyEV"] = float(reaction["minimumEnergyEV"])
        if "maximumEnergyEV" in reaction:
            normalized["maximumEnergyEV"] = float(reaction["maximumEnergyEV"])
        normalized_processes.append(normalized)

    return normalized_processes


def _normalize_reaction_processes(
    reaction: dict[str, Any],
    *,
    particle_family: str,
    context: str,
) -> list[dict[str, Any]]:
    if "processId" not in reaction:
        return _normalize_legacy_reaction_processes(reaction, particle_family=particle_family, context=context)

    process_id = require_string(reaction, "processId", context)
    projectile_kind = require_string(reaction, "projectileKind", context)
    if projectile_kind not in BUILTIN_PARTICLE_KINDS:
        raise WorkflowError(f"{context}.projectileKind is not supported: {projectile_kind}")

    fit_degree = reaction.get("fitDegree")
    if not isinstance(fit_degree, int) or isinstance(fit_degree, bool) or fit_degree < 0 or fit_degree > 6:
        raise WorkflowError(f"{context}.fitDegree must be an integer between 0 and 6")

    coefficients = reaction.get("coefficients")
    if not isinstance(coefficients, list) or len(coefficients) != fit_degree + 1:
        raise WorkflowError(f"{context}.coefficients must contain exactly fitDegree + 1 numeric values")

    normalized_coefficients: list[float] = []
    for coefficient in coefficients:
        if not isinstance(coefficient, (int, float)) or isinstance(coefficient, bool):
            raise WorkflowError(f"{context}.coefficients must contain only numeric values")
        normalized_coefficients.append(float(coefficient))

    projectile_fate = require_string(reaction, "projectileFate", context).strip().lower()
    if projectile_fate not in _SUPPORTED_PROJECTILE_FATES:
        raise WorkflowError(f"{context}.projectileFate must be one of: consume, survive")

    normalized: dict[str, Any] = {
        "processId": process_id,
        "name": require_string(reaction, "name", context) if "name" in reaction else process_id,
        "projectileKind": projectile_kind,
        "sourcePath": require_string(reaction, "sourcePath", context),
        "fitDegree": fit_degree,
        "coefficients": normalized_coefficients,
        "scaleEnergyByIonMass": bool(reaction.get("scaleEnergyByIonMass", True)),
        "projectileFate": projectile_fate,
        "products": _normalize_process_products(reaction.get("products", []), context),
    }
    if "minimumEnergyEV" in reaction:
        normalized["minimumEnergyEV"] = float(reaction["minimumEnergyEV"])
    if "maximumEnergyEV" in reaction:
        normalized["maximumEnergyEV"] = float(reaction["maximumEnergyEV"])
    if "maximumEnergyEV" in normalized and "minimumEnergyEV" in normalized:
        if normalized["maximumEnergyEV"] < normalized["minimumEnergyEV"]:
            raise WorkflowError(f"{context}.maximumEnergyEV must be >= minimumEnergyEV when provided")
    return [normalized]


def _build_cross_section_reactions(authoring: dict[str, Any], *, stripping_enabled: bool) -> list[dict[str, Any]]:
    reactions = authoring.get("reactions", [])
    if not isinstance(reactions, list):
        raise WorkflowError("gasInteractions.reactions must be an array when provided")

    particle_family = _infer_particle_family_from_gas_species(
        str(authoring.get("densityProfile", {}).get("species") or "H2")
    )
    runtime_reactions: list[dict[str, Any]] = []
    seen_process_ids: set[str] = set()
    for index, reaction in enumerate(reactions):
        context = f"gasInteractions.reactions[{index}]"
        if not isinstance(reaction, dict):
            raise WorkflowError(f"{context} must be an object")

        normalized_processes = _normalize_reaction_processes(
            reaction,
            particle_family=particle_family,
            context=context,
        )
        for normalized_process in normalized_processes:
            process_id = str(normalized_process["processId"])
            if process_id in seen_process_ids:
                raise WorkflowError(f"Duplicate gasInteractions.reactions processId: {process_id}")
            seen_process_ids.add(process_id)
            runtime_reactions.append(normalized_process)

    if stripping_enabled:
        if not runtime_reactions:
            raise WorkflowError("gasInteractions.reactions must contain at least one process when stripping is enabled")

    return runtime_reactions


def build_gas_density(authoring: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    density_profile = require_object(authoring, "densityProfile", "gasInteractions")
    profile_name = require_string(density_profile, "name", "gasInteractions.densityProfile")

    profile: dict[str, Any] = {
        "name": profile_name,
        "species": require_string(density_profile, "species", "gasInteractions.densityProfile"),
        "source": {
            "mode": "file",
            "path": require_string(density_profile, "path", "gasInteractions.densityProfile"),
            "format": density_profile.get("format", "z-density-columns"),
        },
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

    reactions_runtime = _build_cross_section_reactions(authoring, stripping_enabled=enabled)

    return {
        "defaultProfile": profile_name,
        "profiles": [profile],
    }, stripping_runtime, reactions_runtime


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
