"""Built-in particle-kind metadata shared by the GUI and workflow."""

from __future__ import annotations

from dataclasses import dataclass


PARTICLE_FAMILIES = ("H", "D")
ELECTRON_PARTICLE_KIND = "e-"


@dataclass(frozen=True)
class ParticleKindDefinition:
    kind: str
    label: str
    charge_state: float
    mass_u: float
    sourceable: bool
    type_id: str
    family: str | None


_BUILTIN_PARTICLE_KIND_DEFINITIONS = (
    ParticleKindDefinition("H-", "H-", -1.0, 1.0, True, "hminus", "H"),
    ParticleKindDefinition("H0", "H0", 0.0, 1.0, True, "h0", "H"),
    ParticleKindDefinition("H+", "H+", 1.0, 1.0, True, "hplus", "H"),
    ParticleKindDefinition("H2+", "H2+", 1.0, 2.0, True, "h2plus", "H"),
    ParticleKindDefinition("H20", "H20", 0.0, 2.0, True, "h20", "H"),
    ParticleKindDefinition("H3+", "H3+", 1.0, 3.0, True, "h3plus", "H"),
    ParticleKindDefinition("D-", "D-", -1.0, 2.0, True, "dminus", "D"),
    ParticleKindDefinition("D0", "D0", 0.0, 2.0, True, "d0", "D"),
    ParticleKindDefinition("D+", "D+", 1.0, 2.0, True, "dplus", "D"),
    ParticleKindDefinition("D2+", "D2+", 1.0, 4.0, True, "d2plus", "D"),
    ParticleKindDefinition("D20", "D20", 0.0, 4.0, True, "d20", "D"),
    ParticleKindDefinition("D3+", "D3+", 1.0, 6.0, True, "d3plus", "D"),
    ParticleKindDefinition(ELECTRON_PARTICLE_KIND, "Electrons", -1.0, 5.48579909065e-4, True, "electrons", None),
)


BUILTIN_PARTICLE_KINDS = {
    definition.kind: definition
    for definition in _BUILTIN_PARTICLE_KIND_DEFINITIONS
}

PARTICLE_TYPE_IDS_TO_KINDS = {
    definition.type_id: definition.kind
    for definition in _BUILTIN_PARTICLE_KIND_DEFINITIONS
}

EXPORTED_PARTICLE_KIND_IDS = (0, 1, 2, 3, 4, 5, 6)

_EXPORTED_PARTICLE_KIND_BASE_KINDS = {
    0: "H-",
    1: "H0",
    2: "H+",
    3: "H2+",
    4: "H20",
    5: "H3+",
    6: ELECTRON_PARTICLE_KIND,
}

_EXPORTED_PARTICLE_KIND_COLORS = {
    0: "#c24d2c",
    1: "#7c9b2e",
    2: "#2a9d8f",
    3: "#3d5a80",
    4: "#8d6a9f",
    5: "#f4a261",
    6: "#4c566a",
}

FAMILY_PARTICLE_KINDS = {
    family: tuple(
        definition.kind
        for definition in _BUILTIN_PARTICLE_KIND_DEFINITIONS
        if definition.family == family
    )
    for family in PARTICLE_FAMILIES
}

SOURCEABLE_PARTICLE_KINDS = {
    definition.kind
    for definition in _BUILTIN_PARTICLE_KIND_DEFINITIONS
    if definition.sourceable
}


def get_particle_kind_definition(kind: str) -> ParticleKindDefinition:
    return BUILTIN_PARTICLE_KINDS[kind]


def get_particle_kind_for_type_id(type_id: str) -> str | None:
    return PARTICLE_TYPE_IDS_TO_KINDS.get(type_id)


def get_particle_kind_family(kind: str) -> str | None:
    definition = BUILTIN_PARTICLE_KINDS.get(kind)
    return definition.family if definition else None


def particle_type_id_for_kind(kind: str) -> str:
    return get_particle_kind_definition(kind).type_id


def get_default_source_particle_type_id(family: str) -> str:
    return particle_type_id_for_kind(family + "-")


def iter_particle_kind_definitions() -> tuple[ParticleKindDefinition, ...]:
    return _BUILTIN_PARTICLE_KIND_DEFINITIONS


def iter_particle_kind_definitions_for_family(
    family: str,
    *,
    include_electrons: bool = True,
) -> tuple[ParticleKindDefinition, ...]:
    kinds = list(FAMILY_PARTICLE_KINDS.get(family, FAMILY_PARTICLE_KINDS[PARTICLE_FAMILIES[0]]))
    if include_electrons:
        kinds.append(ELECTRON_PARTICLE_KIND)
    return tuple(BUILTIN_PARTICLE_KINDS[kind] for kind in kinds)


def build_family_particle_types(
    family: str,
    *,
    include_electrons: bool = True,
) -> list[dict[str, object]]:
    return [
        {
            "id": definition.type_id,
            "name": definition.label,
            "kind": definition.kind,
        }
        for definition in iter_particle_kind_definitions_for_family(
            family,
            include_electrons=include_electrons,
        )
    ]


def detect_particle_family_from_kinds(kinds: list[str]) -> str:
    for kind in kinds:
        family = get_particle_kind_family(kind)
        if family in PARTICLE_FAMILIES:
            return family
    return PARTICLE_FAMILIES[0]


def map_particle_kind_to_family(kind: str, family: str) -> str:
    if kind == ELECTRON_PARTICLE_KIND:
        return kind
    if family not in PARTICLE_FAMILIES:
        return kind

    current_family = get_particle_kind_family(kind)
    if current_family is None or current_family == family:
        return kind

    return family + kind[1:]


def particle_kind_from_export_id(kind_id: int, family: str = "H") -> str | None:
    base_kind = _EXPORTED_PARTICLE_KIND_BASE_KINDS.get(int(kind_id))
    if base_kind is None:
        return None
    return map_particle_kind_to_family(base_kind, family)


def particle_label_from_export_id(kind_id: int, family: str = "H") -> str:
    kind = particle_kind_from_export_id(kind_id, family)
    if kind is None:
        return f"Species {int(kind_id)}"
    return get_particle_kind_definition(kind).label


def particle_color_from_export_id(kind_id: int) -> str:
    return _EXPORTED_PARTICLE_KIND_COLORS.get(int(kind_id), "#52606d")