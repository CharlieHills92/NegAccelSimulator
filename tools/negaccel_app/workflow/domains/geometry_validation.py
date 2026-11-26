"""Shared geometry-aware validation and diagnostics derivation helpers."""

from __future__ import annotations

import math
from typing import Any

from ..common import WorkflowError


_DOMAIN_TOLERANCE = 1.0e-12

# IBSimu reserves boundaries 1-6 for the six faces of the domain; solids start at 7.
FIRST_SOLID_BOUNDARY_ID = 7


def resolve_solid_boundary_ids(solids: Any, context: str = "geometry.solids") -> list[int]:
    """The IBSimu boundary id of each solid, in list order.

    An explicit ``boundaryId`` is honoured exactly as written, including when two solids
    share one: that is how a split electrode -- EG and SG, say -- is held on a single
    potential and shown as one boundary. A solid that declares none is given the lowest id
    from 7 up that no other solid has taken.

    This replaces the rule ``7 + index`` that each caller used to apply on its own. The two
    agree only while every solid has a distinct id. The moment two solids share one, the
    index rule skips a number, and ``Geometry::set_solid`` rejects a non-consecutive id
    ("illegal solid number N") because IBSimu appends solids one at a time -- so a case with
    a shared electrode plus any later unnumbered solid died during geometry construction.
    """
    if not isinstance(solids, list):
        raise WorkflowError(f"{context} must be an array")

    explicit: list[int | None] = []
    for index, solid in enumerate(solids):
        if not isinstance(solid, dict):
            raise WorkflowError(f"{context}[{index}] must be an object")
        raw = solid.get("boundaryId")
        if raw is None:
            explicit.append(None)
            continue
        if isinstance(raw, bool) or not isinstance(raw, int):
            raise WorkflowError(f"{context}[{index}].boundaryId must be an integer")
        if raw < FIRST_SOLID_BOUNDARY_ID:
            raise WorkflowError(
                f"{context}[{index}].boundaryId must be >= {FIRST_SOLID_BOUNDARY_ID}; "
                "1-6 are reserved for the domain faces"
            )
        explicit.append(raw)

    taken = {value for value in explicit if value is not None}
    resolved: list[int] = []
    next_free = FIRST_SOLID_BOUNDARY_ID
    for value in explicit:
        if value is not None:
            resolved.append(value)
            continue
        while next_free in taken:
            next_free += 1
        taken.add(next_free)
        resolved.append(next_free)
    return resolved


def validate_solid_boundary_ids(solids: Any, context: str = "geometry.solids") -> list[int]:
    """Resolve the ids and check IBSimu will accept them.

    ``Geometry::set_solid`` grows its solid table by one per call, so it only accepts ids
    that are consecutive from 7 upwards. A gap is not a warning there, it is a fatal error
    part-way through building the geometry -- after the mesh has been allocated. Catching it
    here costs nothing and names the solid that caused it.
    """
    resolved = resolve_solid_boundary_ids(solids, context)
    distinct = sorted(set(resolved))
    expected = list(range(FIRST_SOLID_BOUNDARY_ID, FIRST_SOLID_BOUNDARY_ID + len(distinct)))
    if distinct != expected:
        missing = sorted(set(expected) - set(distinct))
        names = ", ".join(
            "{} (id {})".format(solids[i].get("name") or f"solids[{i}]", resolved[i])
            for i in range(len(resolved))
        )
        raise WorkflowError(
            "geometry solid boundary ids must be consecutive from "
            f"{FIRST_SOLID_BOUNDARY_ID}; got {distinct} with {missing} unused. "
            "IBSimu adds solids one at a time and refuses a gap. Solids: " + names
        )
    return resolved


def get_domain_z_bounds(geometry: dict[str, Any]) -> tuple[float, float]:
    domain = geometry.get("domain")
    if not isinstance(domain, dict):
        raise WorkflowError("geometry.domain must be an object")

    z_size = _require_number(domain.get("zSizeMeters"), "geometry.domain.zSizeMeters")
    if z_size <= 0.0:
        raise WorkflowError("geometry.domain.zSizeMeters must be greater than zero")

    z_start = _require_number(domain.get("zStartMeters", 0.0), "geometry.domain.zStartMeters")
    return z_start, z_start + z_size


def get_solid_z_extents(geometry: dict[str, Any]) -> dict[int, tuple[float, float]]:
    solids = geometry.get("solids")
    if not isinstance(solids, list):
        raise WorkflowError("geometry.solids must be an array")

    resolved_ids = resolve_solid_boundary_ids(solids)
    extents: dict[int, tuple[float, float]] = {}
    for index, solid in enumerate(solids):
        context = f"geometry.solids[{index}]"
        if not isinstance(solid, dict):
            raise WorkflowError(f"{context} must be an object")
        boundary_id = resolved_ids[index]

        z_profile = solid.get("zProfileMeters")
        if not isinstance(z_profile, list) or len(z_profile) < 2:
            raise WorkflowError(f"{context}.zProfileMeters must contain at least two entries")

        z_values = [_require_number(value, f"{context}.zProfileMeters[{point_index}]") for point_index, value in enumerate(z_profile)]
        solid_min = min(z_values)
        solid_max = max(z_values)
        if boundary_id in extents:
            current_min, current_max = extents[boundary_id]
            extents[boundary_id] = (min(current_min, solid_min), max(current_max, solid_max))
        else:
            extents[boundary_id] = (solid_min, solid_max)
    return extents


def derive_diagnostic_planes(geometry: dict[str, Any]) -> dict[str, float | list[float]]:
    domain_z_min, domain_z_max = get_domain_z_bounds(geometry)
    solid_extents = get_solid_z_extents(geometry)

    sample_planes: list[float] = []
    for boundary_id in sorted(solid_extents):
        solid_min, solid_max = solid_extents[boundary_id]
        sample_planes.append(solid_min)
        if solid_max != solid_min:
            sample_planes.append(solid_max)

    deduplicated_planes = _deduplicate_sorted(sample_planes)

    transmission_plane = domain_z_max
    pg_extent = solid_extents.get(7)
    if pg_extent is not None:
        transmission_plane = pg_extent[1]

    return {
        "sampleZPositionsMeters": deduplicated_planes,
        "summaryZPositionMeters": domain_z_max,
        "emitterExportZPositionMeters": domain_z_max,
        "transmissionPlaneZPositionMeters": transmission_plane,
        "domainZMinMeters": domain_z_min,
        "domainZMaxMeters": domain_z_max,
    }


def validate_diagnostic_planes(geometry: dict[str, Any], diagnostics: dict[str, Any]) -> None:
    domain_z_min, domain_z_max = get_domain_z_bounds(geometry)
    planes = diagnostics.get("planes")
    if not isinstance(planes, dict):
        raise WorkflowError("diagnostics.planes must be an object")
    summary = diagnostics.get("summary")
    if not isinstance(summary, dict):
        raise WorkflowError("diagnostics.summary must be an object")

    sample_z_positions = planes.get("sampleZPositionsMeters")
    if not isinstance(sample_z_positions, list):
        raise WorkflowError("diagnostics.planes.sampleZPositionsMeters must be an array")
    for index, z_value in enumerate(sample_z_positions):
        _require_z_in_domain(
            _require_number(z_value, f"diagnostics.planes.sampleZPositionsMeters[{index}]"),
            domain_z_min,
            domain_z_max,
            f"diagnostics.planes.sampleZPositionsMeters[{index}]",
        )

    _require_z_in_domain(
        _require_number(planes.get("summaryZPositionMeters"), "diagnostics.planes.summaryZPositionMeters"),
        domain_z_min,
        domain_z_max,
        "diagnostics.planes.summaryZPositionMeters",
    )
    _require_z_in_domain(
        _require_number(planes.get("emitterExportZPositionMeters"), "diagnostics.planes.emitterExportZPositionMeters"),
        domain_z_min,
        domain_z_max,
        "diagnostics.planes.emitterExportZPositionMeters",
    )
    _require_z_in_domain(
        _require_number(summary.get("transmissionPlaneZPositionMeters"), "diagnostics.summary.transmissionPlaneZPositionMeters"),
        domain_z_min,
        domain_z_max,
        "diagnostics.summary.transmissionPlaneZPositionMeters",
    )


def validate_uniform_sources_within_domain(geometry: dict[str, Any], particle_sources: list[dict[str, Any]]) -> None:
    domain = geometry.get("domain")
    if not isinstance(domain, dict):
        raise WorkflowError("geometry.domain must be an object")

    domain_x_size = _require_number(domain.get("xSizeMeters"), "geometry.domain.xSizeMeters")
    domain_y_size = _require_number(domain.get("ySizeMeters"), "geometry.domain.ySizeMeters")
    domain_z_min, domain_z_max = get_domain_z_bounds(geometry)

    x_limits = (-0.5 * domain_x_size, 0.5 * domain_x_size)
    y_limits = (-0.5 * domain_y_size, 0.5 * domain_y_size)
    for index, source in enumerate(particle_sources):
        if not isinstance(source, dict):
            continue
        uniform = source.get("uniform")
        if not isinstance(uniform, dict):
            continue
        center = uniform.get("centerMeters")
        if not isinstance(center, list) or len(center) != 3:
            continue
        main_direction = uniform.get("mainDirection")
        reference_direction = uniform.get("inPlaneReferenceDirection")
        if not isinstance(main_direction, list) or len(main_direction) != 3:
            continue
        if not isinstance(reference_direction, list) or len(reference_direction) != 3:
            continue
        height = uniform.get("heightMeters")
        width = uniform.get("widthMeters")
        if isinstance(height, bool) or not isinstance(height, (int, float)):
            continue
        if isinstance(width, bool) or not isinstance(width, (int, float)):
            continue

        center_vector = [
            _require_number(center[0], f"particles.sources[{index}].uniform.centerMeters[0]"),
            _require_number(center[1], f"particles.sources[{index}].uniform.centerMeters[1]"),
            _require_number(center[2], f"particles.sources[{index}].uniform.centerMeters[2]"),
        ]
        normalized_main = _normalize_vector(
            [_require_number(component, f"particles.sources[{index}].uniform.mainDirection[{component_index}]")
             for component_index, component in enumerate(main_direction)],
            f"particles.sources[{index}].uniform.mainDirection",
        )
        projected_reference = _project_vector_onto_plane(
            [_require_number(component, f"particles.sources[{index}].uniform.inPlaneReferenceDirection[{component_index}]")
             for component_index, component in enumerate(reference_direction)],
            normalized_main,
        )
        reference_axis = _normalize_vector(
            projected_reference,
            f"particles.sources[{index}].uniform.inPlaneReferenceDirection",
            zero_message="must not be parallel to mainDirection",
        )
        secondary_axis = _normalize_vector(
            _cross_product(normalized_main, reference_axis),
            f"particles.sources[{index}].uniform.inPlaneReferenceDirection",
        )

        half_width_vector = [0.5 * float(width) * component for component in reference_axis]
        half_height_vector = [0.5 * float(height) * component for component in secondary_axis]

        corner_coordinates = []
        for width_sign in (-1.0, 1.0):
            for height_sign in (-1.0, 1.0):
                corner_coordinates.append(
                    [
                        center_vector[0] + width_sign * half_width_vector[0] + height_sign * half_height_vector[0],
                        center_vector[1] + width_sign * half_width_vector[1] + height_sign * half_height_vector[1],
                        center_vector[2] + width_sign * half_width_vector[2] + height_sign * half_height_vector[2],
                    ]
                )

        source_x_min = min(corner[0] for corner in corner_coordinates)
        source_x_max = max(corner[0] for corner in corner_coordinates)
        source_y_min = min(corner[1] for corner in corner_coordinates)
        source_y_max = max(corner[1] for corner in corner_coordinates)
        source_z_min = min(corner[2] for corner in corner_coordinates)
        source_z_max = max(corner[2] for corner in corner_coordinates)

        if (
            source_x_min < x_limits[0] - _DOMAIN_TOLERANCE
            or source_x_max > x_limits[1] + _DOMAIN_TOLERANCE
            or source_y_min < y_limits[0] - _DOMAIN_TOLERANCE
            or source_y_max > y_limits[1] + _DOMAIN_TOLERANCE
            or source_z_min < domain_z_min - _DOMAIN_TOLERANCE
            or source_z_max > domain_z_max + _DOMAIN_TOLERANCE
        ):
            # Name the source as well as its index: with multi-beamlet extraction the index
            # points into the expanded list, where "sources[7]" means nothing to the author.
            raise WorkflowError(
                "particles.sources[{}]{}.uniform spans x = [{:.6g}, {:.6g}] m, y = [{:.6g}, {:.6g}] m, z = [{:.6g}, {:.6g}] m, outside geometry.domain bounds x = [{:.6g}, {:.6g}] m, y = [{:.6g}, {:.6g}] m, z = [{:.6g}, {:.6g}] m. The current runtime seeds particles exactly at this source plane, so the full source footprint must lie inside geometry.domain; otherwise the particles become bad definitions before tracing begins.".format(
                    index,
                    _source_label(source),
                    source_x_min,
                    source_x_max,
                    source_y_min,
                    source_y_max,
                    source_z_min,
                    source_z_max,
                    x_limits[0],
                    x_limits[1],
                    y_limits[0],
                    y_limits[1],
                    domain_z_min,
                    domain_z_max,
                )
            )


def _require_number(value: Any, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise WorkflowError(f"{context} must be numeric")
    return float(value)


def _deduplicate_sorted(values: list[float], tolerance: float = 1.0e-12) -> list[float]:
    ordered = sorted(values)
    deduplicated: list[float] = []
    for value in ordered:
        if not deduplicated or abs(deduplicated[-1] - value) > tolerance:
            deduplicated.append(value)
    return deduplicated


def _require_z_in_domain(value: float, domain_z_min: float, domain_z_max: float, context: str) -> None:
    if value < domain_z_min - _DOMAIN_TOLERANCE or value > domain_z_max + _DOMAIN_TOLERANCE:
        raise WorkflowError(
            f"{context}={value:g} m is outside geometry.domain z range [{domain_z_min:g}, {domain_z_max:g}] m"
        )


def _normalize_vector(vector: list[float], context: str, zero_message: str = "must not be the zero vector") -> list[float]:
    magnitude = math.sqrt(sum(component * component for component in vector))
    if magnitude <= 1.0e-15:
        raise WorkflowError(f"{context} {zero_message}")
    return [component / magnitude for component in vector]


def _project_vector_onto_plane(vector: list[float], plane_normal: list[float]) -> list[float]:
    projection = sum(component * normal for component, normal in zip(vector, plane_normal))
    return [component - projection * normal for component, normal in zip(vector, plane_normal)]


def _cross_product(left: list[float], right: list[float]) -> list[float]:
    return [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]

def _source_label(source: dict[str, Any]) -> str:
    """Empty for an unnamed source, " (id 'x')" otherwise.

    With multi-beamlet extraction the index points into the expanded list, where
    "sources[7]" on its own tells the author nothing about which beamlet is out of bounds.
    """
    source_id = str(source.get("id") or "").strip()
    return f" (id '{source_id}')" if source_id else ""
