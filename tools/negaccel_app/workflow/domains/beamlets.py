"""Multi-beamlet extraction: the source-side twin of the grid aperture pattern.

The grid solids already tile their bore profile over a lattice (GeometryManager.cpp,
``evaluateGeneratedSolidProfile``). Nothing tiled the *beam*, so a patterned geometry had to
be fed by hand-placed sources, one per aperture, with nothing checking that the two agreed.

This module expands one authored source into one source per beamlet, on a lattice that
follows the same rule the grid geometry uses, so every beamlet is emitted into a real
aperture. The physics of a beamlet is untouched: each copy keeps the current density,
particle count, energy, temperatures and emission rectangle of the source it came from, and
only its centre moves.
"""

from __future__ import annotations

import copy
from typing import Any

from ..common import WorkflowError

SUPPORTED_BEAMLET_LAYOUTS = {"single", "rectangular-grid", "staggered-grid"}


def axis_node_indices(count: int) -> list[int]:
    """Indices of the ``count`` beamlet nodes along one axis.

    Always exactly ``count`` of them, and always whole multiples of the pitch, because that
    is where the grid geometry puts its apertures: ``evaluateGeneratedSolidProfile`` folds a
    point onto the nearest integer multiple of the pitch. Placing a beamlet anywhere else
    fires it at metal.

    For an odd count the group is centred on the lattice origin. For an even count there is
    no whole node at the centre, so the group is floor-centred -- shifted half a pitch toward
    -x -- rather than straddling the origin at half-integer positions where there is no
    aperture. Use the x/y offset to put it where you want it.
    """
    start = -(count // 2)
    return list(range(start, start + count))


def geometry_full_aperture_indices(count: int) -> list[int]:
    """Indices of the apertures the grid geometry actually opens, along one axis.

    The same nodes as :func:`axis_node_indices`, because the geometry centres its patterned
    region on the *group* rather than on the pattern origin: every one of the ``count``
    nodes is therefore a whole aperture, for an even count as well as an odd one.

    It used to differ. The region was measured from the origin, which put its boundary
    exactly on the outermost node of an even group and sliced those apertures in half, so
    only ``count - 1`` were whole. Grid and beam then disagreed for every even count. This
    is kept as a separate name because it answers a different question -- what the metal
    opens, rather than where the beam is aimed -- and a future pattern could again make the
    two differ.
    """
    return axis_node_indices(count)


def beamlet_nodes(pattern: dict[str, Any]) -> list[tuple[float, float, int, int]]:
    """Every beamlet as ``(dx, dy, column_index, row_index)``.

    ``dx``/``dy`` are the transverse offsets in metres; the two indices are the lattice
    node the beamlet sits on, signed and centred the way ``axis_node_indices`` lays them
    out. The indices are carried through to the runtime case and on into the trajectory
    VTK so a row or a column can be selected in ParaView without reconstructing the
    lattice from the centre coordinates.
    """
    layout = str(pattern.get("layout", "single"))
    x_offset = float(pattern.get("xOffsetMeters", 0.0))
    y_offset = float(pattern.get("yOffsetMeters", 0.0))
    if layout == "single":
        return [(x_offset, y_offset, 0, 0)]

    count_x = int(pattern.get("countX", 1))
    count_y = int(pattern.get("countY", 1))
    pitch_x = float(pattern.get("pitchXMeters", 0.0))
    pitch_y = float(pattern.get("pitchYMeters", 0.0))
    row_shift = float(pattern.get("rowShiftXMeters", 0.0)) if layout == "staggered-grid" else 0.0

    nodes: list[tuple[float, float, int, int]] = []
    for row_index in axis_node_indices(count_y):
        # Sign convention copied from the geometry: it shifts the sampling point by
        # +rowShift on even rows, so the apertures themselves sit at -rowShift there.
        row_dx = (-row_shift if row_index % 2 == 0 else row_shift) if layout == "staggered-grid" else 0.0
        for column_index in axis_node_indices(count_x):
            nodes.append(
                (
                    x_offset + row_dx + column_index * pitch_x,
                    y_offset + row_index * pitch_y,
                    column_index,
                    row_index,
                )
            )
    return nodes


def beamlet_offsets(pattern: dict[str, Any]) -> list[tuple[float, float]]:
    """Transverse (dx, dy) offsets of every beamlet, in metres.

    The offsets are absolute lattice positions including the pattern's own x/y offset, and
    are *added* to each authored source centre -- so author the reference source on the
    lattice origin (normally the machine axis) and let this place the copies.
    """
    return [(dx, dy) for dx, dy, _column, _row in beamlet_nodes(pattern)]


def beamlet_count(pattern: dict[str, Any]) -> int:
    return len(beamlet_offsets(pattern))


def normalize_beamlet_pattern(raw: Any) -> dict[str, Any] | None:
    """Validate the authoring block and return it in canonical form.

    Returns None when no multi-beamlet extraction was asked for, so the runtime case is bit
    for bit what it was before this feature existed. A single-beamlet case must never become
    multi-beamlet by defaulting.
    """
    if not isinstance(raw, dict) or not raw.get("enabled", False):
        return None

    layout = str(raw.get("layout", "rectangular-grid") or "rectangular-grid")
    if layout not in SUPPORTED_BEAMLET_LAYOUTS:
        supported = ", ".join(sorted(SUPPORTED_BEAMLET_LAYOUTS))
        raise WorkflowError(
            f"Unsupported particles.beamlets.layout '{layout}'. Supported layouts: {supported}"
        )

    pattern: dict[str, Any] = {
        "enabled": True,
        "layout": layout,
        "xOffsetMeters": _number(raw, "xOffsetMeters", 0.0),
        "yOffsetMeters": _number(raw, "yOffsetMeters", 0.0),
    }

    if layout == "single":
        return pattern

    count_x = _positive_integer(raw, "countX")
    count_y = _positive_integer(raw, "countY")
    pitch_x = _number(raw, "pitchXMeters", None)
    pitch_y = _number(raw, "pitchYMeters", None)
    if pitch_x <= 0.0 or pitch_y <= 0.0:
        raise WorkflowError(
            "particles.beamlets.pitchXMeters and pitchYMeters must be greater than zero "
            "for a patterned layout"
        )

    pattern.update(
        {
            "countX": count_x,
            "countY": count_y,
            "pitchXMeters": pitch_x,
            "pitchYMeters": pitch_y,
        }
    )
    if layout == "staggered-grid":
        pattern["rowShiftXMeters"] = _number(raw, "rowShiftXMeters", 0.0)

    return pattern


def expand_sources_over_beamlets(
    sources: list[dict[str, Any]],
    pattern: dict[str, Any],
) -> list[dict[str, Any]]:
    """Replicate every source once per beamlet, moving only its centre.

    ``particleCount`` and ``currentDensityAm2`` are per beamlet, deliberately: a beamlet in
    an N-beamlet run is then statistically identical to the same beamlet run on its own, so
    an existing single-beamlet result stays the reference. The cost of that is that both the
    total current and the total particle count scale with the number of beamlets, which is
    what the GUI reports before the run starts.

    Each copy is stamped with ``beamletIndex`` (1-based, in this list's order) and its
    lattice node ``beamletColumn``/``beamletRow``. The C++ runtime reads those back and
    labels every macroparticle the source emits with them, so a trajectory can be traced to
    the aperture it came out of. Without the stamp the only link between a particle and its
    beamlet is the emission coordinate, which stops being reliable as soon as stripping
    starts creating secondaries away from the source plane.
    """
    nodes = beamlet_nodes(pattern)
    if len(nodes) <= 1:
        return sources

    digits = len(str(len(nodes)))
    expanded: list[dict[str, Any]] = []
    for source in sources:
        base_id = str(source.get("id") or "source")
        base_name = str(source.get("name") or base_id)
        for index, (dx, dy, column_index, row_index) in enumerate(nodes, start=1):
            clone = copy.deepcopy(source)
            clone["id"] = f"{base_id}-b{index:0{digits}d}"
            clone["name"] = f"{base_name} [beamlet {index}]"
            center = list(clone["uniform"]["centerMeters"])
            clone["uniform"]["centerMeters"] = [center[0] + dx, center[1] + dy, center[2]]
            clone["beamletIndex"] = index
            clone["beamletColumn"] = column_index
            clone["beamletRow"] = row_index
            expanded.append(clone)
    return expanded


def describe_beamlet_pattern(pattern: dict[str, Any] | None) -> str:
    """One line for the GUI and the run log."""
    if pattern is None:
        return "Single beamlet."

    offsets = beamlet_offsets(pattern)
    if len(offsets) <= 1:
        return "Single beamlet."

    xs = [offset[0] for offset in offsets]
    ys = [offset[1] for offset in offsets]
    return (
        "{count} beamlets ({layout}), spanning {span_x:.1f} x {span_y:.1f} mm".format(
            count=len(offsets),
            layout=pattern["layout"],
            span_x=1.0e3 * (max(xs) - min(xs)),
            span_y=1.0e3 * (max(ys) - min(ys)),
        )
    )


def beamlet_pattern_warnings(
    pattern: dict[str, Any] | None,
    aperture_pattern: dict[str, Any] | None = None,
) -> list[str]:
    """Problems worth showing the user that are not hard errors."""
    if pattern is None or str(pattern.get("layout")) == "single":
        return []

    warnings: list[str] = []
    for axis, key, pitch_key, offset_key in (
        ("X", "countX", "pitchXMeters", "xOffsetMeters"),
        ("Y", "countY", "pitchYMeters", "yOffsetMeters"),
    ):
        count = int(pattern.get(key, 1))
        if count <= 1 or count % 2 != 0:
            continue
        pitch = float(pattern.get(pitch_key, 0.0) or 0.0)
        offset = float(pattern.get(offset_key, 0.0) or 0.0)
        # An even group has no node at the pattern origin, so it is floor-centred: its own
        # centre is half a pitch below the origin. Adding half a pitch of offset therefore
        # puts it back on the machine axis -- so say nothing when the author already has.
        # The group centre is only defined modulo the pitch, hence the remainder test.
        centred = pitch > 0.0 and abs((offset % pitch) - 0.5 * pitch) <= 1.0e-9 * pitch
        if centred:
            continue
        warnings.append(
            "count{axis} is even ({count}), so the group has no beamlet on the axis and "
            "sits half a pitch ({half:.1f} mm) toward -{axis}. The grid apertures follow "
            "the same rule, so beam and metal stay aligned, but the whole group is off "
            "the machine axis and the pooled diagnostics assume it is not. Set the "
            "{axis} offset to {half:.1f} mm to centre it.".format(
                axis=axis, count=count, half=0.5e3 * pitch,
            )
        )

    # The nodes must be apertures the grid geometry actually opens, or the beamlet is fired
    # at metal. Nothing downstream checks this: the run completes and the beamlet is simply
    # lost on the plasma grid.
    if isinstance(aperture_pattern, dict) and str(aperture_pattern.get("layout", "single")) != "single":
        for axis, key in (("X", "countX"), ("Y", "countY")):
            open_nodes = set(geometry_full_aperture_indices(int(aperture_pattern.get(key, 1) or 1)))
            missing = [k for k in axis_node_indices(int(pattern.get(key, 1))) if k not in open_nodes]
            if missing:
                warnings.append(
                    "Along {axis} the grid opens apertures only at {open}; beamlets are asked "
                    "for at {missing}, which is solid metal. Widen the grid aperture pattern "
                    "or shrink the beamlet group.".format(
                        axis=axis,
                        open=sorted(open_nodes),
                        missing=missing,
                    )
                )
    return warnings


def _number(document: dict[str, Any], key: str, default: float | None) -> float:
    value = document.get(key, default)
    if value is None:
        raise WorkflowError(f"particles.beamlets requires numeric '{key}'")
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise WorkflowError(f"particles.beamlets.{key} must be numeric")
    return float(value)


def _positive_integer(document: dict[str, Any], key: str) -> int:
    value = document.get(key, 1)
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise WorkflowError(f"particles.beamlets.{key} must be an integer >= 1")
    return int(value)
