"""Visualization helpers for NegAccel GUI tooling."""

from __future__ import annotations

import contextlib
from dataclasses import dataclass
import importlib.util
import io
import math
import re
from pathlib import Path
from typing import Any

try:
    import matplotlib.cm as cm
    import numpy as np
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize
    from matplotlib.figure import Figure
    from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
except ImportError as exc:  # pragma: no cover - depends on the local environment
    raise SystemExit(
        "numpy and matplotlib are required for the GUI. Install the packages listed in tools/gui_requirements.txt."
    ) from exc

from .common import REPO_ROOT, Qt
from negaccel_app.particles import (
    PARTICLE_FAMILIES,
    particle_color_from_export_id,
    particle_label_from_export_id,
)
from negaccel_app.workflow.post_processing import parse_diagnostic_summary_txt, parse_grid_power_summary_txt


PLOT_UTILS_PATH = REPO_ROOT / "plot_trajectories_vtk.py"


def load_plot_utils() -> Any:
    spec = importlib.util.spec_from_file_location("negaccel_plot_utils", PLOT_UTILS_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load plotting helpers from {PLOT_UTILS_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PLOT_UTILS = load_plot_utils()

GEOMETRY_FILL_COLOR = "#b8b8b8"
GEOMETRY_EDGE_COLOR = "#6e6e6e"
DEFAULT_TRAJECTORY_COLOR = "#2f6f7e"
TRAJECTORY_COLOR_MODE_SINGLE = "Single color"
TRAJECTORY_COLOR_MODE_SPECIES = "Particle species"
TRAJECTORY_COLOR_MODE_ENERGY = "Particle energy"
TRAJECTORY_COLOR_MODE_STATUS = "Particle status"


@dataclass(frozen=True)
class TrajectoryRenderOptions:
    color_mode: str = TRAJECTORY_COLOR_MODE_SINGLE
    selected_particle_kinds: tuple[int, ...] | None = None
    energy_min_ev: float | None = None
    energy_max_ev: float | None = None
    particle_family: str = "H"


@dataclass(frozen=True)
class TrajectoryMetadata:
    source_path: Path
    available_particle_kinds: tuple[int, ...]
    trajectory_count: int
    energy_min_ev: float | None
    energy_max_ev: float | None


def _normalize_particle_family(family: str | None) -> str:
    candidate = (family or "").strip()
    return candidate if candidate in PARTICLE_FAMILIES else PARTICLE_FAMILIES[0]


def _normalize_trajectory_options(
    options: TrajectoryRenderOptions | None,
) -> TrajectoryRenderOptions:
    if options is None:
        return TrajectoryRenderOptions()
    return TrajectoryRenderOptions(
        color_mode=str(options.color_mode or TRAJECTORY_COLOR_MODE_SINGLE),
        selected_particle_kinds=options.selected_particle_kinds,
        energy_min_ev=options.energy_min_ev,
        energy_max_ev=options.energy_max_ev,
        particle_family=_normalize_particle_family(options.particle_family),
    )


def _trajectory_plane_crossing(
    points: np.ndarray,
    line_indices: np.ndarray,
    plane_axis: str,
    plane_value: float,
) -> tuple[int, np.ndarray, np.ndarray, float] | None:
    plane_key = plane_axis.lower()
    axis_index = {"x": 0, "y": 1, "z": 2}.get(plane_key, 2)

    trajectory_points = points[line_indices]
    for point_index in range(len(trajectory_points) - 1):
        point_a = trajectory_points[point_index]
        point_b = trajectory_points[point_index + 1]
        coord_a = point_a[axis_index]
        coord_b = point_b[axis_index]
        if not ((coord_a <= plane_value <= coord_b) or (coord_b <= plane_value <= coord_a)):
            continue

        if abs(coord_b - coord_a) < 1.0e-12:
            return point_index, point_a, point_b, 0.5

        return point_index, point_a, point_b, (plane_value - coord_a) / (coord_b - coord_a)

    return None


def _interpolate_scalar_value(value_a: float, value_b: float, fraction: float) -> float:
    return float(value_a + fraction * (value_b - value_a))


def _resolve_particle_kind_from_cell_data(
    cell_data: dict[str, np.ndarray],
    trajectory_index: int,
    particle_family: str,
) -> int | None:
    particle_kinds = cell_data.get("particle_kind")
    if particle_kinds is not None and trajectory_index < len(particle_kinds):
        try:
            return int(particle_kinds[trajectory_index])
        except (TypeError, ValueError):
            pass

    charge_values = cell_data.get("charge_e")
    mass_values = cell_data.get("mass_amu")
    if charge_values is None or mass_values is None:
        return None
    if trajectory_index >= len(charge_values) or trajectory_index >= len(mass_values):
        return None

    try:
        charge_e = float(charge_values[trajectory_index])
        mass_amu = float(mass_values[trajectory_index])
    except (TypeError, ValueError):
        return None

    if not math.isfinite(charge_e) or not math.isfinite(mass_amu):
        return None
    if abs(mass_amu) < 0.05:
        return 6

    ion_mass = 2.0 if _normalize_particle_family(particle_family) == "D" else 1.0
    mass_ratio = mass_amu / ion_mass
    if 0.75 <= mass_ratio <= 1.25:
        if charge_e < -0.5:
            return 0
        if charge_e > 0.5:
            return 2
        return 1
    if 1.75 <= mass_ratio <= 2.25:
        return 3 if charge_e > 0.5 else 4
    if 2.75 <= mass_ratio <= 3.25 and charge_e > 0.5:
        return 5
    return None


def _trajectory_energy_values(
    point_data: dict[str, np.ndarray],
    line_indices: np.ndarray,
) -> np.ndarray | None:
    raw_energy = point_data.get("kinetic_energy_eV")
    if raw_energy is None:
        return None
    energy_array = np.asarray(raw_energy, dtype=float)
    if energy_array.size == 0:
        return None
    return np.asarray(energy_array[line_indices], dtype=float)


def _trajectory_visibility_mask(
    line_length: int,
    energies: np.ndarray | None,
    options: TrajectoryRenderOptions,
) -> np.ndarray:
    visible = np.ones(line_length, dtype=bool)
    if energies is None:
        return visible

    visible &= np.isfinite(energies)
    if options.energy_min_ev is not None:
        visible &= energies >= float(options.energy_min_ev)
    if options.energy_max_ev is not None:
        visible &= energies <= float(options.energy_max_ev)
    return visible


def _split_visible_polylines(points: np.ndarray, visible: np.ndarray) -> list[np.ndarray]:
    segments: list[np.ndarray] = []
    start_index: int | None = None
    for index, is_visible in enumerate(np.asarray(visible, dtype=bool)):
        if is_visible:
            if start_index is None:
                start_index = index
            continue
        if start_index is not None and index - start_index >= 2:
            segments.append(np.asarray(points[start_index:index], dtype=float))
        start_index = None
    if start_index is not None and len(points) - start_index >= 2:
        segments.append(np.asarray(points[start_index:], dtype=float))
    return segments


def _visible_segment_pairs(
    points: np.ndarray,
    visible: np.ndarray,
    values: np.ndarray | None = None,
) -> tuple[list[np.ndarray], list[float]]:
    segments: list[np.ndarray] = []
    segment_values: list[float] = []
    for index in range(len(points) - 1):
        if not (visible[index] and visible[index + 1]):
            continue
        segments.append(np.asarray([points[index], points[index + 1]], dtype=float))
        if values is not None:
            segment_values.append(float(0.5 * (values[index] + values[index + 1])))
    return segments, segment_values


def _trajectory_color_mode_key(color_mode: str) -> str:
    return color_mode.strip().lower()


def _build_normalize(vmin: float, vmax: float) -> Normalize:
    lower = float(vmin)
    upper = float(vmax)
    if upper <= lower:
        upper = lower + max(abs(lower) * 1.0e-12, 1.0e-12)
    return Normalize(vmin=lower, vmax=upper)


def _add_species_legend(
    axis: Any,
    particle_kind_ids: set[int],
    particle_family: str,
    *,
    marker: bool,
    three_dimensional: bool = False,
) -> None:
    for kind_id in sorted(particle_kind_ids):
        label = particle_label_from_export_id(kind_id, particle_family)
        color = particle_color_from_export_id(kind_id)
        if marker:
            axis.scatter([], [], color=color, s=34, alpha=0.9, label=label)
        elif three_dimensional:
            axis.plot([], [], [], color=color, linewidth=1.6, alpha=0.9, label=label)
        else:
            axis.plot([], [], color=color, linewidth=1.6, alpha=0.9, label=label)


def configure_matplotlib_canvas(widget: FigureCanvasQTAgg) -> None:
    widget.setProperty("matplotlibCanvas", True)
    widget.setAutoFillBackground(False)
    widget.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, False)
    widget.setStyleSheet("background: transparent;")


def configure_matplotlib_toolbar(toolbar: NavigationToolbar2QT) -> None:
    toolbar.setProperty("matplotlibToolbar", True)
    toolbar.setAutoFillBackground(False)
    toolbar.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, False)
    toolbar.setStyleSheet(
        """
        QToolBar {
            background: transparent;
            border: none;
            spacing: 2px;
        }
        QToolButton {
            background: transparent;
            border: 1px solid transparent;
            border-radius: 4px;
            padding: 4px;
        }
        QToolButton:hover {
            background: rgba(194, 77, 44, 0.10);
            border-color: rgba(194, 77, 44, 0.30);
        }
        QToolButton:pressed {
            background: rgba(194, 77, 44, 0.18);
        }
        """
    )


def _set_adaptive_3d_bounds(
    axis: Any,
    point_cloud: np.ndarray,
    scale: float = 1.0,
    x_limits: tuple[float, float] | None = None,
    y_limits: tuple[float, float] | None = None,
    z_limits: tuple[float, float] | None = None,
) -> None:
    if point_cloud.size == 0 and x_limits is None and y_limits is None and z_limits is None:
        return

    if point_cloud.size == 0:
        scaled_points = np.zeros((1, 3), dtype=float)
    else:
        scaled_points = np.asarray(point_cloud, dtype=float) * scale

    minima = scaled_points.min(axis=0)
    maxima = scaled_points.max(axis=0)
    ranges = maxima - minima
    padding = np.maximum(ranges * 0.05, 0.5)

    computed_x_limits = x_limits or (minima[0] - padding[0], maxima[0] + padding[0])
    computed_y_limits = y_limits or (minima[1] - padding[1], maxima[1] + padding[1])
    computed_z_limits = z_limits or (minima[2] - padding[2], maxima[2] + padding[2])

    axis.set_xlim(*computed_x_limits)
    axis.set_ylim(*computed_y_limits)
    axis.set_zlim(*computed_z_limits)

    if hasattr(axis, "set_box_aspect"):
        axis.set_box_aspect(
            np.maximum(
                np.asarray(
                    [
                        computed_x_limits[1] - computed_x_limits[0],
                        computed_y_limits[1] - computed_y_limits[0],
                        computed_z_limits[1] - computed_z_limits[0],
                    ],
                    dtype=float,
                ),
                1.0,
            )
        )


def _set_equal_spatial_aspect(axis: Any) -> None:
    axis.set_aspect("equal", adjustable="box")


def _interpolate_trajectory_at_plane(
    points: np.ndarray,
    line_indices: np.ndarray,
    plane_axis: str,
    plane_value: float,
) -> tuple[float, float] | None:
    plane_key = plane_axis.lower()
    projected_axes = {
        "x": (2, 1),
        "y": (2, 0),
        "z": (0, 1),
    }[plane_key]

    crossing = _trajectory_plane_crossing(points, line_indices, plane_axis, plane_value)
    if crossing is None:
        return None

    _point_index, point_a, point_b, fraction = crossing
    return (
        _interpolate_scalar_value(point_a[projected_axes[0]], point_b[projected_axes[0]], fraction),
        _interpolate_scalar_value(point_a[projected_axes[1]], point_b[projected_axes[1]], fraction),
    )


class GeometryCanvas(FigureCanvasQTAgg):
    def __init__(self) -> None:
        self.figure = Figure(figsize=(8, 6), tight_layout=True)
        super().__init__(self.figure)
        configure_matplotlib_canvas(self)
        self.plot_placeholder("Select a geometry preview mode to inspect the current solid setup.")

    def _redraw(self) -> None:
        self.draw_idle()

    def plot_placeholder(self, message: str) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.axis("off")
        axis.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
        self._redraw()

    def _solid_profile_arrays(self, solid: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        z_values = np.asarray(solid.get("zProfileMeters") or [], dtype=float)
        r_values = np.asarray(solid.get("rProfileMeters") or [], dtype=float)
        rounding_values = np.asarray(solid.get("roundingRadiiMeters") or [], dtype=float)
        if rounding_values.size != z_values.size:
            rounding_values = np.zeros_like(z_values)
        return z_values, r_values, rounding_values

    def _normalize_angle(self, angle: float) -> float:
        while angle <= -np.pi:
            angle += 2.0 * np.pi
        while angle > np.pi:
            angle -= 2.0 * np.pi
        return float(angle)

    def _positive_angle_distance(self, start_angle: float, end_angle: float) -> float:
        delta = self._normalize_angle(end_angle - start_angle)
        if delta < 0.0:
            delta += 2.0 * np.pi
        return float(delta)

    def _angle_on_minor_arc(self, angle: float, start_angle: float, sweep: float) -> bool:
        tolerance = 1.0e-9
        if sweep >= 0.0:
            return self._positive_angle_distance(start_angle, angle) <= sweep + tolerance
        return self._positive_angle_distance(angle, start_angle) <= -sweep + tolerance

    def _segment_boundary_radius(
        self,
        z_point: float,
        z0: float,
        r0: float,
        z1: float,
        r1: float,
        rounding_value: float,
    ) -> float | None:
        radius = abs(float(rounding_value))
        if radius <= 0.0:
            return None

        dz = z1 - z0
        dr = r1 - r0
        chord = float(np.hypot(dz, dr))
        if chord <= 1.0e-12 or chord > 2.0 * radius + 1.0e-12:
            return None

        mid_z = 0.5 * (z0 + z1)
        mid_r = 0.5 * (r0 + r1)
        center_offset = float(np.sqrt(max(0.0, radius * radius - 0.25 * chord * chord)))
        perp_z = -dr / chord
        perp_r = dz / chord

        chosen = None
        for sign in (-1.0, 1.0):
            center_z = mid_z + sign * center_offset * perp_z
            center_r = mid_r + sign * center_offset * perp_r
            start_angle = float(np.arctan2(r0 - center_r, z0 - center_z))
            end_angle = float(np.arctan2(r1 - center_r, z1 - center_z))
            sweep = self._normalize_angle(end_angle - start_angle)
            mid_angle = start_angle + 0.5 * sweep
            mid_r_value = center_r + radius * np.sin(mid_angle)
            if chosen is None:
                chosen = (center_z, center_r, start_angle, sweep, mid_r_value)
                continue
            if rounding_value > 0.0 and mid_r_value < chosen[4]:
                chosen = (center_z, center_r, start_angle, sweep, mid_r_value)
            if rounding_value < 0.0 and mid_r_value > chosen[4]:
                chosen = (center_z, center_r, start_angle, sweep, mid_r_value)

        if chosen is None:
            return None

        center_z, center_r, start_angle, sweep, _mid_r = chosen
        normalized_z = (z_point - center_z) / radius
        if normalized_z < -1.0 - 1.0e-12 or normalized_z > 1.0 + 1.0e-12:
            return None

        clamped_z = max(-1.0, min(1.0, normalized_z))
        base_angle = float(np.arccos(clamped_z))
        selected_radius = None
        for angle in (base_angle, -base_angle):
            if not self._angle_on_minor_arc(angle, start_angle, sweep):
                continue
            candidate_radius = center_r + radius * np.sin(angle)
            if selected_radius is None:
                selected_radius = candidate_radius
            elif rounding_value > 0.0 and candidate_radius < selected_radius:
                selected_radius = candidate_radius
            elif rounding_value < 0.0 and candidate_radius > selected_radius:
                selected_radius = candidate_radius

        return None if selected_radius is None else float(selected_radius)

    def _boundary_radius_at_z(
        self,
        z_point: float,
        z_profile: np.ndarray,
        r_profile: np.ndarray,
        rounding_radii: np.ndarray,
    ) -> float:
        for index in range(1, z_profile.size):
            z0 = float(z_profile[index - 1])
            z1 = float(z_profile[index])
            if abs(z1 - z0) <= 1.0e-15:
                continue
            if z_point < z0 or z_point > z1:
                continue
            r0 = float(r_profile[index - 1])
            r1 = float(r_profile[index])
            rounding_value = float(rounding_radii[index]) if rounding_radii.size == z_profile.size else 0.0
            rounded_boundary = self._segment_boundary_radius(z_point, z0, r0, z1, r1, rounding_value)
            if rounded_boundary is not None:
                return rounded_boundary
            return float(r0 + (z_point - z0) * (r1 - r0) / (z1 - z0))
        return 0.0

    def _sample_profile_boundary(
        self,
        z_profile: np.ndarray,
        r_profile: np.ndarray,
        rounding_radii: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        if z_profile.size < 2 or r_profile.size < 2:
            return z_profile, r_profile

        if not np.any(np.abs(rounding_radii) > 0.0):
            return z_profile, r_profile

        z_min = float(z_profile[0])
        z_max = float(z_profile[-1])
        total_span = z_max - z_min
        if total_span <= 0.0:
            return z_profile, r_profile

        positive_spans = np.diff(z_profile)
        positive_spans = positive_spans[positive_spans > 0.0]
        min_positive_span = float(np.min(positive_spans)) if positive_spans.size else total_span
        epsilon = max(min(total_span, min_positive_span) * 1.0e-4, 1.0e-9)

        sample_positions = set(np.linspace(z_min, z_max, max(220, 48 * z_profile.size)))
        for z_value, radius in zip(z_profile, rounding_radii):
            sample_positions.add(float(z_value))
            radius_value = abs(float(radius))
            if radius_value > 0.0:
                sample_positions.add(max(z_min, float(z_value) - radius_value))
                sample_positions.add(min(z_max, float(z_value) + radius_value))

        step_radii: dict[float, tuple[float, float]] = {}
        for index in range(z_profile.size - 1):
            if abs(float(z_profile[index + 1] - z_profile[index])) > 1.0e-15:
                continue
            step_z = float(z_profile[index])
            sample_positions.add(max(z_min, step_z - epsilon))
            sample_positions.add(min(z_max, step_z + epsilon))
            # Preserve explicit repeated-z profile edges as-authored. Rounded segments on either side
            # may start from the step endpoints, but they should not rewrite the vertical step itself.
            step_radii[step_z] = (float(r_profile[index]), float(r_profile[index + 1]))

        sampled_z: list[float] = []
        sampled_r: list[float] = []
        for z_value in sorted(sample_positions):
            if z_value in step_radii:
                sampled_z.append(z_value)
                sampled_r.append(step_radii[z_value][0])
                sampled_z.append(z_value)
                sampled_r.append(step_radii[z_value][1])
                continue
            sampled_z.append(float(z_value))
            sampled_r.append(self._boundary_radius_at_z(float(z_value), z_profile, r_profile, rounding_radii))

        return np.asarray(sampled_z, dtype=float), np.asarray(sampled_r, dtype=float)

    def plot_solid_profile(self, solid: dict[str, Any]) -> str:
        z_values, r_values, rounding_values = self._solid_profile_arrays(solid)
        if z_values.size < 2 or r_values.size < 2:
            self.plot_placeholder("The selected solid has no profile points.")
            return "No profile points available"

        display_z, display_r = self._sample_profile_boundary(z_values, r_values, rounding_values)

        self.figure.clear()
        axis = self.figure.add_subplot(111)
        solid_name = str(solid.get("name", "Selected solid"))

        axis.plot(display_z * 1000.0, display_r * 1000.0, linewidth=2.0)
        axis.plot(display_z * 1000.0, -display_r * 1000.0, linewidth=1.2, linestyle="--", alpha=0.7)
        axis.scatter(z_values * 1000.0, r_values * 1000.0, s=18)
        axis.scatter(z_values * 1000.0, -r_values * 1000.0, s=18)

        axis.axhline(0.0, color="#52606d", linewidth=0.8, alpha=0.6)
        axis.set_title(f"{solid_name} profile")
        axis.set_xlabel("z [mm]")
        axis.set_ylabel("radius [mm]")
        axis.grid(alpha=0.25)
        self._redraw()
        return f"Profile view for {solid_name}"

    def _solid_aperture_centers(self, solid: dict[str, Any]) -> np.ndarray:
        aperture_pattern = solid.get("aperturePattern")
        if not isinstance(aperture_pattern, dict):
            return np.asarray([[0.0, 0.0]], dtype=float)

        layout = str(aperture_pattern.get("layout", "single"))
        x_offset = float(aperture_pattern.get("xOffsetMeters", 0.0))
        y_offset = float(aperture_pattern.get("yOffsetMeters", 0.0))
        if layout == "single":
            return np.asarray([[x_offset, y_offset]], dtype=float)

        count_x = int(aperture_pattern.get("countX", 1))
        count_y = int(aperture_pattern.get("countY", 1))
        pitch_x = float(aperture_pattern.get("pitchXMeters", 0.0))
        pitch_y = float(aperture_pattern.get("pitchYMeters", 0.0))
        if count_x < 1 or count_y < 1 or pitch_x <= 0.0 or pitch_y <= 0.0:
            return np.asarray([[x_offset, y_offset]], dtype=float)

        row_shift = float(aperture_pattern.get("rowShiftXMeters", 0.0)) if layout == "staggered-grid" else 0.0
        x_extent = count_x // 2
        y_extent = count_y // 2
        centers: list[tuple[float, float]] = []
        for row_index in range(-y_extent, y_extent + 1):
            row_y = y_offset + row_index * pitch_y
            row_x_offset = 0.0
            if layout == "staggered-grid":
                row_x_offset = -row_shift if row_index % 2 == 0 else row_shift
            for column_index in range(-x_extent, x_extent + 1):
                centers.append((x_offset + row_x_offset + column_index * pitch_x, row_y))

        return np.asarray(centers, dtype=float)

    def _pattern_parameters(self, solid: dict[str, Any]) -> dict[str, float | str | bool]:
        aperture_pattern = solid.get("aperturePattern")
        if not isinstance(aperture_pattern, dict):
            return {
                "layout": "single",
                "count_x": 1,
                "count_y": 1,
                "pitch_x": 0.0,
                "pitch_y": 0.0,
                "margin": 0.0,
                "x_offset": 0.0,
                "y_offset": 0.0,
                "row_shift": 0.0,
                "outside_pattern_is_solid": True,
            }

        return {
            "layout": str(aperture_pattern.get("layout", "single")),
            "count_x": int(aperture_pattern.get("countX", 1) or 1),
            "count_y": int(aperture_pattern.get("countY", 1) or 1),
            "pitch_x": float(aperture_pattern.get("pitchXMeters", 0.0) or 0.0),
            "pitch_y": float(aperture_pattern.get("pitchYMeters", 0.0) or 0.0),
            "margin": float(aperture_pattern.get("marginMeters", 0.0) or 0.0),
            "x_offset": float(aperture_pattern.get("xOffsetMeters", 0.0) or 0.0),
            "y_offset": float(aperture_pattern.get("yOffsetMeters", 0.0) or 0.0),
            "row_shift": float(aperture_pattern.get("rowShiftXMeters", 0.0) or 0.0),
            "outside_pattern_is_solid": bool(aperture_pattern.get("outsidePatternIsSolid", True)),
        }

    def _solid_render_limits(
        self,
        solid: dict[str, Any],
        domain: dict[str, Any] | None,
        z_values: np.ndarray,
        r_values: np.ndarray,
    ) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]] | None:
        if z_values.size < 2 or r_values.size < 2:
            return None

        domain_x_limits, domain_y_limits, domain_z_limits = self._domain_clip_limits(domain)
        pattern = self._pattern_parameters(solid)
        layout = str(pattern["layout"])
        outside_pattern_is_solid = bool(pattern["outside_pattern_is_solid"])
        max_radius = float(np.max(np.abs(r_values))) if r_values.size else 0.0

        z_limits = (float(z_values[0]), float(z_values[-1]))
        if domain_z_limits is not None:
            z_limits = (max(z_limits[0], domain_z_limits[0]), min(z_limits[1], domain_z_limits[1]))
        if z_limits[1] <= z_limits[0]:
            return None

        if domain_x_limits is not None and (layout == "single" or outside_pattern_is_solid):
            x_limits = domain_x_limits
        elif layout != "single":
            half_x = 0.5 * float(pattern["count_x"]) * float(pattern["pitch_x"])
            margin = float(pattern["margin"])
            x_offset = float(pattern["x_offset"])
            row_shift = abs(float(pattern["row_shift"])) if layout == "staggered-grid" else 0.0
            x_limits = (x_offset - half_x - margin - row_shift, x_offset + half_x + margin + row_shift)
            if domain_x_limits is not None:
                x_limits = (max(x_limits[0], domain_x_limits[0]), min(x_limits[1], domain_x_limits[1]))
        else:
            x_offset = float(pattern["x_offset"])
            fallback_span = max(2.5 * max_radius, 1.0e-3)
            x_limits = (x_offset - fallback_span, x_offset + fallback_span)
            if domain_x_limits is not None:
                x_limits = (max(x_limits[0], domain_x_limits[0]), min(x_limits[1], domain_x_limits[1]))

        if domain_y_limits is not None and (layout == "single" or outside_pattern_is_solid):
            y_limits = domain_y_limits
        elif layout != "single":
            half_y = 0.5 * float(pattern["count_y"]) * float(pattern["pitch_y"])
            margin = float(pattern["margin"])
            y_offset = float(pattern["y_offset"])
            y_limits = (y_offset - half_y - margin, y_offset + half_y + margin)
            if domain_y_limits is not None:
                y_limits = (max(y_limits[0], domain_y_limits[0]), min(y_limits[1], domain_y_limits[1]))
        else:
            y_offset = float(pattern["y_offset"])
            fallback_span = max(2.5 * max_radius, 1.0e-3)
            y_limits = (y_offset - fallback_span, y_offset + fallback_span)
            if domain_y_limits is not None:
                y_limits = (max(y_limits[0], domain_y_limits[0]), min(y_limits[1], domain_y_limits[1]))

        if x_limits[1] <= x_limits[0] or y_limits[1] <= y_limits[0]:
            return None
        return x_limits, y_limits, z_limits

    def _solid_voxel_resolution(
        self,
        solid: dict[str, Any],
        x_limits: tuple[float, float],
        y_limits: tuple[float, float],
        z_limits: tuple[float, float],
        z_values: np.ndarray,
        r_values: np.ndarray,
    ) -> tuple[int, int, int]:
        pattern = self._pattern_parameters(solid)
        max_radius = float(np.max(np.abs(r_values))) if r_values.size else 1.0e-3
        x_span = max(x_limits[1] - x_limits[0], 1.0e-6)
        y_span = max(y_limits[1] - y_limits[0], 1.0e-6)
        z_span = max(z_limits[1] - z_limits[0], 1.0e-6)

        positive_spans = np.diff(z_values)
        positive_spans = positive_spans[positive_spans > 1.0e-9]
        min_z_feature = float(np.min(positive_spans)) if positive_spans.size else z_span
        x_feature = max_radius / 4.0
        y_feature = max_radius / 4.0
        if float(pattern["pitch_x"]) > 0.0:
            x_feature = min(x_feature, float(pattern["pitch_x"]) / 6.0)
        if float(pattern["pitch_y"]) > 0.0:
            y_feature = min(y_feature, float(pattern["pitch_y"]) / 6.0)
        z_feature = min(min_z_feature / 6.0, max_radius / 4.0)

        x_feature = max(x_feature, x_span / 72.0, 2.5e-4)
        y_feature = max(y_feature, y_span / 72.0, 2.5e-4)
        z_feature = max(z_feature, z_span / 96.0, 2.0e-4)

        x_count = int(np.clip(np.ceil(x_span / x_feature), 10, 72))
        y_count = int(np.clip(np.ceil(y_span / y_feature), 10, 72))
        z_count = int(np.clip(np.ceil(z_span / z_feature), 8, 96))
        return x_count, y_count, z_count

    def _sample_solid_volume(
        self,
        solid: dict[str, Any],
        domain: dict[str, Any] | None,
        z_values: np.ndarray,
        r_values: np.ndarray,
        rounding_values: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
        render_limits = self._solid_render_limits(solid, domain, z_values, r_values)
        if render_limits is None:
            return None

        x_limits, y_limits, z_limits = render_limits
        x_count, y_count, z_count = self._solid_voxel_resolution(solid, x_limits, y_limits, z_limits, z_values, r_values)
        x_edges = np.linspace(x_limits[0], x_limits[1], x_count + 1)
        y_edges = np.linspace(y_limits[0], y_limits[1], y_count + 1)
        z_edges = np.linspace(z_limits[0], z_limits[1], z_count + 1)
        x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
        y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
        z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])

        x_grid, y_grid = np.meshgrid(x_centers, y_centers, indexing="ij")
        filled = np.zeros((x_count, y_count, z_count), dtype=bool)
        pattern = self._pattern_parameters(solid)
        layout = str(pattern["layout"])
        x_offset = float(pattern["x_offset"])
        y_offset = float(pattern["y_offset"])
        row_shift = float(pattern["row_shift"])
        count_x = int(pattern["count_x"])
        count_y = int(pattern["count_y"])
        pitch_x = float(pattern["pitch_x"])
        pitch_y = float(pattern["pitch_y"])
        margin = float(pattern["margin"])
        outside_pattern_is_solid = bool(pattern["outside_pattern_is_solid"])
        half_x = 0.5 * count_x * pitch_x
        half_y = 0.5 * count_y * pitch_y

        for z_index, z_center in enumerate(z_centers):
            boundary_radius = self._boundary_radius_at_z(float(z_center), z_values, r_values, rounding_values)
            if boundary_radius <= 0.0 and (z_center < z_values[0] or z_center > z_values[-1]):
                continue

            local_x = x_grid - x_offset
            local_y = y_grid - y_offset
            if layout == "staggered-grid" and pitch_y > 0.0:
                rows = np.floor((local_y + pitch_y / 2.0) / pitch_y).astype(int)
                local_x = local_x + np.where(rows % 2 == 0, row_shift, -row_shift)

            if layout != "single" and pitch_x > 0.0 and pitch_y > 0.0:
                outside_mask = (np.abs(local_x) > half_x + margin) | (np.abs(local_y) > half_y + margin)
                fold_mask = (np.abs(local_x) < half_x) & (np.abs(local_y) < half_y)
                local_x_folded = np.array(local_x, copy=True)
                local_y_folded = np.array(local_y, copy=True)
                local_x_folded[fold_mask] -= np.round(local_x_folded[fold_mask] / pitch_x) * pitch_x
                local_y_folded[fold_mask] -= np.round(local_y_folded[fold_mask] / pitch_y) * pitch_y
                radius = np.hypot(local_x_folded, local_y_folded)
                solid_mask = radius >= boundary_radius
                solid_mask[outside_mask] = outside_pattern_is_solid
                filled[:, :, z_index] = solid_mask
            else:
                filled[:, :, z_index] = np.hypot(local_x, local_y) >= boundary_radius

        if not np.any(filled):
            return None

        return x_edges, y_edges, z_edges, filled

    def _domain_clip_limits(
        self,
        domain: dict[str, Any] | None,
    ) -> tuple[tuple[float, float] | None, tuple[float, float] | None, tuple[float, float] | None]:
        if not isinstance(domain, dict):
            return None, None, None

        x_size = float(domain.get("xSizeMeters", 0.0) or 0.0)
        y_size = float(domain.get("ySizeMeters", 0.0) or 0.0)
        z_size = float(domain.get("zSizeMeters", 0.0) or 0.0)
        x_limits = (-0.5 * x_size, 0.5 * x_size) if x_size > 0.0 else None
        y_limits = (-0.5 * y_size, 0.5 * y_size) if y_size > 0.0 else None
        z_limits = (0.0, z_size) if z_size > 0.0 else None
        return x_limits, y_limits, z_limits

    def _clip_segment_to_box(
        self,
        start_point: np.ndarray,
        end_point: np.ndarray,
        clip_limits: tuple[tuple[float, float] | None, tuple[float, float] | None, tuple[float, float] | None],
    ) -> tuple[np.ndarray, np.ndarray] | None:
        t_min = 0.0
        t_max = 1.0
        delta = end_point - start_point

        for axis_index, axis_limits in enumerate(clip_limits):
            if axis_limits is None:
                continue
            lower, upper = axis_limits
            axis_delta = float(delta[axis_index])
            axis_start = float(start_point[axis_index])
            if abs(axis_delta) <= 1.0e-15:
                if axis_start < lower or axis_start > upper:
                    return None
                continue

            entry_t = (lower - axis_start) / axis_delta
            exit_t = (upper - axis_start) / axis_delta
            if entry_t > exit_t:
                entry_t, exit_t = exit_t, entry_t
            t_min = max(t_min, entry_t)
            t_max = min(t_max, exit_t)
            if t_min > t_max + 1.0e-12:
                return None

        return start_point + t_min * delta, start_point + t_max * delta

    def _clip_polyline_to_box(
        self,
        polyline: np.ndarray,
        clip_limits: tuple[tuple[float, float] | None, tuple[float, float] | None, tuple[float, float] | None],
    ) -> list[np.ndarray]:
        if polyline.shape[0] < 2:
            return []
        if all(axis_limits is None for axis_limits in clip_limits):
            return [polyline]

        clipped_polylines: list[np.ndarray] = []
        current_points: list[np.ndarray] = []
        tolerance = 1.0e-12
        for point_index in range(polyline.shape[0] - 1):
            clipped_segment = self._clip_segment_to_box(
                polyline[point_index],
                polyline[point_index + 1],
                clip_limits,
            )
            if clipped_segment is None:
                if len(current_points) >= 2:
                    clipped_polylines.append(np.vstack(current_points))
                current_points = []
                continue

            clipped_start, clipped_end = clipped_segment
            if not current_points:
                current_points = [clipped_start, clipped_end]
                continue

            if float(np.linalg.norm(current_points[-1] - clipped_start)) > tolerance:
                if len(current_points) >= 2:
                    clipped_polylines.append(np.vstack(current_points))
                current_points = [clipped_start, clipped_end]
                continue

            current_points.append(clipped_end)

        if len(current_points) >= 2:
            clipped_polylines.append(np.vstack(current_points))
        return clipped_polylines

    def _segment_plane_intersection(
        self,
        start_point: np.ndarray,
        end_point: np.ndarray,
        axis_index: int,
        boundary: float,
    ) -> np.ndarray:
        delta = float(end_point[axis_index] - start_point[axis_index])
        if abs(delta) <= 1.0e-15:
            return np.asarray(start_point, dtype=float)

        fraction = (boundary - float(start_point[axis_index])) / delta
        fraction = max(0.0, min(1.0, fraction))
        return np.asarray(start_point + fraction * (end_point - start_point), dtype=float)

    def _clip_polygon_to_halfspace(
        self,
        polygon: np.ndarray,
        axis_index: int,
        boundary: float,
        keep_greater: bool,
    ) -> np.ndarray:
        if polygon.shape[0] < 3:
            return np.empty((0, 3), dtype=float)

        clipped: list[np.ndarray] = []
        tolerance = 1.0e-12
        previous_point = polygon[-1]
        previous_inside = (
            float(previous_point[axis_index]) >= boundary - tolerance
            if keep_greater
            else float(previous_point[axis_index]) <= boundary + tolerance
        )

        for current_point in polygon:
            current_inside = (
                float(current_point[axis_index]) >= boundary - tolerance
                if keep_greater
                else float(current_point[axis_index]) <= boundary + tolerance
            )
            if current_inside:
                if not previous_inside:
                    clipped.append(
                        self._segment_plane_intersection(previous_point, current_point, axis_index, boundary)
                    )
                clipped.append(np.asarray(current_point, dtype=float))
            elif previous_inside:
                clipped.append(self._segment_plane_intersection(previous_point, current_point, axis_index, boundary))

            previous_point = current_point
            previous_inside = current_inside

        if len(clipped) < 3:
            return np.empty((0, 3), dtype=float)
        return np.asarray(clipped, dtype=float)

    def _clip_polygon_to_box(
        self,
        polygon: np.ndarray,
        clip_limits: tuple[tuple[float, float] | None, tuple[float, float] | None, tuple[float, float] | None],
    ) -> np.ndarray:
        clipped_polygon = np.asarray(polygon, dtype=float)
        for axis_index, axis_limits in enumerate(clip_limits):
            if axis_limits is None:
                continue
            lower, upper = axis_limits
            clipped_polygon = self._clip_polygon_to_halfspace(clipped_polygon, axis_index, lower, keep_greater=True)
            if clipped_polygon.shape[0] < 3:
                return np.empty((0, 3), dtype=float)
            clipped_polygon = self._clip_polygon_to_halfspace(clipped_polygon, axis_index, upper, keep_greater=False)
            if clipped_polygon.shape[0] < 3:
                return np.empty((0, 3), dtype=float)
        return clipped_polygon

    def _sanitize_polygon(self, polygon: np.ndarray) -> np.ndarray:
        if polygon.shape[0] < 3:
            return np.empty((0, 3), dtype=float)

        sanitized_points: list[np.ndarray] = []
        tolerance = 1.0e-12
        for point in polygon:
            point_array = np.asarray(point, dtype=float)
            if not np.all(np.isfinite(point_array)):
                continue
            if sanitized_points and float(np.linalg.norm(point_array - sanitized_points[-1])) <= tolerance:
                continue
            sanitized_points.append(point_array)

        if len(sanitized_points) >= 2 and float(np.linalg.norm(sanitized_points[0] - sanitized_points[-1])) <= tolerance:
            sanitized_points.pop()
        if len(sanitized_points) < 3:
            return np.empty((0, 3), dtype=float)

        area_vector = np.zeros(3, dtype=float)
        for index in range(len(sanitized_points)):
            current_point = sanitized_points[index]
            next_point = sanitized_points[(index + 1) % len(sanitized_points)]
            area_vector += np.cross(current_point, next_point)
        if float(np.linalg.norm(area_vector)) <= 1.0e-14:
            return np.empty((0, 3), dtype=float)

        return np.asarray(sanitized_points, dtype=float)

    def plot_solid_3d(
        self,
        solid: dict[str, Any],
        domain: dict[str, Any] | None = None,
        show_solid: bool = False,
        view_direction: str | None = None,
    ) -> str:
        solid_name = str(solid.get("name", "Selected solid"))
        if not self._plot_solids_3d(
            [solid],
            title=f"{solid_name} 3D preview",
            domain=domain,
            show_solid=show_solid,
            view_direction=view_direction,
        ):
            self.plot_placeholder("The selected solid has no profile points to preview.")
            return "No solid geometry available"
        if show_solid:
            return f"Solid 3D view for {solid_name}"
        return f"Wireframe 3D view for {solid_name}"

    def plot_geometry_3d(
        self,
        solids: list[dict[str, Any]],
        domain: dict[str, Any] | None = None,
        show_solid: bool = False,
        view_direction: str | None = None,
    ) -> str:
        if not self._plot_solids_3d(
            solids,
            title="Structured geometry 3D preview",
            domain=domain,
            show_solid=show_solid,
            view_direction=view_direction,
        ):
            self.plot_placeholder("Add structured solids to preview the authored geometry.")
            return "No structured geometry available"
        if show_solid:
            return f"Solid 3D view for {len(solids)} solids"
        return f"Wireframe 3D view for {len(solids)} solids"

    def plot_geometry_sources_3d(
        self,
        solids: list[dict[str, Any]],
        sources: Any,
        domain: dict[str, Any] | None = None,
        show_solid: bool = True,
        view_direction: str | None = None,
    ) -> str:
        if not self._plot_solids_3d(
            solids,
            title="Generated geometry and sources",
            domain=domain,
            show_solid=show_solid,
            view_direction=view_direction,
            include_domain_x=True,
            include_domain_y=True,
            include_domain_z=True,
        ):
            self.plot_placeholder("No generated geometry is available for the selected case.")
            return "No generated geometry available"

        axis = self.figure.axes[0]
        source_count = self._plot_source_overlays_3d(axis, sources)
        self._redraw()
        if source_count == 0:
            return f"Displayed {len(solids)} solids in 3D; explicit source markers are unavailable"
        return f"Displayed {len(solids)} solids and {source_count} source marker(s) in 3D"

    def plot_geometry_sources_section(
        self,
        solids: list[dict[str, Any]],
        sources: Any,
        domain: dict[str, Any] | None = None,
    ) -> str:
        if not solids:
            self.plot_placeholder("No generated geometry is available for the selected case.")
            return "No generated geometry available"

        self.figure.clear()
        axis = self.figure.add_subplot(111)
        plotted_any = False
        geometry_legend_added = False
        for solid in solids:
            z_values, r_values, rounding_values = self._solid_profile_arrays(solid)
            if z_values.size < 2 or r_values.size < 2:
                continue
            display_z, display_r = self._sample_profile_boundary(z_values, r_values, rounding_values)
            display_z_mm = display_z * 1000.0
            display_r_mm = display_r * 1000.0
            polygon_z = np.concatenate([display_z_mm, display_z_mm[::-1]])
            polygon_r = np.concatenate([display_r_mm, -display_r_mm[::-1]])
            axis.fill(
                polygon_z,
                polygon_r,
                facecolor=GEOMETRY_FILL_COLOR,
                edgecolor=GEOMETRY_EDGE_COLOR,
                linewidth=0.9,
                alpha=0.75,
            )
            axis.plot(display_z_mm, display_r_mm, color=GEOMETRY_EDGE_COLOR, linewidth=1.0)
            axis.plot(display_z_mm, -display_r_mm, color=GEOMETRY_EDGE_COLOR, linewidth=1.0)
            if not geometry_legend_added:
                axis.plot([], [], color=GEOMETRY_EDGE_COLOR, linewidth=1.0, label="Geometry solids")
                geometry_legend_added = True
            plotted_any = True

        if not plotted_any:
            self.plot_placeholder("No generated geometry is available for the selected case.")
            return "No generated geometry available"

        source_count = self._plot_source_overlays_section(axis, sources)
        if isinstance(domain, dict):
            z_start = float(domain.get("zStartMeters", 0.0) or 0.0)
            z_size = float(domain.get("zSizeMeters", 0.0) or 0.0)
            x_size = float(domain.get("xSizeMeters", 0.0) or 0.0)
            if z_size > 0.0:
                axis.set_xlim(z_start * 1000.0, (z_start + z_size) * 1000.0)
            if x_size > 0.0:
                axis.set_ylim(-0.5 * x_size * 1000.0, 0.5 * x_size * 1000.0)

        axis.axhline(0.0, color="#52606d", linewidth=0.8, alpha=0.6)
        axis.set_title("Generated geometry and source section (z-x)")
        axis.set_xlabel("z [mm]")
        axis.set_ylabel("x / radial envelope [mm]")
        axis.grid(alpha=0.25)
        _set_equal_spatial_aspect(axis)
        axis.legend(loc="upper left", fontsize=8)
        self._redraw()
        if source_count == 0:
            return f"Displayed section view for {len(solids)} solids; explicit source markers are unavailable"
        return f"Displayed section view for {len(solids)} solids and {source_count} source marker(s)"

    def _runtime_sources(self, sources: Any) -> list[dict[str, Any]]:
        if not isinstance(sources, list):
            return []
        return [source for source in sources if isinstance(source, dict)]

    def _plot_source_overlays_3d(self, axis: Any, sources: Any) -> int:
        source_count = 0
        for source_index, source in enumerate(self._runtime_sources(sources), start=1):
            uniform = source.get("uniform") if isinstance(source.get("uniform"), dict) else {}
            center = uniform.get("centerMeters")
            direction = uniform.get("mainDirection")
            if not isinstance(center, list) or len(center) != 3:
                continue
            center_mm = np.asarray(center, dtype=float) * 1000.0
            axis.scatter(center_mm[0], center_mm[1], center_mm[2], color="#c24d2c", s=42, marker="o")
            if isinstance(direction, list) and len(direction) == 3:
                direction_vector = np.asarray(direction, dtype=float)
                norm = float(np.linalg.norm(direction_vector))
                if norm > 1.0e-12:
                    scale = max(
                        float(uniform.get("widthMeters", 0.0) or 0.0),
                        float(uniform.get("heightMeters", 0.0) or 0.0),
                        1.0e-3,
                    ) * 1000.0
                    direction_unit = direction_vector / norm
                    axis.quiver(
                        center_mm[0],
                        center_mm[1],
                        center_mm[2],
                        direction_unit[0] * scale,
                        direction_unit[1] * scale,
                        direction_unit[2] * scale,
                        color="#c24d2c",
                        linewidth=1.0,
                        arrow_length_ratio=0.25,
                    )
            axis.text(
                center_mm[0],
                center_mm[1],
                center_mm[2],
                str(source.get("name") or source.get("id") or f"Source {source_index}"),
                fontsize=8,
            )
            source_count += 1
        return source_count

    def _plot_source_overlays_section(self, axis: Any, sources: Any) -> int:
        source_count = 0
        for source_index, source in enumerate(self._runtime_sources(sources), start=1):
            uniform = source.get("uniform") if isinstance(source.get("uniform"), dict) else {}
            center = uniform.get("centerMeters")
            if not isinstance(center, list) or len(center) != 3:
                continue
            center_x = float(center[0]) * 1000.0
            center_z = float(center[2]) * 1000.0
            axis.scatter(center_z, center_x, color="#c24d2c", s=42, marker="x")
            axis.text(
                center_z,
                center_x,
                f" {source.get('name') or source.get('id') or f'Source {source_index}'}",
                fontsize=8,
                va="bottom",
            )
            source_count += 1
        return source_count

    def _plot_solids_3d(
        self,
        solids: list[dict[str, Any]],
        title: str,
        domain: dict[str, Any] | None = None,
        show_solid: bool = False,
        view_direction: str | None = None,
        include_domain_x: bool = False,
        include_domain_y: bool = False,
        include_domain_z: bool = False,
    ) -> bool:
        self.figure.clear()
        axis = self.figure.add_subplot(111, projection="3d")
        axis.set_proj_type("ortho")
        theta = np.linspace(0.0, 2.0 * np.pi, 28)
        palette = cm.tab10(np.linspace(0.0, 1.0, max(len(solids), 1)))
        plotted_points: list[np.ndarray] = []
        plotted_any = False
        clip_limits = self._domain_clip_limits(domain)

        for solid_index, solid in enumerate(solids):
            color = palette[solid_index % len(palette)]
            solid_name = str(solid.get("name", f"Solid {solid_index + 1}"))
            z_values, r_values, rounding_values = self._solid_profile_arrays(solid)
            if z_values.size < 2 or r_values.size < 2:
                continue
            display_z, display_r = self._sample_profile_boundary(z_values, r_values, rounding_values)
            solid_segments: list[np.ndarray] = []
            if show_solid:
                for center_x, center_y in self._solid_aperture_centers(solid):
                    solid_segments.extend(
                        self._plot_revolved_surface(
                            axis,
                            display_z,
                            display_r,
                            theta,
                            color,
                            center_x,
                            center_y,
                            clip_limits,
                        )
                    )
            else:
                for center_x, center_y in self._solid_aperture_centers(solid):
                    solid_segments.extend(
                        self._plot_revolved_section(
                            axis,
                            display_z,
                            display_r,
                            theta,
                            color,
                            center_x,
                            center_y,
                            clip_limits,
                        )
                    )

            if solid_segments:
                plotted_any = True
                plotted_points.extend(solid_segments)
                axis.plot([], [], [], color=color, label=solid_name)

        if not plotted_any:
            return False

        point_cloud = np.vstack(plotted_points)
        self._set_3d_bounds(
            axis,
            point_cloud,
            domain=domain,
            include_domain_x=include_domain_x,
            include_domain_y=include_domain_y,
            include_domain_z=include_domain_z,
        )
        axis.set_title(title)
        axis.set_xlabel("x [mm]")
        axis.set_ylabel("y [mm]")
        axis.set_zlabel("z [mm]")
        self._apply_3d_view(axis, view_direction)
        axis.legend(loc="upper left", fontsize=8)
        self._redraw()
        return True

        view_angles = {
            "+x": (0.0, 0.0),
            "-x": (0.0, 180.0),
            "+y": (0.0, 90.0),
            "-y": (0.0, -90.0),
            "+z": (90.0, -90.0),
            "-z": (-90.0, -90.0),
        }
        elev, azim = view_angles.get(view_direction, (None, None))
        if elev is None or azim is None:
            return
        axis.view_init(elev=elev, azim=azim)

    def _plot_revolved_section(
        self,
        axis: Any,
        z_values: np.ndarray,
        r_values: np.ndarray,
        theta: np.ndarray,
        color: Any,
        center_x: float = 0.0,
        center_y: float = 0.0,
        clip_limits: tuple[tuple[float, float] | None, tuple[float, float] | None, tuple[float, float] | None] = (None, None, None),
    ) -> list[np.ndarray]:
        plotted_segments: list[np.ndarray] = []
        for theta_value in theta[::4]:
            polyline = np.column_stack(
                (
                    center_x + r_values * np.cos(theta_value),
                    center_y + r_values * np.sin(theta_value),
                    z_values,
                )
            )
            for clipped_polyline in self._clip_polyline_to_box(polyline, clip_limits):
                axis.plot(
                    clipped_polyline[:, 0] * 1000.0,
                    clipped_polyline[:, 1] * 1000.0,
                    clipped_polyline[:, 2] * 1000.0,
                    color=color,
                    linewidth=0.9,
                    alpha=0.7,
                )
                plotted_segments.append(clipped_polyline)

        ring_indices = np.unique(np.linspace(0, len(z_values) - 1, min(len(z_values), 16), dtype=int))
        for ring_index in ring_indices:
            z_value = z_values[ring_index]
            r_value = r_values[ring_index]
            ring_polyline = np.column_stack(
                (
                    center_x + r_value * np.cos(theta),
                    center_y + r_value * np.sin(theta),
                    np.full_like(theta, z_value),
                )
            )
            for clipped_polyline in self._clip_polyline_to_box(ring_polyline, clip_limits):
                axis.plot(
                    clipped_polyline[:, 0] * 1000.0,
                    clipped_polyline[:, 1] * 1000.0,
                    clipped_polyline[:, 2] * 1000.0,
                    color=color,
                    linewidth=0.7,
                    alpha=0.35,
                )
                plotted_segments.append(clipped_polyline)

        return plotted_segments

    def _plot_revolved_surface(
        self,
        axis: Any,
        z_values: np.ndarray,
        r_values: np.ndarray,
        theta: np.ndarray,
        color: Any,
        center_x: float = 0.0,
        center_y: float = 0.0,
        clip_limits: tuple[tuple[float, float] | None, tuple[float, float] | None, tuple[float, float] | None] = (None, None, None),
    ) -> list[np.ndarray]:
        x_grid = center_x + np.outer(r_values, np.cos(theta))
        y_grid = center_y + np.outer(r_values, np.sin(theta))
        z_grid = np.outer(z_values, np.ones_like(theta))
        clipped_polygons: list[np.ndarray] = []

        for z_index in range(len(z_values) - 1):
            for theta_index in range(len(theta) - 1):
                quad = np.asarray(
                    [
                        [x_grid[z_index, theta_index], y_grid[z_index, theta_index], z_grid[z_index, theta_index]],
                        [x_grid[z_index + 1, theta_index], y_grid[z_index + 1, theta_index], z_grid[z_index + 1, theta_index]],
                        [x_grid[z_index + 1, theta_index + 1], y_grid[z_index + 1, theta_index + 1], z_grid[z_index + 1, theta_index + 1]],
                        [x_grid[z_index, theta_index + 1], y_grid[z_index, theta_index + 1], z_grid[z_index, theta_index + 1]],
                    ],
                    dtype=float,
                )
                clipped_quad = self._clip_polygon_to_box(quad, clip_limits)
                sanitized_polygon = self._sanitize_polygon(clipped_quad)
                if sanitized_polygon.shape[0] >= 3:
                    clipped_polygons.append(sanitized_polygon)

        if clipped_polygons:
            color_values = np.asarray(color, dtype=float)
            face_color = tuple(float(component) for component in color_values[:3]) + (0.45,)
            axis.add_collection3d(
                Poly3DCollection(
                    [polygon * 1000.0 for polygon in clipped_polygons],
                    facecolors=face_color,
                    edgecolors="none",
                    linewidths=0.0,
                )
            )

        return clipped_polygons

    def _set_3d_bounds(
        self,
        axis: Any,
        point_cloud: np.ndarray,
        domain: dict[str, Any] | None = None,
        include_domain_x: bool = False,
        include_domain_y: bool = False,
        include_domain_z: bool = False,
    ) -> None:
        x_limits = None
        y_limits = None
        z_limits = None
        domain_x_limits, domain_y_limits, domain_z_limits = self._domain_clip_limits(domain)
        if include_domain_x and domain_x_limits is not None:
            x_limits = (domain_x_limits[0] * 1000.0, domain_x_limits[1] * 1000.0)
        if include_domain_y and domain_y_limits is not None:
            y_limits = (domain_y_limits[0] * 1000.0, domain_y_limits[1] * 1000.0)
        if include_domain_z and domain_z_limits is not None:
            z_limits = (domain_z_limits[0] * 1000.0, domain_z_limits[1] * 1000.0)

        _set_adaptive_3d_bounds(
            axis,
            point_cloud,
            scale=1000.0,
            x_limits=x_limits,
            y_limits=y_limits,
            z_limits=z_limits,
        )


class TrajectoryCanvas(FigureCanvasQTAgg):
    def __init__(self) -> None:
        self.figure = Figure(figsize=(8, 6), tight_layout=True)
        super().__init__(self.figure)
        configure_matplotlib_canvas(self)
        self.plot_placeholder("Materialize or run a case, then choose a VTK file to visualize.")

    def _draw_stats_box(self, axis: Any, text: str) -> None:
        axis.text(
            0.03,
            0.97,
            text,
            transform=axis.transAxes,
            ha="left",
            va="top",
            multialignment="left",
            linespacing=1.35,
            fontsize=9,
            color="#102a43",
            fontfamily="DejaVu Sans Mono",
            bbox={
                "boxstyle": "round,pad=0.45",
                "facecolor": "#fff8e8",
                "edgecolor": "#d9b86c",
                "linewidth": 0.9,
                "alpha": 0.96,
            },
        )

    def plot_placeholder(self, message: str) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.axis("off")
        axis.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
        self.draw_idle()

    def _read_vtk(self, vtk_path: Path) -> tuple[np.ndarray, list[np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray]]:
        with contextlib.redirect_stdout(io.StringIO()):
            return PLOT_UTILS.read_vtk_polydata(str(vtk_path))

    def plot_cross_section(
        self,
        vtk_path: Path,
        z_plane: float,
        trajectory_options: TrajectoryRenderOptions | None = None,
    ) -> str:
        points, lines, point_data, cell_data = self._read_vtk(vtk_path)
        options = _normalize_trajectory_options(trajectory_options)
        selected_particle_kinds = None
        if options.selected_particle_kinds is not None:
            selected_particle_kinds = {int(value) for value in options.selected_particle_kinds}

        crossings: list[tuple[float, float]] = []
        particle_status: list[float | None] = []
        particle_kinds: list[int | None] = []
        crossing_energies: list[float | None] = []

        for index, line_indices in enumerate(lines):
            particle_kind = _resolve_particle_kind_from_cell_data(cell_data, index, options.particle_family)
            if selected_particle_kinds is not None and particle_kind not in selected_particle_kinds:
                continue

            crossing = _trajectory_plane_crossing(points, line_indices, "z", z_plane)
            if crossing is None:
                continue
            segment_index, point_a, point_b, fraction = crossing

            crossing_energy: float | None = None
            energies = _trajectory_energy_values(point_data, line_indices)
            if energies is not None and segment_index + 1 < len(energies):
                interpolated_energy = _interpolate_scalar_value(
                    float(energies[segment_index]),
                    float(energies[segment_index + 1]),
                    fraction,
                )
                if math.isfinite(interpolated_energy):
                    crossing_energy = interpolated_energy

            if crossing_energy is not None:
                if options.energy_min_ev is not None and crossing_energy < float(options.energy_min_ev):
                    continue
                if options.energy_max_ev is not None and crossing_energy > float(options.energy_max_ev):
                    continue

            crossings.append(
                (
                    _interpolate_scalar_value(float(point_a[0]), float(point_b[0]), fraction),
                    _interpolate_scalar_value(float(point_a[1]), float(point_b[1]), fraction),
                )
            )
            particle_kinds.append(particle_kind)
            crossing_energies.append(crossing_energy)
            if "particle_status" in cell_data and index < len(cell_data["particle_status"]):
                particle_status.append(float(cell_data["particle_status"][index]))
            else:
                particle_status.append(None)

        if not crossings:
            if points.size == 0:
                self.plot_placeholder("The selected VTK file does not contain trajectory lines.")
                return "No trajectories available"
            self.plot_placeholder(
                f"No trajectory crossings found at z = {z_plane:.4f} m.\n"
                f"The file spans z = [{points[:, 2].min():.4f}, {points[:, 2].max():.4f}] m."
            )
            return "No crossings found"

        samples = np.array(crossings)
        x_values = samples[:, 0]
        y_values = samples[:, 1]
        x_mean = float(x_values.mean())
        y_mean = float(y_values.mean())
        x_std = float(x_values.std())
        y_std = float(y_values.std())

        self.figure.clear()
        scatter_axis = self.figure.add_subplot(121)
        density_axis = self.figure.add_subplot(122)

        color_mode = _trajectory_color_mode_key(options.color_mode)
        status_available = len(particle_status) == len(crossings) and all(value is not None for value in particle_status)
        energy_available = len(crossing_energies) == len(crossings) and all(value is not None for value in crossing_energies)
        species_used = {kind for kind in particle_kinds if kind is not None}

        if color_mode == TRAJECTORY_COLOR_MODE_ENERGY.lower() and energy_available:
            scatter = scatter_axis.scatter(
                x_values * 1000.0,
                y_values * 1000.0,
                c=np.asarray(crossing_energies, dtype=float),
                cmap="plasma",
                s=20,
                alpha=0.75,
            )
            self.figure.colorbar(scatter, ax=scatter_axis, label="Kinetic energy [eV]")
        elif color_mode == TRAJECTORY_COLOR_MODE_SPECIES.lower() and particle_kinds:
            scatter_axis.scatter(
                x_values * 1000.0,
                y_values * 1000.0,
                c=[
                    particle_color_from_export_id(kind) if kind is not None else DEFAULT_TRAJECTORY_COLOR
                    for kind in particle_kinds
                ],
                s=20,
                alpha=0.75,
            )
            if species_used:
                _add_species_legend(
                    scatter_axis,
                    species_used,
                    options.particle_family,
                    marker=True,
                )
                scatter_axis.legend(loc="upper right")
        elif color_mode == TRAJECTORY_COLOR_MODE_STATUS.lower() and status_available:
            scatter = scatter_axis.scatter(
                x_values * 1000.0,
                y_values * 1000.0,
                c=np.asarray(particle_status, dtype=float),
                cmap="viridis",
                s=20,
                alpha=0.75,
            )
            self.figure.colorbar(scatter, ax=scatter_axis, label="Particle status")
        else:
            scatter_axis.scatter(
                x_values * 1000.0,
                y_values * 1000.0,
                s=20,
                alpha=0.75,
                color="#c24d2c",
            )

        scatter_axis.axvline(x_mean * 1000.0, linestyle="--", color="#334e68", alpha=0.7)
        scatter_axis.axhline(y_mean * 1000.0, linestyle="--", color="#334e68", alpha=0.7)
        scatter_axis.set_title(f"Beam Slice at z = {z_plane * 1000.0:.1f} mm")
        scatter_axis.set_xlabel("x [mm]")
        scatter_axis.set_ylabel("y [mm]")
        scatter_axis.grid(alpha=0.25)
        scatter_axis.set_aspect("equal", adjustable="box")

        hist, x_edges, y_edges = np.histogram2d(x_values * 1000.0, y_values * 1000.0, bins=48)
        image = density_axis.imshow(
            hist.T,
            origin="lower",
            extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
            aspect="equal",
            cmap="magma",
            interpolation="bilinear",
        )
        self.figure.colorbar(image, ax=density_axis, label="Particle count")
        density_axis.set_title("Cross-section density")
        density_axis.set_xlabel("x [mm]")
        density_axis.set_ylabel("y [mm]")
        density_axis.grid(alpha=0.2, color="white")
        _set_equal_spatial_aspect(density_axis)

        stats_text = (
            f"N = {len(crossings)}\n"
            f"<x> = {x_mean * 1000.0:.3f} mm\n"
            f"<y> = {y_mean * 1000.0:.3f} mm\n"
            f"sigma_x = {x_std * 1000.0:.3f} mm\n"
            f"sigma_y = {y_std * 1000.0:.3f} mm"
        )
        self._draw_stats_box(density_axis, stats_text)

        self.draw_idle()
        return f"Plotted {len(crossings)} crossings at z = {z_plane:.4f} m"

    def plot_3d(
        self,
        vtk_path: Path,
        max_trajectories: int,
        trajectory_options: TrajectoryRenderOptions | None = None,
    ) -> str:
        points, lines, point_data, cell_data = self._read_vtk(vtk_path)
        options = _normalize_trajectory_options(trajectory_options)
        selected_particle_kinds = None
        if options.selected_particle_kinds is not None:
            selected_particle_kinds = {int(value) for value in options.selected_particle_kinds}

        available_trajectories: list[tuple[int, np.ndarray, int | None]] = []
        for index, line_indices in enumerate(lines):
            particle_kind = _resolve_particle_kind_from_cell_data(cell_data, index, options.particle_family)
            if selected_particle_kinds is not None and particle_kind not in selected_particle_kinds:
                continue
            available_trajectories.append((index, line_indices, particle_kind))

        count = min(max_trajectories, len(available_trajectories))
        if count == 0:
            self.plot_placeholder("No trajectories matched the selected filters.")
            return "No trajectories available"

        self.figure.clear()
        axis = self.figure.add_subplot(111, projection="3d")
        axis.set_proj_type("ortho")

        color_mode = _trajectory_color_mode_key(options.color_mode)
        status_values = np.asarray(cell_data.get("particle_status", []), dtype=float)
        status_max = max(float(np.max(status_values)), 1.0) if status_values.size else 1.0
        visible_clouds: list[np.ndarray] = []
        visible_trajectory_count = 0
        species_used: set[int] = set()

        if color_mode == TRAJECTORY_COLOR_MODE_ENERGY.lower():
            energy_segments: list[np.ndarray] = []
            energy_values: list[float] = []
            for _trajectory_index, line_indices, _particle_kind in available_trajectories[:count]:
                trajectory = np.asarray(points[line_indices], dtype=float)
                energies = _trajectory_energy_values(point_data, line_indices)
                if energies is None:
                    continue
                visible = _trajectory_visibility_mask(len(trajectory), energies, options)
                segments, segment_energies = _visible_segment_pairs(trajectory, visible, energies)
                if not segments:
                    continue
                energy_segments.extend([segment * 1000.0 for segment in segments])
                energy_values.extend(segment_energies)
                visible_clouds.extend(segments)
                visible_trajectory_count += 1

            if not energy_segments:
                self.plot_placeholder("No trajectory segments remained inside the selected energy range.")
                return "No trajectories available"

            norm = _build_normalize(min(energy_values), max(energy_values))
            collection = Line3DCollection(
                energy_segments,
                cmap="plasma",
                norm=norm,
                linewidths=0.85,
                alpha=0.72,
            )
            collection.set_array(np.asarray(energy_values, dtype=float))
            axis.add_collection(collection)
            self.figure.colorbar(collection, ax=axis, label="Kinetic energy [eV]")
        else:
            for trajectory_index, line_indices, particle_kind in available_trajectories[:count]:
                trajectory = np.asarray(points[line_indices], dtype=float)
                energies = _trajectory_energy_values(point_data, line_indices)
                visible = _trajectory_visibility_mask(len(trajectory), energies, options)
                visible_segments = _split_visible_polylines(trajectory, visible)
                if not visible_segments:
                    continue

                if color_mode == TRAJECTORY_COLOR_MODE_SPECIES.lower() and particle_kind is not None:
                    color = particle_color_from_export_id(particle_kind)
                    species_used.add(particle_kind)
                elif (
                    color_mode == TRAJECTORY_COLOR_MODE_STATUS.lower()
                    and "particle_status" in cell_data
                    and trajectory_index < len(cell_data["particle_status"])
                ):
                    mapped = float(cell_data["particle_status"][trajectory_index]) / status_max
                    color = cm.viridis(mapped)
                else:
                    color = DEFAULT_TRAJECTORY_COLOR

                for segment in visible_segments:
                    axis.plot(
                        segment[:, 0] * 1000.0,
                        segment[:, 1] * 1000.0,
                        segment[:, 2] * 1000.0,
                        color=color,
                        linewidth=0.8,
                        alpha=0.65,
                    )
                    visible_clouds.append(segment)
                visible_trajectory_count += 1

            if visible_trajectory_count == 0:
                self.plot_placeholder("No trajectory segments remained inside the selected filters.")
                return "No trajectories available"

            if color_mode == TRAJECTORY_COLOR_MODE_STATUS.lower() and status_values.size:
                status_mappable = cm.ScalarMappable(norm=_build_normalize(0.0, status_max), cmap="viridis")
                status_mappable.set_array(np.asarray([0.0, status_max], dtype=float))
                self.figure.colorbar(status_mappable, ax=axis, label="Particle status")
            if color_mode == TRAJECTORY_COLOR_MODE_SPECIES.lower() and species_used:
                _add_species_legend(
                    axis,
                    species_used,
                    options.particle_family,
                    marker=False,
                    three_dimensional=True,
                )
                axis.legend(loc="upper right")

        visible_point_cloud = (
            np.vstack([np.asarray(segment, dtype=float) for segment in visible_clouds])
            if visible_clouds
            else points
        )
        _set_adaptive_3d_bounds(axis, visible_point_cloud, scale=1000.0)
        axis.set_title(f"3D Trajectories ({visible_trajectory_count} shown)")
        axis.set_xlabel("x [mm]")
        axis.set_ylabel("y [mm]")
        axis.set_zlabel("z [mm]")
        self.draw_idle()
        return f"Plotted {visible_trajectory_count} trajectories in 3D"

    def plot_phase_space(
        self,
        vtk_path: Path,
        z_plane: float,
        transverse_axis: str,
        trajectory_options: TrajectoryRenderOptions | None = None,
    ) -> str:
        points, lines, point_data, cell_data = self._read_vtk(vtk_path)
        options = _normalize_trajectory_options(trajectory_options)
        selected_particle_kinds = None
        if options.selected_particle_kinds is not None:
            selected_particle_kinds = {int(value) for value in options.selected_particle_kinds}

        coordinates_mm: list[float] = []
        angles_mrad: list[float] = []
        particle_status: list[float | None] = []
        particle_kinds: list[int | None] = []
        crossing_energies: list[float | None] = []
        discarded_parallel = 0

        for index, line_indices in enumerate(lines):
            particle_kind = _resolve_particle_kind_from_cell_data(cell_data, index, options.particle_family)
            if selected_particle_kinds is not None and particle_kind not in selected_particle_kinds:
                continue

            crossing = _trajectory_plane_crossing(points, line_indices, "z", z_plane)
            if crossing is None:
                continue
            segment_index, point_a, point_b, fraction = crossing
            dz = float(point_b[2] - point_a[2])
            dx = float(point_b[0] - point_a[0])
            dy = float(point_b[1] - point_a[1])
            if abs(dz) < 1.0e-12:
                discarded_parallel += 1
                continue

            crossing_energy: float | None = None
            energies = _trajectory_energy_values(point_data, line_indices)
            if energies is not None and segment_index + 1 < len(energies):
                interpolated_energy = _interpolate_scalar_value(
                    float(energies[segment_index]),
                    float(energies[segment_index + 1]),
                    fraction,
                )
                if math.isfinite(interpolated_energy):
                    crossing_energy = interpolated_energy

            if crossing_energy is not None:
                if options.energy_min_ev is not None and crossing_energy < float(options.energy_min_ev):
                    continue
                if options.energy_max_ev is not None and crossing_energy > float(options.energy_max_ev):
                    continue

            x_crossing = _interpolate_scalar_value(float(point_a[0]), float(point_b[0]), fraction)
            y_crossing = _interpolate_scalar_value(float(point_a[1]), float(point_b[1]), fraction)
            if transverse_axis == "y":
                coordinates_mm.append(y_crossing * 1000.0)
                angles_mrad.append(float(np.arctan2(dy, dz)) * 1000.0)
            else:
                coordinates_mm.append(x_crossing * 1000.0)
                angles_mrad.append(float(np.arctan2(dx, dz)) * 1000.0)

            particle_kinds.append(particle_kind)
            crossing_energies.append(crossing_energy)
            if "particle_status" in cell_data and index < len(cell_data["particle_status"]):
                particle_status.append(float(cell_data["particle_status"][index]))
            else:
                particle_status.append(None)

        if not coordinates_mm:
            if points.size == 0:
                self.plot_placeholder("The selected VTK file does not contain trajectory lines.")
                return "No trajectories available"
            self.plot_placeholder(
                f"No phase-space samples found at z = {z_plane:.4f} m.\n"
                f"The file spans z = [{points[:, 2].min():.4f}, {points[:, 2].max():.4f}] m."
            )
            return "No phase-space samples found"

        coord_array = np.asarray(coordinates_mm, dtype=float)
        angle_array = np.asarray(angles_mrad, dtype=float)

        self.figure.clear()
        scatter_axis = self.figure.add_subplot(121)
        density_axis = self.figure.add_subplot(122)

        color_mode = _trajectory_color_mode_key(options.color_mode)
        status_available = len(particle_status) == len(coord_array) and all(value is not None for value in particle_status)
        energy_available = len(crossing_energies) == len(coord_array) and all(value is not None for value in crossing_energies)
        species_used = {kind for kind in particle_kinds if kind is not None}

        if color_mode == TRAJECTORY_COLOR_MODE_ENERGY.lower() and energy_available:
            scatter = scatter_axis.scatter(
                coord_array,
                angle_array,
                c=np.asarray(crossing_energies, dtype=float),
                cmap="plasma",
                s=16,
                marker=".",
                alpha=0.75,
                linewidths=0.0,
                edgecolors="none",
            )
            self.figure.colorbar(scatter, ax=scatter_axis, label="Kinetic energy [eV]")
        elif color_mode == TRAJECTORY_COLOR_MODE_SPECIES.lower() and particle_kinds:
            scatter_axis.scatter(
                coord_array,
                angle_array,
                c=[
                    particle_color_from_export_id(kind) if kind is not None else DEFAULT_TRAJECTORY_COLOR
                    for kind in particle_kinds
                ],
                s=16,
                marker=".",
                alpha=0.75,
                linewidths=0.0,
                edgecolors="none",
            )
            if species_used:
                _add_species_legend(
                    scatter_axis,
                    species_used,
                    options.particle_family,
                    marker=True,
                )
                scatter_axis.legend(loc="upper right")
        elif color_mode == TRAJECTORY_COLOR_MODE_STATUS.lower() and status_available:
            scatter = scatter_axis.scatter(
                coord_array,
                angle_array,
                c=np.asarray(particle_status, dtype=float),
                cmap="viridis",
                s=16,
                marker=".",
                alpha=0.75,
                linewidths=0.0,
                edgecolors="none",
            )
            self.figure.colorbar(scatter, ax=scatter_axis, label="Particle status")
        else:
            scatter_axis.scatter(
                coord_array,
                angle_array,
                color=DEFAULT_TRAJECTORY_COLOR,
                s=16,
                marker=".",
                alpha=0.75,
                linewidths=0.0,
                edgecolors="none",
            )

        coord_label = "x" if transverse_axis != "y" else "y"
        angle_label = "x'" if transverse_axis != "y" else "y'"
        scatter_axis.set_title(f"{coord_label}-{angle_label} at z = {z_plane * 1000.0:.1f} mm")
        scatter_axis.set_xlabel(f"{coord_label} [mm]")
        scatter_axis.set_ylabel(f"{angle_label} [mrad]")
        scatter_axis.grid(alpha=0.25)

        hist, x_edges, y_edges = np.histogram2d(coord_array, angle_array, bins=48)
        masked_hist = np.ma.masked_less_equal(hist.T, 0.0)
        density_cmap = cm.get_cmap("magma").copy()
        density_cmap.set_bad("#ffffff")
        image = density_axis.imshow(
            masked_hist,
            origin="lower",
            extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
            aspect="auto",
            cmap=density_cmap,
            interpolation="nearest",
        )
        self.figure.colorbar(image, ax=density_axis, label="Particle count")
        density_axis.set_title("Phase-space density")
        density_axis.set_xlabel(f"{coord_label} [mm]")
        density_axis.set_ylabel(f"{angle_label} [mrad]")
        density_axis.grid(alpha=0.2, color="white")

        mean_coord = float(coord_array.mean())
        mean_angle = float(angle_array.mean())
        std_coord = float(coord_array.std())
        std_angle = float(angle_array.std())
        stats_text = (
            f"N = {len(coord_array)}\n"
            f"<{coord_label}> = {mean_coord:.3f} mm\n"
            f"<{angle_label}> = {mean_angle:.3f} mrad\n"
            f"sigma_{coord_label} = {std_coord:.3f} mm\n"
            f"sigma_{angle_label} = {std_angle:.3f} mrad"
        )
        if discarded_parallel > 0:
            stats_text += f"\ndiscarded parallel: {discarded_parallel}"
        self._draw_stats_box(density_axis, stats_text)

        self.draw_idle()
        return f"Plotted {coord_label}-{angle_label} phase space with {len(coord_array)} samples"


class OutputCanvas(FigureCanvasQTAgg):
    def __init__(self) -> None:
        self.figure = Figure(figsize=(8, 6), tight_layout=True)
        super().__init__(self.figure)
        configure_matplotlib_canvas(self)
        self.plot_placeholder("Select an output artifact to visualize.")

    def plot_placeholder(self, message: str) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.axis("off")
        axis.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
        self.draw_idle()

    def plot_vtk_cross_section(
        self,
        vtk_path: Path,
        z_plane: float,
        trajectory_options: TrajectoryRenderOptions | None = None,
    ) -> str:
        return self._trajectory_canvas().plot_cross_section(vtk_path, z_plane, trajectory_options)

    def plot_vtk_3d(
        self,
        vtk_path: Path,
        max_trajectories: int,
        trajectory_options: TrajectoryRenderOptions | None = None,
    ) -> str:
        return self._trajectory_canvas().plot_3d(vtk_path, max_trajectories, trajectory_options)

    def plot_vtk_phase_space(
        self,
        vtk_path: Path,
        z_plane: float,
        transverse_axis: str,
        trajectory_options: TrajectoryRenderOptions | None = None,
    ) -> str:
        return self._trajectory_canvas().plot_phase_space(
            vtk_path,
            z_plane,
            transverse_axis,
            trajectory_options,
        )

    def resolve_trajectory_source_path(self, vtk_path: Path) -> Path | None:
        with contextlib.redirect_stdout(io.StringIO()):
            dataset_type = PLOT_UTILS.read_vtk_dataset_type(str(vtk_path))
        if dataset_type == "DATASET POLYDATA":
            return vtk_path
        if dataset_type == "DATASET STRUCTURED_POINTS":
            return self._resolve_field_overlay_path(vtk_path, "trajectories")
        return None

    def get_vtk_trajectory_metadata(
        self,
        vtk_path: Path,
        particle_family: str = "H",
    ) -> TrajectoryMetadata | None:
        source_path = self.resolve_trajectory_source_path(vtk_path)
        if source_path is None or not source_path.is_file() or source_path.stat().st_size == 0:
            return None

        points, lines, point_data, cell_data = self._trajectory_canvas()._read_vtk(source_path)
        particle_family = _normalize_particle_family(particle_family)
        available_particle_kinds = sorted(
            {
                particle_kind
                for index in range(len(lines))
                for particle_kind in [_resolve_particle_kind_from_cell_data(cell_data, index, particle_family)]
                if particle_kind is not None
            }
        )
        energy_values = np.asarray(point_data.get("kinetic_energy_eV", []), dtype=float)
        finite_energy_values = energy_values[np.isfinite(energy_values)]
        energy_min_ev = float(finite_energy_values.min()) if finite_energy_values.size else None
        energy_max_ev = float(finite_energy_values.max()) if finite_energy_values.size else None
        return TrajectoryMetadata(
            source_path=source_path,
            available_particle_kinds=tuple(available_particle_kinds),
            trajectory_count=len(lines),
            energy_min_ev=energy_min_ev,
            energy_max_ev=energy_max_ev,
        )

    def get_vtk_structured_scalar_range(self, vtk_path: Path) -> tuple[float, float]:
        with contextlib.redirect_stdout(io.StringIO()):
            dataset = PLOT_UTILS.read_vtk_structured_points(str(vtk_path))

        values = np.asarray(dataset["values"], dtype=float)
        finite_values = values[np.isfinite(values)]
        if finite_values.size == 0:
            raise ValueError(f"No finite scalar values found in {vtk_path.name}")
        return float(finite_values.min()), float(finite_values.max())

    def get_vtk_trajectory_count(self, vtk_path: Path) -> int:
        with contextlib.redirect_stdout(io.StringIO()):
            _points, lines, _point_data, _cell_data = PLOT_UTILS.read_vtk_polydata(str(vtk_path))
        return len(lines)

    def get_vtk_overlay_trajectory_count(self, vtk_path: Path) -> int:
        trajectories_path = self._resolve_field_overlay_path(vtk_path, "trajectories")
        if trajectories_path is None or not trajectories_path.is_file() or trajectories_path.stat().st_size == 0:
            return 0
        return self.get_vtk_trajectory_count(trajectories_path)

    def plot_vtk_structured_slice(
        self,
        vtk_path: Path,
        plane_value: float,
        plane_axis: str = "z",
        show_geometry: bool = False,
        show_trajectories: bool = False,
        trajectories_count: int = -1,
        scalar_min: float | None = None,
        scalar_max: float | None = None,
        scalar_display_mode: str = "Colormap plot",
        contour_lines: int = 10,
        trajectory_options: TrajectoryRenderOptions | None = None,
    ) -> str:
        with contextlib.redirect_stdout(io.StringIO()):
            dataset = PLOT_UTILS.read_vtk_structured_points(str(vtk_path))

        values = np.asarray(dataset["values"], dtype=float)
        finite_values = values[np.isfinite(values)]
        if finite_values.size == 0:
            raise ValueError(f"No finite scalar values found in {vtk_path.name}")

        full_scalar_min = float(finite_values.min())
        full_scalar_max = float(finite_values.max())
        if scalar_min is None:
            scalar_min = full_scalar_min
        if scalar_max is None:
            scalar_max = full_scalar_max
        if scalar_max < scalar_min:
            raise ValueError("Scalar maximum must be greater than or equal to scalar minimum")
        if scalar_max == scalar_min:
            epsilon = max(abs(float(scalar_min)) * 1.0e-12, 1.0e-12)
            scalar_max = float(scalar_min) + epsilon

        origin = tuple(float(v) for v in dataset["origin"])
        spacing = tuple(float(v) for v in dataset["spacing"])
        scalar_name = str(dataset.get("scalar_name") or "value")
        (
            slice_2d,
            horizontal_coords_mm,
            vertical_coords_mm,
            extent,
            x_label,
            y_label,
            plane_selected,
        ) = self._extract_structured_slice(values, origin, spacing, plane_axis, plane_value)

        self.figure.clear()
        axis = self.figure.add_subplot(111)
        scalar_slice = np.ma.masked_invalid(np.asarray(slice_2d, dtype=float))
        display_mode_key = scalar_display_mode.strip().lower()
        if display_mode_key == "contour lines":
            contour_count = max(2, int(contour_lines))
            contour_levels = np.linspace(float(scalar_min), float(scalar_max), num=contour_count)
            contour_set = axis.contour(
                horizontal_coords_mm,
                vertical_coords_mm,
                scalar_slice,
                levels=contour_levels,
                cmap="viridis",
                linewidths=1.0,
            )
            self.figure.colorbar(contour_set, ax=axis, label=scalar_name)
            scalar_plot_description = f"contours, {contour_count} lines"
        else:
            image = axis.imshow(
                scalar_slice,
                origin="lower",
                extent=extent,
                aspect="equal",
                cmap="viridis",
                interpolation="nearest",
                vmin=float(scalar_min),
                vmax=float(scalar_max),
            )
            self.figure.colorbar(image, ax=axis, label=scalar_name)
            scalar_plot_description = "colormap"
        axis.set_title(
            f"{scalar_name} slice at {plane_axis} = {plane_selected * 1000.0:.2f} mm "
            f"(requested {plane_value * 1000.0:.2f} mm)"
        )
        axis.set_xlabel(x_label)
        axis.set_ylabel(y_label)
        _set_equal_spatial_aspect(axis)

        overlay_descriptions: list[str] = []
        if show_geometry and self._overlay_geometry_slice(axis, vtk_path, plane_axis, plane_selected):
            overlay_descriptions.append("geometry")
        if show_trajectories:
            trajectory_count = self._overlay_projected_trajectories(
                axis,
                vtk_path,
                plane_axis,
                trajectories_count,
                trajectory_options,
            )
            if trajectory_count > 0:
                overlay_descriptions.append(f"{trajectory_count} projected trajectories")
        if overlay_descriptions:
            axis.legend(loc="upper right")
        
        axis.grid(alpha=0.15, color="white")
        self.draw_idle()

        if overlay_descriptions:
            return (
                f"Plotted {scalar_name} slice ({scalar_plot_description}) at {plane_axis} = {plane_selected:.6f} m "
                f"with {' and '.join(overlay_descriptions)} and scalar range "
                f"[{scalar_min:.6g}, {scalar_max:.6g}]"
            )
        return (
            f"Plotted {scalar_name} slice ({scalar_plot_description}) at {plane_axis} = {plane_selected:.6f} m "
            f"with scalar range [{scalar_min:.6g}, {scalar_max:.6g}]"
        )

    def plot_vtk_line_plot(
        self,
        vtk_path: Path,
        start_point: tuple[float, float, float],
        direction: tuple[float, float, float],
        num_points: int,
        step_length: float,
    ) -> str:
        """Plot values along a line in the structured field VTK."""
        with contextlib.redirect_stdout(io.StringIO()):
            dataset = PLOT_UTILS.read_vtk_structured_points(str(vtk_path))

        values = np.asarray(dataset["values"], dtype=float)
        origin = tuple(float(v) for v in dataset["origin"])
        spacing = tuple(float(v) for v in dataset["spacing"])
        scalar_name = str(dataset.get("scalar_name") or "value")
        
        # Normalize direction
        dir_array = np.asarray(direction, dtype=float)
        dir_norm = np.linalg.norm(dir_array)
        if dir_norm == 0:
            self.plot_placeholder("Line direction vector has zero magnitude")
            return "Invalid line direction"
        dir_unit = dir_array / dir_norm
        
        # Sample points along the line
        start = np.asarray(start_point, dtype=float)
        line_values = []
        distances = []
        
        for i in range(num_points):
            point = start + dir_unit * step_length * i
            distance = step_length * i

            x_index = int(round((point[0] - origin[0]) / spacing[0]))
            y_index = int(round((point[1] - origin[1]) / spacing[1]))
            z_index = int(round((point[2] - origin[2]) / spacing[2]))
            if (
                x_index < 0
                or y_index < 0
                or z_index < 0
                or x_index >= values.shape[2]
                or y_index >= values.shape[1]
                or z_index >= values.shape[0]
            ):
                break

            line_values.append(float(values[z_index, y_index, x_index]))
            distances.append(distance)
        
        if not line_values:
            self.plot_placeholder("Line path goes outside the field domain")
            return "No values sampled along line"
        
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.plot(distances, line_values, color="#c24d2c", marker="o", linewidth=2, markersize=5)
        axis.set_title(f"{scalar_name} along line from {start_point}")
        axis.set_xlabel("Distance along line [m]")
        axis.set_ylabel(f"{scalar_name}")
        axis.grid(alpha=0.25)
        self.draw_idle()
        
        return f"Plotted {scalar_name} along line ({num_points} points sampled)"

    def _extract_structured_slice(
        self,
        values: np.ndarray,
        origin: tuple[float, float, float],
        spacing: tuple[float, float, float],
        plane_axis: str,
        plane_value: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[float], str, str, float]:
        x_coords_m = origin[0] + spacing[0] * np.arange(values.shape[2], dtype=float)
        y_coords_m = origin[1] + spacing[1] * np.arange(values.shape[1], dtype=float)
        z_coords_m = origin[2] + spacing[2] * np.arange(values.shape[0], dtype=float)
        plane_key = plane_axis.lower()

        if plane_key == "x":
            if len(x_coords_m) == 0:
                raise ValueError("Structured field has no x slices")
            nearest_index = int(np.argmin(np.abs(x_coords_m - plane_value)))
            plane_selected = float(x_coords_m[nearest_index])
            slice_2d = values[:, :, nearest_index].T
            horizontal_coords_mm = z_coords_m * 1000.0
            vertical_coords_mm = y_coords_m * 1000.0
            x_label = "z [mm]"
            y_label = "y [mm]"
        elif plane_key == "y":
            if len(y_coords_m) == 0:
                raise ValueError("Structured field has no y slices")
            nearest_index = int(np.argmin(np.abs(y_coords_m - plane_value)))
            plane_selected = float(y_coords_m[nearest_index])
            slice_2d = values[:, nearest_index, :].T
            horizontal_coords_mm = z_coords_m * 1000.0
            vertical_coords_mm = x_coords_m * 1000.0
            x_label = "z [mm]"
            y_label = "x [mm]"
        else:
            if len(z_coords_m) == 0:
                raise ValueError("Structured field has no z slices")
            nearest_index = int(np.argmin(np.abs(z_coords_m - plane_value)))
            plane_selected = float(z_coords_m[nearest_index])
            slice_2d = values[nearest_index, :, :]
            horizontal_coords_mm = x_coords_m * 1000.0
            vertical_coords_mm = y_coords_m * 1000.0
            x_label = "x [mm]"
            y_label = "y [mm]"

        extent = [
            float(horizontal_coords_mm[0]),
            float(horizontal_coords_mm[-1]),
            float(vertical_coords_mm[0]),
            float(vertical_coords_mm[-1]),
        ]
        return slice_2d, horizontal_coords_mm, vertical_coords_mm, extent, x_label, y_label, plane_selected

    def _resolve_field_overlay_path(self, vtk_path: Path, overlay_kind: str) -> Path | None:
        field_stem = re.sub(r"_(potential|scharge|simulation_state)$", "", vtk_path.stem)
        base_variants: list[str] = []
        is_iteration_field = re.search(r"_it\d+$", field_stem) is not None
        if overlay_kind == "trajectories" and is_iteration_field:
            base_variants.append(field_stem)
        else:
            for candidate in (
                field_stem,
                re.sub(r"_it\d+$", "", field_stem),
                re.sub(r"_simulation$", "", field_stem),
                re.sub(r"_simulation$", "", re.sub(r"_it\d+$", "", field_stem)),
            ):
                if candidate and candidate not in base_variants:
                    base_variants.append(candidate)

        candidate_names: list[str] = []
        for stem in base_variants:
            if overlay_kind == "geometry":
                candidate_names.extend([f"{stem}_geometry_solids.vtk", f"{stem}_simulation_solids.vtk"])
            else:
                candidate_names.extend([f"{stem}_simulation_trajectories.vtk", f"{stem}_trajectories.vtk"])

        empty_match: Path | None = None
        for candidate_name in candidate_names:
            candidate_path = vtk_path.with_name(candidate_name)
            if not candidate_path.is_file():
                continue
            if candidate_path.stat().st_size > 0:
                return candidate_path
            if empty_match is None:
                empty_match = candidate_path
        return empty_match

    def _overlay_geometry_slice(
        self,
        axis: Any,
        vtk_path: Path,
        plane_axis: str,
        plane_value: float,
    ) -> bool:
        geometry_path = self._resolve_field_overlay_path(vtk_path, "geometry")
        if geometry_path is None or not geometry_path.is_file() or geometry_path.stat().st_size == 0:
            return False

        with contextlib.redirect_stdout(io.StringIO()):
            dataset = PLOT_UTILS.read_vtk_structured_points(str(geometry_path))

        geometry_values = np.asarray(dataset["values"], dtype=float)
        geometry_origin = tuple(float(v) for v in dataset["origin"])
        geometry_spacing = tuple(float(v) for v in dataset["spacing"])
        (
            geometry_slice,
            horizontal_coords_mm,
            vertical_coords_mm,
            _extent,
            _x_label,
            _y_label,
            _plane_selected,
        ) = self._extract_structured_slice(geometry_values, geometry_origin, geometry_spacing, plane_axis, plane_value)

        solid_mask = np.not_equal(geometry_slice, 0.0)
        if not np.any(solid_mask) or np.all(solid_mask):
            return False

        axis.contourf(
            horizontal_coords_mm,
            vertical_coords_mm,
            solid_mask.astype(float),
            levels=[0.5, 1.5],
            colors=[GEOMETRY_FILL_COLOR],
            alpha=0.55,
            antialiased=True,
        )
        axis.contour(
            horizontal_coords_mm,
            vertical_coords_mm,
            solid_mask.astype(float),
            levels=[0.5],
            colors=GEOMETRY_EDGE_COLOR,
            linewidths=1.0,
            alpha=0.95,
        )
        axis.plot([], [], color=GEOMETRY_EDGE_COLOR, linewidth=1.0, label="Geometry solids")
        return True

    def _overlay_projected_trajectories(
        self,
        axis: Any,
        vtk_path: Path,
        plane_axis: str,
        trajectories_count: int,
        trajectory_options: TrajectoryRenderOptions | None = None,
    ) -> int:
        trajectories_path = self._resolve_field_overlay_path(vtk_path, "trajectories")
        if trajectories_path is None or not trajectories_path.is_file() or trajectories_path.stat().st_size == 0:
            return 0

        with contextlib.redirect_stdout(io.StringIO()):
            points, lines, point_data, cell_data = PLOT_UTILS.read_vtk_polydata(str(trajectories_path))

        if not lines or points.size == 0:
            return 0

        options = _normalize_trajectory_options(trajectory_options)
        selected_particle_kinds = None
        if options.selected_particle_kinds is not None:
            selected_particle_kinds = {int(value) for value in options.selected_particle_kinds}

        available_indices: list[int] = []
        particle_kinds_by_index: dict[int, int | None] = {}
        for line_index in range(len(lines)):
            particle_kind = _resolve_particle_kind_from_cell_data(cell_data, line_index, options.particle_family)
            if selected_particle_kinds is not None and particle_kind not in selected_particle_kinds:
                continue
            available_indices.append(line_index)
            particle_kinds_by_index[line_index] = particle_kind

        if not available_indices:
            return 0

        selected_indices = np.asarray(available_indices, dtype=int)
        if trajectories_count > 0 and trajectories_count < len(selected_indices):
            rng = np.random.default_rng(0)
            selected_indices = np.sort(rng.choice(selected_indices, size=trajectories_count, replace=False))

        color_mode = _trajectory_color_mode_key(options.color_mode)
        visible_trajectory_count = 0

        if color_mode == TRAJECTORY_COLOR_MODE_ENERGY.lower():
            projected_segments: list[np.ndarray] = []
            segment_energies: list[float] = []
            for line_index in selected_indices:
                trajectory = np.asarray(points[lines[int(line_index)]], dtype=float)
                energies = _trajectory_energy_values(point_data, lines[int(line_index)])
                if energies is None:
                    continue
                visible = _trajectory_visibility_mask(len(trajectory), energies, options)
                segments_3d, energy_values = _visible_segment_pairs(trajectory, visible, energies)
                if not segments_3d:
                    continue
                local_projected: list[np.ndarray] = []
                for segment in segments_3d:
                    projected = self._project_trajectory_to_plane(segment, plane_axis)
                    if projected is not None:
                        local_projected.append(projected)
                if not local_projected:
                    continue
                projected_segments.extend([segment * 1000.0 for segment in local_projected])
                segment_energies.extend(energy_values[: len(local_projected)])
                visible_trajectory_count += 1

            if not projected_segments:
                return 0

            collection = LineCollection(
                projected_segments,
                cmap="plasma",
                norm=_build_normalize(min(segment_energies), max(segment_energies)),
                linewidths=0.7,
                alpha=0.6,
            )
            collection.set_array(np.asarray(segment_energies, dtype=float))
            axis.add_collection(collection)
            self.figure.colorbar(collection, ax=axis, label="Kinetic energy [eV]")
            axis.plot([], [], color="#ff9f1c", linewidth=1.2, alpha=0.7, label="Trajectories")
            return visible_trajectory_count

        projected_lines: list[np.ndarray] = []
        line_colors: list[Any] = []
        line_statuses: list[float] = []
        species_used: set[int] = set()
        status_values = np.asarray(cell_data.get("particle_status", []), dtype=float)
        status_max = max(float(np.max(status_values)), 1.0) if status_values.size else 1.0

        for line_index in selected_indices:
            trajectory = np.asarray(points[lines[int(line_index)]], dtype=float)
            energies = _trajectory_energy_values(point_data, lines[int(line_index)])
            visible = _trajectory_visibility_mask(len(trajectory), energies, options)
            visible_segments = _split_visible_polylines(trajectory, visible)
            if not visible_segments:
                continue

            particle_kind = particle_kinds_by_index.get(int(line_index))
            local_projected: list[np.ndarray] = []
            for segment in visible_segments:
                projected = self._project_trajectory_to_plane(segment, plane_axis)
                if projected is not None:
                    local_projected.append(projected)
            if not local_projected:
                continue

            visible_trajectory_count += 1
            projected_lines.extend(local_projected)
            if color_mode == TRAJECTORY_COLOR_MODE_SPECIES.lower() and particle_kind is not None:
                species_used.add(particle_kind)
                line_colors.extend([particle_color_from_export_id(particle_kind)] * len(local_projected))
            elif (
                color_mode == TRAJECTORY_COLOR_MODE_STATUS.lower()
                and int(line_index) < len(status_values)
                and math.isfinite(float(status_values[int(line_index)]))
            ):
                line_statuses.extend([float(status_values[int(line_index)])] * len(local_projected))
            else:
                line_colors.extend(["#ff9f1c"] * len(local_projected))

        if not projected_lines:
            return 0

        if color_mode == TRAJECTORY_COLOR_MODE_STATUS.lower() and line_statuses:
            collection = LineCollection(
                [line * 1000.0 for line in projected_lines],
                cmap="viridis",
                norm=_build_normalize(0.0, status_max),
                linewidths=0.7,
                alpha=0.55,
            )
            collection.set_array(np.asarray(line_statuses, dtype=float))
            axis.add_collection(collection)
            status_mappable = cm.ScalarMappable(norm=_build_normalize(0.0, status_max), cmap="viridis")
            status_mappable.set_array(np.asarray([0.0, status_max], dtype=float))
            self.figure.colorbar(status_mappable, ax=axis, label="Particle status")
            axis.plot([], [], color="#ff9f1c", linewidth=1.2, alpha=0.7, label="Trajectories")
            return visible_trajectory_count

        collection = LineCollection(
            [line * 1000.0 for line in projected_lines],
            colors=line_colors if line_colors else "#ff9f1c",
            linewidths=0.7,
            alpha=0.55,
        )
        axis.add_collection(collection)
        if color_mode == TRAJECTORY_COLOR_MODE_SPECIES.lower() and species_used:
            _add_species_legend(axis, species_used, options.particle_family, marker=False)
        else:
            axis.plot([], [], color="#ff9f1c", linewidth=1.2, alpha=0.7, label="Trajectories")
        return visible_trajectory_count

    def _project_trajectory_to_plane(self, trajectory: np.ndarray, plane_axis: str) -> np.ndarray | None:
        if trajectory.ndim != 2 or trajectory.shape[0] < 2 or trajectory.shape[1] != 3:
            return None

        plane_key = plane_axis.lower()
        if plane_key == "x":
            projected = trajectory[:, [2, 1]]
        elif plane_key == "y":
            projected = trajectory[:, [2, 0]]
        else:
            projected = trajectory[:, [0, 1]]

        return np.asarray(projected, dtype=float)

    def plot_runtime_geometry_3d(self, runtime_case: dict[str, Any]) -> str:
        geometry = runtime_case.get("geometry") if isinstance(runtime_case.get("geometry"), dict) else {}
        domain = geometry.get("domain") if isinstance(geometry.get("domain"), dict) else None
        solids = geometry.get("solids") if isinstance(geometry.get("solids"), list) else []
        sources = runtime_case.get("particleSources")
        return self._geometry_canvas().plot_geometry_sources_3d(solids, sources, domain=domain, show_solid=True)

    def plot_runtime_geometry_section(self, runtime_case: dict[str, Any]) -> str:
        geometry = runtime_case.get("geometry") if isinstance(runtime_case.get("geometry"), dict) else {}
        domain = geometry.get("domain") if isinstance(geometry.get("domain"), dict) else None
        solids = geometry.get("solids") if isinstance(geometry.get("solids"), list) else []
        sources = runtime_case.get("particleSources")
        return self._geometry_canvas().plot_geometry_sources_section(solids, sources, domain=domain)

    def plot_diagnostic_summary(self, summary_path: Path, metric: str) -> str:
        parsed = parse_diagnostic_summary_txt(summary_path)
        rows = parsed.get("rows", [])
        if not rows:
            self.plot_placeholder(f"No diagnostic rows found in\n{summary_path.name}")
            return "No diagnostic rows found"

        z_values = np.asarray([float(row["z[mm]"]) for row in rows], dtype=float)
        
        # Handle grouped metrics
        if metric == "centroids_grouped":
            centroid_x, _, x_ylabel = self._diagnostic_metric_series(rows, "centroid_x")
            centroid_y, _, y_ylabel = self._diagnostic_metric_series(rows, "centroid_y")
            self._plot_grouped_metrics(
                (1, 2),
                z_values,
                [("x", centroid_x, "Beam Centroid x", x_ylabel),
                 ("y", centroid_y, "Beam Centroid y", y_ylabel)],
                "vs z",
                "z [mm]",
            )
            return f"Plotted beam centroids from {summary_path.name}"
        
        if metric == "divergences_grouped":
            divergence_x, _, x_ylabel = self._diagnostic_metric_series(rows, "divergence_x")
            divergence_y, _, y_ylabel = self._diagnostic_metric_series(rows, "divergence_y")
            self._plot_grouped_metrics(
                (1, 2),
                z_values,
                [("Dx", divergence_x, "Beam Divergence Dx", x_ylabel),
                 ("Dy", divergence_y, "Beam Divergence Dy", y_ylabel)],
                "vs z",
                "z [mm]",
            )
            return f"Plotted beam divergences from {summary_path.name}"
        
        if metric == "deviations_grouped":
            deviation_x, _, x_ylabel = self._diagnostic_metric_series(rows, "deviation_x")
            deviation_y, _, y_ylabel = self._diagnostic_metric_series(rows, "deviation_y")
            self._plot_grouped_metrics(
                (1, 2),
                z_values,
                [("x'", deviation_x, "Average Deviation x'", x_ylabel),
                 ("y'", deviation_y, "Average Deviation y'", y_ylabel)],
                "vs z",
                "z [mm]",
            )
            return f"Plotted average deviations from {summary_path.name}"
        
        y_values, title, ylabel = self._diagnostic_metric_series(rows, metric)

        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.plot(z_values, y_values, color="#c24d2c", marker="o", linewidth=1.8, markersize=4)
        axis.set_title(title)
        axis.set_xlabel("z [mm]")
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.25)
        if len(z_values) > 1:
            axis.set_xlim(float(np.min(z_values)), float(np.max(z_values)))
        self.draw_idle()
        return f"Plotted {metric} from {summary_path.name}"

    def plot_iteration_summary(self, summary_path: Path, metric: str) -> str:
        parsed = parse_diagnostic_summary_txt(summary_path)
        rows = parsed.get("rows", [])
        if not rows:
            self.plot_placeholder(f"No iteration rows found in\n{summary_path.name}")
            return "No iteration rows found"

        iteration_values = np.asarray([float(row["it"]) for row in rows], dtype=float)
        
        # Handle grouped metrics
        if metric == "centroids_grouped":
            centroid_x, _, x_ylabel = self._diagnostic_metric_series(rows, "centroid_x")
            centroid_y, _, y_ylabel = self._diagnostic_metric_series(rows, "centroid_y")
            self._plot_grouped_metrics(
                (1, 2),
                iteration_values,
                [("x", centroid_x, "Beam Centroid x", x_ylabel),
                 ("y", centroid_y, "Beam Centroid y", y_ylabel)],
                "vs iteration",
                "iteration",
            )
            return f"Plotted beam centroids from {summary_path.name}"
        
        if metric == "divergences_grouped":
            divergence_x, _, x_ylabel = self._diagnostic_metric_series(rows, "divergence_x")
            divergence_y, _, y_ylabel = self._diagnostic_metric_series(rows, "divergence_y")
            self._plot_grouped_metrics(
                (1, 2),
                iteration_values,
                [("Dx", divergence_x, "Beam Divergence Dx", x_ylabel),
                 ("Dy", divergence_y, "Beam Divergence Dy", y_ylabel)],
                "vs iteration",
                "iteration",
            )
            return f"Plotted beam divergences from {summary_path.name}"
        
        if metric == "deviations_grouped":
            deviation_x, _, x_ylabel = self._diagnostic_metric_series(rows, "deviation_x")
            deviation_y, _, y_ylabel = self._diagnostic_metric_series(rows, "deviation_y")
            self._plot_grouped_metrics(
                (1, 2),
                iteration_values,
                [("x'", deviation_x, "Average Deviation x'", x_ylabel),
                 ("y'", deviation_y, "Average Deviation y'", y_ylabel)],
                "vs iteration",
                "iteration",
            )
            return f"Plotted average deviations from {summary_path.name}"
        
        if metric == "convergence":
            self._plot_convergence_metrics(iteration_values, rows, "iteration")
            return f"Plotted convergence metrics from {summary_path.name}"
        
        y_values, title, ylabel = self._diagnostic_metric_series(rows, metric)

        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.plot(iteration_values, y_values, color="#2f6f7e", marker="o", linewidth=1.8, markersize=4)
        axis.set_title(title.replace(" vs z", " vs iteration"))
        axis.set_xlabel("iteration")
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.25)
        if len(iteration_values) > 1:
            axis.set_xlim(float(np.min(iteration_values)), float(np.max(iteration_values)))
        self.draw_idle()
        return f"Plotted {metric} from {summary_path.name}"

    def plot_grid_power_summary(self, summary_path: Path, metric: str) -> str:
        parsed = parse_grid_power_summary_txt(summary_path)
        rows = parsed.get("rows", [])
        if not rows:
            self.plot_placeholder(f"No grid-power rows found in\n{summary_path.name}")
            return "No grid-power rows found"

        descriptions = [str(row["Description"]) for row in rows]
        if metric == "current":
            values = np.asarray([float(row["Current[A]"]) for row in rows], dtype=float)
            ylabel = "Current [A]"
            title = "Grid Current Summary (gross, arriving)"
            color = "#2f6f7e"
        elif metric == "particles":
            values = np.asarray([float(row["Particles"]) for row in rows], dtype=float)
            ylabel = "Particles"
            title = "Grid Particle Counts"
            color = "#8d6a9f"
        elif metric == "net_power":
            # Gross minus what the surface emitted back. parse_grid_power_summary_txt
            # merges these in from the file's "# NetRow" lines and falls back to gross
            # where no emission ledger exists, so the key is always present.
            values = np.asarray([float(row["NetPower[W]"]) for row in rows], dtype=float)
            ylabel = "Net power [W]"
            title = "Grid Power Summary (net of emitted secondaries)"
            color = "#a8452a"
        elif metric == "net_current":
            values = np.asarray([float(row["NetCurrent[A]"]) for row in rows], dtype=float)
            ylabel = "Net current [A]"
            title = "Grid Current Summary (net drain)"
            color = "#25525d"
        else:
            values = np.asarray([float(row["Power[W]"]) for row in rows], dtype=float)
            ylabel = "Power [W]"
            title = "Grid Power Summary (gross, arriving)"
            color = "#c24d2c"

        positions = np.arange(len(descriptions), dtype=float)
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.bar(positions, values, color=color, alpha=0.85)
        axis.set_title(title)
        axis.set_ylabel(ylabel)
        axis.set_xticks(positions)
        axis.set_xticklabels(descriptions, rotation=35, ha="right")
        axis.grid(axis="y", alpha=0.25)
        self.draw_idle()
        return f"Plotted {metric} from {summary_path.name}"

    def plot_emitter_map(self, emitter_path: Path) -> str:
        coordinates: list[tuple[float, float]] = []
        with emitter_path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                tokens = line.split()
                if len(tokens) < 2:
                    continue
                try:
                    coordinates.append((float(tokens[0]) * 1000.0, float(tokens[1]) * 1000.0))
                except ValueError:
                    continue

        if not coordinates:
            self.plot_placeholder(f"No emitter coordinates found in\n{emitter_path.name}")
            return "No emitter coordinates found"

        samples = np.asarray(coordinates, dtype=float)
        x_values = samples[:, 0]
        y_values = samples[:, 1]
        bins = max(24, min(72, int(math.sqrt(len(samples)))))

        self.figure.clear()
        scatter_axis = self.figure.add_subplot(121)
        density_axis = self.figure.add_subplot(122)

        scatter_axis.scatter(x_values, y_values, s=6, alpha=0.4, color="#2f6f7e")
        scatter_axis.set_title("Emitter Exit Footprint")
        scatter_axis.set_xlabel("x [mm]")
        scatter_axis.set_ylabel("y [mm]")
        scatter_axis.grid(alpha=0.2)
        scatter_axis.set_aspect("equal", adjustable="box")

        hist, x_edges, y_edges = np.histogram2d(x_values, y_values, bins=bins)
        image = density_axis.imshow(
            hist.T,
            origin="lower",
            extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
            aspect="auto",
            cmap="magma",
            interpolation="bilinear",
        )
        self.figure.colorbar(image, ax=density_axis, label="Particle count")
        density_axis.set_title("Emitter Density")
        density_axis.set_xlabel("x [mm]")
        density_axis.set_ylabel("y [mm]")
        density_axis.grid(alpha=0.15, color="white")

        self.draw_idle()
        return f"Plotted emitter footprint from {emitter_path.name}"

    def _diagnostic_metric_series(
        self,
        rows: list[dict[str, int | float]],
        metric: str,
    ) -> tuple[np.ndarray, str, str]:
        if metric == "beam_current":
            return (
                np.asarray([float(row["I[mA]"]) for row in rows], dtype=float),
                "Beam Current vs z",
                "Current [mA]",
            )
        if metric == "centroid_x":
            return (
                np.asarray([float(row["<x>[mm]"]) for row in rows], dtype=float),
                "Beam Centroid x vs z",
                "<x> [mm]",
            )
        if metric == "centroid_y":
            return (
                np.asarray([float(row["<y>[mm]"]) for row in rows], dtype=float),
                "Beam Centroid y vs z",
                "<y> [mm]",
            )
        if metric == "divergence_x":
            return (
                np.asarray([float(row["Dx[mrad]"]) for row in rows], dtype=float),
                "Beam Divergence Dx vs z",
                "Dx [mrad]",
            )
        if metric == "divergence_y":
            return (
                np.asarray([float(row["Dy[mrad]"]) for row in rows], dtype=float),
                "Beam Divergence Dy vs z",
                "Dy [mrad]",
            )
        if metric == "deviation_x":
            return (
                np.asarray([float(row.get("<x'>[mrad]", 0.0)) for row in rows], dtype=float),
                "Average Deviation x' vs z",
                "<x'> [mrad]",
            )
        if metric == "deviation_y":
            return (
                np.asarray([float(row.get("<y'>[mrad]", 0.0)) for row in rows], dtype=float),
                "Average Deviation y' vs z",
                "<y'> [mrad]",
            )
        return (
            np.asarray([float(row["<V>[V]"]) for row in rows], dtype=float),
            "Average Voltage vs z",
            "<V> [V]",
        )

    def _plot_grouped_metrics(
        self,
        axis_indices: tuple[int, int],
        x_values: np.ndarray,
        metric_pairs: list[tuple[str, np.ndarray, str, str]],
        title_suffix: str,
        x_label: str,
    ) -> None:
        """Plot two metrics side by side."""
        self.figure.clear()
        
        for plot_idx, (axis_idx, (metric_label, y_values, metric_title, y_label)) in enumerate(
            zip(axis_indices, metric_pairs)
        ):
            axis = self.figure.add_subplot(1, 2, axis_idx)
            axis.plot(x_values, y_values, color="#c24d2c", marker="o", linewidth=1.8, markersize=4)
            axis.set_title(f"{metric_title} {title_suffix}")
            axis.set_xlabel(x_label)
            axis.set_ylabel(y_label)
            axis.grid(alpha=0.25)
            if len(x_values) > 1:
                axis.set_xlim(float(np.min(x_values)), float(np.max(x_values)))
        
        self.draw_idle()

    def _plot_convergence_metrics(
        self,
        iteration_values: np.ndarray,
        rows: list[dict[str, int | float]],
        x_label: str,
    ) -> None:
        """Plot convergence metrics (ratio i/i-1 with log scale)."""
        self.figure.clear()
        
        # Plot 1: Beam current
        ax1 = self.figure.add_subplot(2, 2, 1)
        current_values = np.asarray([float(row.get("I[mA]", 0.0)) for row in rows], dtype=float)
        ax1.plot(iteration_values, current_values, color="#c24d2c", marker="o", linewidth=1.8, markersize=4)
        ax1.set_title("Beam Current vs Iteration")
        ax1.set_xlabel(x_label)
        ax1.set_ylabel("Current [mA]")
        ax1.grid(alpha=0.25)
        
        # Plot 2: Deviations
        ax2 = self.figure.add_subplot(2, 2, 2)
        deviation_x = np.asarray([float(row.get("<x'>[mrad]", 0.0)) for row in rows], dtype=float)
        deviation_y = np.asarray([float(row.get("<y'>[mrad]", 0.0)) for row in rows], dtype=float)
        ax2.plot(iteration_values, deviation_x, color="#2f6f7e", marker="o", linewidth=1.8, markersize=4, label="x'")
        ax2.plot(iteration_values, deviation_y, color="#8d6a9f", marker="s", linewidth=1.8, markersize=4, label="y'")
        ax2.set_title("Beam Deviations vs Iteration")
        ax2.set_xlabel(x_label)
        ax2.set_ylabel("Deviation [mrad]")
        ax2.grid(alpha=0.25)
        ax2.legend()
        
        # Plot 3: Divergences
        ax3 = self.figure.add_subplot(2, 2, 3)
        divergence_x = np.asarray([float(row.get("Dx[mrad]", 0.0)) for row in rows], dtype=float)
        divergence_y = np.asarray([float(row.get("Dy[mrad]", 0.0)) for row in rows], dtype=float)
        ax3.plot(iteration_values, divergence_x, color="#2f6f7e", marker="o", linewidth=1.8, markersize=4, label="Dx")
        ax3.plot(iteration_values, divergence_y, color="#8d6a9f", marker="s", linewidth=1.8, markersize=4, label="Dy")
        ax3.set_title("Beam Divergences vs Iteration")
        ax3.set_xlabel(x_label)
        ax3.set_ylabel("Divergence [mrad]")
        ax3.grid(alpha=0.25)
        ax3.legend()
        
        # Plot 4: Convergence ratios (i/i-1 - 1) with log scale
        ax4 = self.figure.add_subplot(2, 2, 4)
        
        if len(iteration_values) > 1:
            convergence_current = []
            convergence_div_x = []
            convergence_div_y = []
            iteration_conv = []
            
            for i in range(1, len(rows)):
                prev_current = float(rows[i - 1].get("I[mA]", 1.0)) or 1.0
                curr_current = float(rows[i].get("I[mA]", 1.0)) or 1.0
                convergence_current.append(abs(curr_current / prev_current - 1.0))
                
                prev_div_x = float(rows[i - 1].get("Dx[mrad]", 1.0)) or 1.0
                curr_div_x = float(rows[i].get("Dx[mrad]", 1.0)) or 1.0
                convergence_div_x.append(abs(curr_div_x / prev_div_x - 1.0))
                
                prev_div_y = float(rows[i - 1].get("Dy[mrad]", 1.0)) or 1.0
                curr_div_y = float(rows[i].get("Dy[mrad]", 1.0)) or 1.0
                convergence_div_y.append(abs(curr_div_y / prev_div_y - 1.0))
                
                iteration_conv.append(float(iteration_values[i]))
            
            if convergence_current:
                convergence_current = np.asarray(convergence_current, dtype=float)
                convergence_div_x = np.asarray(convergence_div_x, dtype=float)
                convergence_div_y = np.asarray(convergence_div_y, dtype=float)
                iteration_conv = np.asarray(iteration_conv, dtype=float)
                
                ax4.semilogy(iteration_conv, convergence_current, color="#c24d2c", marker="o", linewidth=1.8, markersize=4, label="Current")
                ax4.semilogy(iteration_conv, convergence_div_x, color="#2f6f7e", marker="s", linewidth=1.8, markersize=4, label="Dx")
                ax4.semilogy(iteration_conv, convergence_div_y, color="#8d6a9f", marker="^", linewidth=1.8, markersize=4, label="Dy")
                ax4.set_title("Convergence: |x(i)/x(i-1) - 1|")
                ax4.set_xlabel(x_label)
                ax4.set_ylabel("Relative Change (log scale)")
                ax4.grid(alpha=0.25, which="both")
                ax4.legend()
        
        self.draw_idle()

    def _trajectory_canvas(self) -> TrajectoryCanvas:
        proxy = TrajectoryCanvas.__new__(TrajectoryCanvas)
        proxy.figure = self.figure
        proxy.draw_idle = self.draw_idle
        return proxy

    def _geometry_canvas(self) -> GeometryCanvas:
        proxy = GeometryCanvas.__new__(GeometryCanvas)
        proxy.figure = self.figure
        proxy._redraw = self.draw_idle
        return proxy
