"""Visualization helpers for NegAccel GUI tooling."""

from __future__ import annotations

import contextlib
import importlib.util
import io
from pathlib import Path
from typing import Any

try:
    import matplotlib.cm as cm
    import numpy as np
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
    from matplotlib.figure import Figure
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
except ImportError as exc:  # pragma: no cover - depends on the local environment
    raise SystemExit(
        "numpy and matplotlib are required for the GUI. Install the packages listed in tools/gui_requirements.txt."
    ) from exc

from .common import REPO_ROOT


PLOT_UTILS_PATH = REPO_ROOT / "plot_trajectories_vtk.py"


def load_plot_utils() -> Any:
    spec = importlib.util.spec_from_file_location("negaccel_plot_utils", PLOT_UTILS_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load plotting helpers from {PLOT_UTILS_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PLOT_UTILS = load_plot_utils()


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


class GeometryCanvas(FigureCanvasQTAgg):
    def __init__(self) -> None:
        self.figure = Figure(figsize=(8, 6), tight_layout=True)
        super().__init__(self.figure)
        self.plot_placeholder("Select a geometry preview mode to inspect the current solid setup.")

    def _redraw(self) -> None:
        self.draw()

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

    def _apply_3d_view(self, axis: Any, view_direction: str | None) -> None:
        if not view_direction:
            return

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
        self.plot_placeholder("Materialize or run a case, then choose a VTK file to visualize.")

    def plot_placeholder(self, message: str) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.axis("off")
        axis.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
        self.draw_idle()

    def _read_vtk(self, vtk_path: Path) -> tuple[np.ndarray, list[np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray]]:
        with contextlib.redirect_stdout(io.StringIO()):
            return PLOT_UTILS.read_vtk_polydata(str(vtk_path))

    def plot_cross_section(self, vtk_path: Path, z_plane: float) -> str:
        points, lines, _, cell_data = self._read_vtk(vtk_path)
        crossings: list[tuple[float, float]] = []
        particle_status: list[float] = []

        for index, line_indices in enumerate(lines):
            crossing = PLOT_UTILS.interpolate_trajectory_at_z(points, line_indices, z_plane)
            if crossing is None:
                continue
            crossings.append(crossing)
            if "particle_status" in cell_data and index < len(cell_data["particle_status"]):
                particle_status.append(float(cell_data["particle_status"][index]))

        if not crossings:
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

        if particle_status and len(particle_status) == len(crossings):
            scatter = scatter_axis.scatter(
                x_values * 1000.0,
                y_values * 1000.0,
                c=np.array(particle_status),
                cmap="viridis",
                s=20,
                alpha=0.75,
            )
            self.figure.colorbar(scatter, ax=scatter_axis, label="Particle status")
        else:
            scatter_axis.scatter(x_values * 1000.0, y_values * 1000.0, s=20, alpha=0.75, color="#c24d2c")

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
            aspect="auto",
            cmap="magma",
            interpolation="bilinear",
        )
        self.figure.colorbar(image, ax=density_axis, label="Particle count")
        density_axis.set_title("Cross-section density")
        density_axis.set_xlabel("x [mm]")
        density_axis.set_ylabel("y [mm]")
        density_axis.grid(alpha=0.2, color="white")

        stats_text = (
            f"N = {len(crossings)}\n"
            f"<x> = {x_mean * 1000.0:.3f} mm\n"
            f"<y> = {y_mean * 1000.0:.3f} mm\n"
            f"sigma_x = {x_std * 1000.0:.3f} mm\n"
            f"sigma_y = {y_std * 1000.0:.3f} mm"
        )
        density_axis.text(
            0.03,
            0.97,
            stats_text,
            transform=density_axis.transAxes,
            va="top",
            fontsize=9,
            bbox={"boxstyle": "round", "facecolor": "#fff5e1", "alpha": 0.9},
        )

        self.draw_idle()
        return f"Plotted {len(crossings)} crossings at z = {z_plane:.4f} m"

    def plot_3d(self, vtk_path: Path, max_trajectories: int) -> str:
        points, lines, _, cell_data = self._read_vtk(vtk_path)
        count = min(max_trajectories, len(lines))
        if count == 0:
            self.plot_placeholder("The selected VTK file does not contain trajectory lines.")
            return "No trajectories available"

        self.figure.clear()
        axis = self.figure.add_subplot(111, projection="3d")
        axis.set_proj_type("ortho")

        status_max = 1.0
        if "particle_status" in cell_data and len(cell_data["particle_status"]) > 0:
            status_max = max(float(np.max(cell_data["particle_status"])), 1.0)

        for index in range(count):
            trajectory = points[lines[index]]
            if "particle_status" in cell_data and index < len(cell_data["particle_status"]):
                mapped = float(cell_data["particle_status"][index]) / status_max
                color = cm.viridis(mapped)
            else:
                color = "#2f6f7e"
            axis.plot(
                trajectory[:, 0] * 1000.0,
                trajectory[:, 1] * 1000.0,
                trajectory[:, 2] * 1000.0,
                color=color,
                linewidth=0.8,
                alpha=0.65,
            )

        _set_adaptive_3d_bounds(axis, points, scale=1000.0)
        axis.set_title(f"3D Trajectories ({count} shown)")
        axis.set_xlabel("x [mm]")
        axis.set_ylabel("y [mm]")
        axis.set_zlabel("z [mm]")
        self.draw_idle()
        return f"Plotted {count} trajectories in 3D"
