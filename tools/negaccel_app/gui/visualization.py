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

        ranges = np.array([
            points[:, 0].max() - points[:, 0].min(),
            points[:, 1].max() - points[:, 1].min(),
            points[:, 2].max() - points[:, 2].min(),
        ])
        max_range = float(ranges.max()) / 2.0 if len(points) else 1.0
        mid_x = float(points[:, 0].max() + points[:, 0].min()) / 2.0
        mid_y = float(points[:, 1].max() + points[:, 1].min()) / 2.0
        mid_z = float(points[:, 2].max() + points[:, 2].min()) / 2.0

        axis.set_xlim((mid_x - max_range) * 1000.0, (mid_x + max_range) * 1000.0)
        axis.set_ylim((mid_y - max_range) * 1000.0, (mid_y + max_range) * 1000.0)
        axis.set_zlim((mid_z - max_range) * 1000.0, (mid_z + max_range) * 1000.0)
        axis.set_title(f"3D Trajectories ({count} shown)")
        axis.set_xlabel("x [mm]")
        axis.set_ylabel("y [mm]")
        axis.set_zlabel("z [mm]")
        self.draw_idle()
        return f"Plotted {count} trajectories in 3D"
