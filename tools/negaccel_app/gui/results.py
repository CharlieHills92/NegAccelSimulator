"""VTK result browsing and plotting helpers for the NegAccel GUI."""

from __future__ import annotations

from pathlib import Path

from .common import (
    QComboBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)
from .visualization import NavigationToolbar2QT, TrajectoryCanvas


class ResultsMixin:
    def _build_detail_tabs(self) -> QTabWidget:
        self.output_tabs = QTabWidget()

        self.authoring_preview = QPlainTextEdit()
        self.authoring_preview.setReadOnly(True)
        self.output_tabs.addTab(self.authoring_preview, "Authoring JSON")

        self.runtime_preview = QPlainTextEdit()
        self.runtime_preview.setReadOnly(True)
        self.output_tabs.addTab(self.runtime_preview, "Runtime JSON")

        self.process_output = QPlainTextEdit()
        self.process_output.setReadOnly(True)
        self.output_tabs.addTab(self.process_output, "Run Log")

        visualization_tab = QWidget()
        visualization_layout = QVBoxLayout(visualization_tab)
        controls_row = QHBoxLayout()
        self.vtk_file_combo = QComboBox()
        self.vtk_file_combo.setMinimumWidth(320)
        self.vtk_mode_combo = QComboBox()
        self.vtk_mode_combo.addItems(["Cross-section", "3D trajectories"])
        self.z_plane_spin = QDoubleSpinBox()
        self.z_plane_spin.setRange(-10.0, 10.0)
        self.z_plane_spin.setDecimals(4)
        self.z_plane_spin.setSingleStep(0.001)
        self.z_plane_spin.setValue(0.565)
        self.max_trajectories_spin = QSpinBox()
        self.max_trajectories_spin.setRange(1, 5000)
        self.max_trajectories_spin.setValue(200)
        refresh_vtk_button = QPushButton("Refresh VTK Files")
        refresh_vtk_button.clicked.connect(self.refresh_vtk_files)
        plot_button = QPushButton("Plot")
        plot_button.clicked.connect(self.plot_selected_vtk)

        controls_row.addWidget(QLabel("VTK file"))
        controls_row.addWidget(self.vtk_file_combo, 1)
        controls_row.addWidget(QLabel("Mode"))
        controls_row.addWidget(self.vtk_mode_combo)
        controls_row.addWidget(QLabel("z-plane [m]"))
        controls_row.addWidget(self.z_plane_spin)
        controls_row.addWidget(QLabel("Max trajectories"))
        controls_row.addWidget(self.max_trajectories_spin)
        controls_row.addWidget(refresh_vtk_button)
        controls_row.addWidget(plot_button)
        visualization_layout.addLayout(controls_row)

        self.visualization_canvas = TrajectoryCanvas()
        self.visualization_toolbar = NavigationToolbar2QT(self.visualization_canvas, visualization_tab)
        self.visualization_status = QLabel("No VTK file selected yet.")
        visualization_layout.addWidget(self.visualization_toolbar)
        visualization_layout.addWidget(self.visualization_canvas, 1)
        visualization_layout.addWidget(self.visualization_status)
        self.output_tabs.addTab(visualization_tab, "Visualization")

        return self.output_tabs

    def refresh_vtk_files(self) -> None:
        current_selection = self.vtk_file_combo.currentData()
        self.vtk_file_combo.clear()
        vtk_dir = self.resolve_vtk_directory()
        if vtk_dir is None or not vtk_dir.exists():
            self.vtk_file_combo.addItem("No VTK directory available", "")
            self.visualization_status.setText("No VTK directory available yet.")
            return

        vtk_files = sorted(vtk_dir.glob("*.vtk"))
        if not vtk_files:
            self.vtk_file_combo.addItem(f"No VTK files in {vtk_dir}", "")
            self.visualization_status.setText(f"No VTK files found in {vtk_dir}.")
            return

        selected_index = 0
        for index, path in enumerate(vtk_files):
            self.vtk_file_combo.addItem(path.name, str(path))
            if current_selection and str(path) == str(current_selection):
                selected_index = index
        self.vtk_file_combo.setCurrentIndex(selected_index)
        self.visualization_status.setText(f"Found {len(vtk_files)} VTK files in {vtk_dir}.")

    def resolve_vtk_directory(self) -> Path | None:
        if self.current_runtime_case is None:
            try:
                runtime_case, runtime_path = self.build_runtime_case_preview()
            except Exception:
                return None
        else:
            runtime_case = self.current_runtime_case
            runtime_path = self.current_runtime_path

        if runtime_path is None:
            return None

        outputs = runtime_case.get("outputs", {})
        vtk_config = outputs.get("vtk", {})
        vtk_directory = vtk_config.get("directory", "VTK")
        root_directory = outputs.get("rootDirectory") or str(runtime_path.parent)

        vtk_path = Path(vtk_directory)
        if not vtk_path.is_absolute():
            vtk_path = Path(root_directory) / vtk_path
        return vtk_path

    def plot_selected_vtk(self) -> None:
        selected_path = self.vtk_file_combo.currentData()
        if not selected_path:
            self.visualization_canvas.plot_placeholder("No VTK file is currently available.")
            self.visualization_status.setText("No VTK file is currently available.")
            return

        vtk_path = Path(str(selected_path))
        if not vtk_path.exists():
            self.visualization_canvas.plot_placeholder(f"VTK file not found:\n{vtk_path}")
            self.visualization_status.setText(f"VTK file not found: {vtk_path}")
            return

        try:
            if self.vtk_mode_combo.currentText() == "Cross-section":
                message = self.visualization_canvas.plot_cross_section(vtk_path, float(self.z_plane_spin.value()))
            else:
                message = self.visualization_canvas.plot_3d(vtk_path, int(self.max_trajectories_spin.value()))
        except Exception as exc:  # pragma: no cover - visualization depends on data quality
            QMessageBox.critical(self, "Visualization failed", str(exc))
            return

        self.visualization_status.setText(message)
        self.statusBar().showMessage(message)