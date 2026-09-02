"""VTK result browsing and plotting helpers for the NegAccel GUI."""

from __future__ import annotations

from pathlib import Path

from .common import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QSpinBox,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
    Qt,
    nested_get,
)
from .visualization import (
    NavigationToolbar2QT,
    OutputCanvas,
    TRAJECTORY_COLOR_MODE_ENERGY,
    TRAJECTORY_COLOR_MODE_SINGLE,
    TRAJECTORY_COLOR_MODE_SPECIES,
    TRAJECTORY_COLOR_MODE_STATUS,
    TrajectoryRenderOptions,
    configure_matplotlib_toolbar,
)
from negaccel_app.particles import (
    PARTICLE_FAMILIES,
    detect_particle_family_from_kinds,
    particle_label_from_export_id,
)
from negaccel_app.workflow.post_processing import (
    ARTIFACT_TYPE_PNG,
    ARTIFACT_TYPE_TEXT,
    ARTIFACT_TYPE_VTK,
    OUTPUT_CATEGORY_FINAL,
    OUTPUT_CATEGORY_GENERAL,
    OUTPUT_CATEGORY_ITERATION,
    build_output_catalog,
)


OUTPUT_CATEGORIES = [
    OUTPUT_CATEGORY_GENERAL,
    OUTPUT_CATEGORY_ITERATION,
    OUTPUT_CATEGORY_FINAL,
]
ARTIFACT_TYPES = [ARTIFACT_TYPE_VTK, ARTIFACT_TYPE_PNG, ARTIFACT_TYPE_TEXT]
TRAJECTORY_VTK_MODES = [
    "Cross-section",
    "3D trajectories",
    "Phase space x-x'",
    "Phase space y-y'",
]
GEOMETRY_VTK_MODES = ["Section view", "3D preview"]
FIELD_VTK_MODES = ["Scalar slice", "Line plot"]
SCALAR_FIELD_DISPLAY_MODES = ["Colormap plot", "Contour lines"]
TRAJECTORY_COLOR_MODES = [
    TRAJECTORY_COLOR_MODE_SINGLE,
    TRAJECTORY_COLOR_MODE_SPECIES,
    TRAJECTORY_COLOR_MODE_ENERGY,
    TRAJECTORY_COLOR_MODE_STATUS,
]
DEFAULT_TRAJECTORY_SPINBOX_MAX = 5000
_BEAM_SUMMARY_QUANTITIES = [
    "Beam current",
    "Beam centroids (x, y)",
    "Beam divergences (Dx, Dy)",
    "Average deviations (x', y')",
    "Average voltage",
    "Convergence",
    "Raw text preview",
]
_GRID_POWER_QUANTITIES = [
    "Grid power",
    "Grid current",
    "Grid net power",
    "Grid net current",
    "Grid particles",
    "Raw text preview",
]
# The breakdown file has no dedicated plotter yet, so it is offered as raw text; the
# quantity list is separate so a stacked species/origin plot can be added here later.
_GRID_POWER_BREAKDOWN_QUANTITIES = ["Raw text preview"]
_EMITTER_QUANTITIES = ["Emitter footprint", "Raw text preview"]
_RAW_TEXT_QUANTITIES = ["Raw text preview"]

TEXT_QUANTITY_MAP = {
    "Beam current": "beam_current",
    "Beam centroids (x, y)": "centroids_grouped",
    "Beam centroid x": "centroid_x",
    "Beam centroid y": "centroid_y",
    "Beam divergences (Dx, Dy)": "divergences_grouped",
    "Beam divergence Dx": "divergence_x",
    "Beam divergence Dy": "divergence_y",
    "Average deviations (x', y')": "deviations_grouped",
    "Average voltage": "average_voltage",
    "Convergence": "convergence",
    "Grid power": "power",
    "Grid current": "current",
    "Grid net power": "net_power",
    "Grid net current": "net_current",
    "Grid particles": "particles",
    "Emitter footprint": "emitter",
    "Raw text preview": "raw",
}


class ResultsMixin:
    def _ensure_result_widgets(self) -> None:
        if hasattr(self, "authoring_preview"):
            return

        self.authoring_preview = QPlainTextEdit()
        self.authoring_preview.setReadOnly(True)

        self.runtime_preview = QPlainTextEdit()
        self.runtime_preview.setReadOnly(True)

        self.process_output = QPlainTextEdit()
        self.process_output.setReadOnly(True)
        self.visualization_entries_by_path: dict[str, dict[str, object]] = {}
        self._default_scalar_slice_limits: tuple[float, float] | None = None
        self._default_scalar_slice_limits_path: str | None = None
        self._default_trajectory_energy_limits: tuple[float, float] | None = None
        self._default_trajectory_energy_limits_path: str | None = None
        self._trajectory_metadata_source_path: str | None = None
        self._trajectory_metadata_family: str = PARTICLE_FAMILIES[0]
        self._trajectory_species_checkboxes: dict[int, QCheckBox] = {}

        self.visualization_tab = QWidget()
        visualization_layout = QHBoxLayout(self.visualization_tab)

        controls_panel = QWidget()
        controls_panel.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        controls_panel.setMaximumWidth(360)
        controls_layout = QVBoxLayout(controls_panel)
        controls_layout.setContentsMargins(0, 0, 0, 0)

        self.output_category_combo = QComboBox()
        self.output_category_combo.addItems(OUTPUT_CATEGORIES)
        self.output_category_combo.setCurrentText(OUTPUT_CATEGORY_FINAL)
        self.output_category_combo.currentTextChanged.connect(self._on_visualization_category_changed)

        self.artifact_type_combo = QComboBox()
        self.artifact_type_combo.addItems(ARTIFACT_TYPES)
        self.artifact_type_combo.currentTextChanged.connect(self._on_visualization_artifact_changed)

        self.visualization_file_combo = QComboBox()
        self.visualization_file_combo.setMinimumWidth(240)
        self.visualization_file_combo.currentIndexChanged.connect(self._on_visualization_file_changed)

        self.vtk_mode_combo = QComboBox()
        self.vtk_mode_combo.addItems(TRAJECTORY_VTK_MODES)
        self.vtk_mode_combo.currentTextChanged.connect(self._refresh_visualization_options)

        self.quantity_combo = QComboBox()
        self.quantity_combo.setEnabled(False)

        self.z_plane_spin = QDoubleSpinBox()
        self.z_plane_spin.setRange(-10.0, 10.0)
        self.z_plane_spin.setDecimals(4)
        self.z_plane_spin.setSingleStep(0.001)
        self.z_plane_spin.setValue(0.565)

        self.scalar_min_spin = QDoubleSpinBox()
        self.scalar_min_spin.setRange(-1.0e30, 1.0e30)
        self.scalar_min_spin.setDecimals(12)
        self.scalar_min_spin.setSingleStep(0.001)
        self.scalar_min_spin.setKeyboardTracking(False)

        self.scalar_max_spin = QDoubleSpinBox()
        self.scalar_max_spin.setRange(-1.0e30, 1.0e30)
        self.scalar_max_spin.setDecimals(12)
        self.scalar_max_spin.setSingleStep(0.001)
        self.scalar_max_spin.setKeyboardTracking(False)

        self.scalar_display_mode_combo = QComboBox()
        self.scalar_display_mode_combo.addItems(SCALAR_FIELD_DISPLAY_MODES)
        self.scalar_display_mode_combo.currentTextChanged.connect(self._refresh_visualization_options)

        self.contour_lines_spin = QSpinBox()
        self.contour_lines_spin.setRange(2, 100)
        self.contour_lines_spin.setValue(10)

        self.scalar_range_reset_button = QPushButton("Reset scalar range")
        self.scalar_range_reset_button.clicked.connect(self._reset_scalar_slice_limits)
        
        self.plane_combo = QComboBox()
        self.plane_combo.addItems(["x", "y", "z"])
        self.plane_combo.setCurrentText("z")
        self.plane_combo.currentTextChanged.connect(self._refresh_visualization_options)
        
        self.show_geometry_checkbox = QCheckBox("Show geometry solids")
        self.show_geometry_checkbox.setChecked(False)
        
        self.show_trajectories_checkbox = QCheckBox("Show trajectories")
        self.show_trajectories_checkbox.setChecked(False)
        self.show_trajectories_checkbox.toggled.connect(self._refresh_visualization_options)
        
        self.trajectories_count_spin = QSpinBox()
        self.trajectories_count_spin.setRange(-1, 5000)
        self.trajectories_count_spin.setSpecialValueText("All")
        self.trajectories_count_spin.setValue(-1)
        self.trajectories_count_spin.setToolTip("-1: all trajectories, >0: random sample")

        self.trajectory_color_mode_combo = QComboBox()
        self.trajectory_color_mode_combo.addItems(TRAJECTORY_COLOR_MODES)

        self.trajectory_energy_min_spin = QDoubleSpinBox()
        self.trajectory_energy_min_spin.setRange(-1.0e30, 1.0e30)
        self.trajectory_energy_min_spin.setDecimals(6)
        self.trajectory_energy_min_spin.setSingleStep(1.0)
        self.trajectory_energy_min_spin.setKeyboardTracking(False)

        self.trajectory_energy_max_spin = QDoubleSpinBox()
        self.trajectory_energy_max_spin.setRange(-1.0e30, 1.0e30)
        self.trajectory_energy_max_spin.setDecimals(6)
        self.trajectory_energy_max_spin.setSingleStep(1.0)
        self.trajectory_energy_max_spin.setKeyboardTracking(False)

        self.trajectory_energy_reset_button = QPushButton("Reset energy range")
        self.trajectory_energy_reset_button.clicked.connect(self._reset_trajectory_energy_limits)

        self.trajectory_species_widget = QWidget()
        trajectory_species_layout = QVBoxLayout(self.trajectory_species_widget)
        trajectory_species_layout.setContentsMargins(0, 0, 0, 0)
        trajectory_species_layout.setSpacing(6)

        trajectory_species_buttons_row = QHBoxLayout()
        trajectory_species_buttons_row.setContentsMargins(0, 0, 0, 0)
        self.trajectory_species_all_button = QPushButton("All")
        self.trajectory_species_all_button.clicked.connect(lambda: self._set_all_trajectory_species(True))
        self.trajectory_species_none_button = QPushButton("None")
        self.trajectory_species_none_button.clicked.connect(lambda: self._set_all_trajectory_species(False))
        trajectory_species_buttons_row.addWidget(self.trajectory_species_all_button)
        trajectory_species_buttons_row.addWidget(self.trajectory_species_none_button)
        trajectory_species_buttons_row.addStretch(1)
        trajectory_species_layout.addLayout(trajectory_species_buttons_row)

        self.trajectory_species_list_widget = QWidget()
        self.trajectory_species_list_layout = QVBoxLayout(self.trajectory_species_list_widget)
        self.trajectory_species_list_layout.setContentsMargins(0, 0, 0, 0)
        self.trajectory_species_list_layout.setSpacing(3)
        trajectory_species_layout.addWidget(self.trajectory_species_list_widget)

        self.trajectory_species_empty_label = QLabel("No trajectory metadata available.")
        self.trajectory_species_empty_label.setWordWrap(True)
        self.trajectory_species_list_layout.addWidget(self.trajectory_species_empty_label)
        
        self.line_start_x_spin = QDoubleSpinBox()
        self.line_start_x_spin.setRange(-10.0, 10.0)
        self.line_start_x_spin.setDecimals(4)
        self.line_start_x_spin.setValue(0.0)
        
        self.line_start_y_spin = QDoubleSpinBox()
        self.line_start_y_spin.setRange(-10.0, 10.0)
        self.line_start_y_spin.setDecimals(4)
        self.line_start_y_spin.setValue(0.0)
        
        self.line_start_z_spin = QDoubleSpinBox()
        self.line_start_z_spin.setRange(-10.0, 10.0)
        self.line_start_z_spin.setDecimals(4)
        self.line_start_z_spin.setValue(0.01)
        
        self.line_dir_x_spin = QDoubleSpinBox()
        self.line_dir_x_spin.setRange(-1.0, 1.0)
        self.line_dir_x_spin.setDecimals(4)
        self.line_dir_x_spin.setValue(0.0)
        
        self.line_dir_y_spin = QDoubleSpinBox()
        self.line_dir_y_spin.setRange(-1.0, 1.0)
        self.line_dir_y_spin.setDecimals(4)
        self.line_dir_y_spin.setValue(0.0)
        
        self.line_dir_z_spin = QDoubleSpinBox()
        self.line_dir_z_spin.setRange(-1.0, 1.0)
        self.line_dir_z_spin.setDecimals(4)
        self.line_dir_z_spin.setValue(1.0)
        
        self.line_num_points_spin = QSpinBox()
        self.line_num_points_spin.setRange(2, 10000)
        self.line_num_points_spin.setValue(100)
        
        self.line_step_length_spin = QDoubleSpinBox()
        self.line_step_length_spin.setRange(0.0001, 1.0)
        self.line_step_length_spin.setDecimals(5)
        self.line_step_length_spin.setValue(0.001)
        
        self.max_trajectories_spin = QSpinBox()
        self.max_trajectories_spin.setRange(1, 5000)
        self.max_trajectories_spin.setValue(200)
        self.refresh_visualization_button = QPushButton("Refresh Files")
        self.refresh_visualization_button.clicked.connect(self.refresh_vtk_files)
        self.plot_button = QPushButton("Plot")
        self.plot_button.clicked.connect(self.plot_selected_output)

        controls_box = QGroupBox("Visualization controls")
        controls_box_layout = QVBoxLayout(controls_box)

        artifact_form = QFormLayout()
        artifact_form.addRow("Category", self.output_category_combo)
        artifact_form.addRow("Artifact", self.artifact_type_combo)
        artifact_form.addRow("Output file", self.visualization_file_combo)
        artifact_form.addRow("Quantity", self.quantity_combo)
        artifact_form.addRow("Display", self.vtk_mode_combo)
        controls_box_layout.addLayout(artifact_form)
        self.visualization_artifact_form = artifact_form

        options_form = QFormLayout()
        options_form.addRow("Slice plane", self.plane_combo)
        options_form.addRow("z-plane [m]", self.z_plane_spin)
        options_form.addRow("Scalar display", self.scalar_display_mode_combo)
        options_form.addRow("Contour lines", self.contour_lines_spin)
        options_form.addRow("Scalar min", self.scalar_min_spin)
        options_form.addRow("Scalar max", self.scalar_max_spin)
        options_form.addRow(self.scalar_range_reset_button)
        options_form.addRow(self.show_geometry_checkbox)
        options_form.addRow(self.show_trajectories_checkbox)
        options_form.addRow("Trajectories count", self.trajectories_count_spin)
        options_form.addRow("Trajectory colors", self.trajectory_color_mode_combo)
        options_form.addRow("Energy min [eV]", self.trajectory_energy_min_spin)
        options_form.addRow("Energy max [eV]", self.trajectory_energy_max_spin)
        options_form.addRow(self.trajectory_energy_reset_button)
        options_form.addRow("Particle species", self.trajectory_species_widget)
        
        # Line plot section
        self.line_start_widget = QWidget()
        line_start_layout = QHBoxLayout(self.line_start_widget)
        line_start_layout.setContentsMargins(0, 0, 0, 0)
        line_start_layout.addWidget(self.line_start_x_spin)
        line_start_layout.addWidget(self.line_start_y_spin)
        line_start_layout.addWidget(self.line_start_z_spin)
        options_form.addRow("Line start (x,y,z)", self.line_start_widget)
        
        self.line_dir_widget = QWidget()
        line_dir_layout = QHBoxLayout(self.line_dir_widget)
        line_dir_layout.setContentsMargins(0, 0, 0, 0)
        line_dir_layout.addWidget(self.line_dir_x_spin)
        line_dir_layout.addWidget(self.line_dir_y_spin)
        line_dir_layout.addWidget(self.line_dir_z_spin)
        options_form.addRow("Line direction (x,y,z)", self.line_dir_widget)
        
        options_form.addRow("Line num points", self.line_num_points_spin)
        options_form.addRow("Line step length [m]", self.line_step_length_spin)
        
        options_form.addRow("Max trajectories", self.max_trajectories_spin)
        controls_box_layout.addLayout(options_form)
        self.visualization_options_form = options_form

        actions_row = QHBoxLayout()
        actions_row.addWidget(self.refresh_visualization_button)
        actions_row.addWidget(self.plot_button)
        actions_row.addStretch(1)
        controls_box_layout.addLayout(actions_row)
        controls_layout.addWidget(controls_box)

        self.visualization_metadata = QLabel("Select an output artifact to populate the viewer.")
        self.visualization_metadata.setWordWrap(True)
        controls_layout.addWidget(self.visualization_metadata)
        controls_layout.addStretch(1)

        viewer_panel = QWidget()
        viewer_layout = QVBoxLayout(viewer_panel)
        viewer_layout.setContentsMargins(0, 0, 0, 0)

        self.visualization_stack = QStackedWidget()

        plot_page = QWidget()
        plot_page_layout = QVBoxLayout(plot_page)
        plot_page_layout.setContentsMargins(0, 0, 0, 0)
        self.visualization_canvas = OutputCanvas()
        self.visualization_toolbar = NavigationToolbar2QT(self.visualization_canvas, plot_page)
        configure_matplotlib_toolbar(self.visualization_toolbar)
        plot_page_layout.addWidget(self.visualization_toolbar)
        plot_page_layout.addWidget(self.visualization_canvas, 1)

        image_page = QWidget()
        image_page_layout = QVBoxLayout(image_page)
        image_page_layout.setContentsMargins(0, 0, 0, 0)
        self.image_scroll = self._build_image_viewer()
        image_page_layout.addWidget(self.image_scroll, 1)

        text_page = QWidget()
        text_page_layout = QVBoxLayout(text_page)
        text_page_layout.setContentsMargins(0, 0, 0, 0)
        self.visualization_text_preview = QPlainTextEdit()
        self.visualization_text_preview.setReadOnly(True)
        text_page_layout.addWidget(self.visualization_text_preview, 1)

        self.visualization_stack.addWidget(plot_page)
        self.visualization_stack.addWidget(image_page)
        self.visualization_stack.addWidget(text_page)
        viewer_layout.addWidget(self.visualization_stack, 1)

        self.visualization_status = QLabel("No VTK file selected yet.")
        viewer_layout.addWidget(self.visualization_status)

        self.visualization_controls_scroll_area = QScrollArea()
        self.visualization_controls_scroll_area.setWidgetResizable(True)
        self.visualization_controls_scroll_area.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self.visualization_controls_scroll_area.setSizePolicy(
            QSizePolicy.Policy.Preferred,
            QSizePolicy.Policy.Expanding,
        )
        self.visualization_controls_scroll_area.setMaximumWidth(380)
        self.visualization_controls_scroll_area.setMinimumHeight(0)
        self.visualization_controls_scroll_area.setWidget(controls_panel)

        visualization_layout.addWidget(self.visualization_controls_scroll_area)
        visualization_layout.addWidget(viewer_panel, 1)
        self._refresh_visualization_options()
        self.refresh_vtk_files()

    def _build_json_preview_panel(self) -> QWidget:
        self._ensure_result_widgets()
        if hasattr(self, "metadata_preview_panel"):
            return self.metadata_preview_panel

        def _wrap_preview(title: str, editor: QPlainTextEdit) -> QGroupBox:
            box = QGroupBox(title)
            layout = QVBoxLayout(box)
            layout.setContentsMargins(8, 8, 8, 8)
            layout.addWidget(editor)
            return box

        self.metadata_preview_panel = QSplitter(Qt.Orientation.Horizontal)
        self.metadata_preview_panel.addWidget(_wrap_preview("Authoring JSON", self.authoring_preview))
        self.metadata_preview_panel.addWidget(_wrap_preview("Generated Runtime JSON", self.runtime_preview))
        self.metadata_preview_panel.setStretchFactor(0, 1)
        self.metadata_preview_panel.setStretchFactor(1, 1)
        self.metadata_preview_panel.setSizes([520, 520])
        return self.metadata_preview_panel

    def _build_run_log_tab(self) -> QWidget:
        self._ensure_result_widgets()
        return self.process_output

    def _build_visualization_tab(self) -> QWidget:
        self._ensure_result_widgets()
        return self.visualization_tab

    def refresh_vtk_files(self) -> None:
        current_selection = self.visualization_file_combo.currentData()
        self.visualization_file_combo.blockSignals(True)
        self.visualization_file_combo.clear()
        entries = self._resolve_visualization_entries()
        self.visualization_entries_by_path = {
            str(entry["path"]): entry for entry in entries
        }
        if not entries:
            label = self._empty_visualization_message()
            self.visualization_file_combo.addItem(label, "")
            self.visualization_status.setText(label)
            self.visualization_metadata.setText(label)
            self.visualization_file_combo.blockSignals(False)
            self._on_visualization_file_changed()
            return

        selected_index = 0
        for index, entry in enumerate(entries):
            path = Path(str(entry["path"]))
            self.visualization_file_combo.addItem(path.name, str(path))
            if current_selection and str(path) == str(current_selection):
                selected_index = index
        self.visualization_file_combo.setCurrentIndex(selected_index)
        self.visualization_file_combo.blockSignals(False)
        artifact_type = self.artifact_type_combo.currentText()
        category = self.output_category_combo.currentText()
        self.visualization_status.setText(
            f"Found {len(entries)} {artifact_type.lower()} file(s) in {category.lower()}."
        )
        self._on_visualization_file_changed()

    def resolve_vtk_directory(self) -> Path | None:
        return self._resolve_output_directory("vtk", default_directory="VTK")

    def plot_selected_output(self) -> None:
        selected_path = self.visualization_file_combo.currentData()
        if not selected_path:
            self.visualization_canvas.plot_placeholder("No output artifact is currently available.")
            self.visualization_status.setText("No output artifact is currently available.")
            return

        output_path = Path(str(selected_path))
        if not output_path.exists():
            self.visualization_canvas.plot_placeholder(f"Output file not found:\n{output_path}")
            self.visualization_status.setText(f"Output file not found: {output_path}")
            return

        try:
            artifact_type = self.artifact_type_combo.currentText()
            entry = self.visualization_entries_by_path.get(str(output_path))
            trajectory_options = self._build_trajectory_render_options(output_path, entry)
            if artifact_type == ARTIFACT_TYPE_VTK:
                vtk_group = self._resolve_vtk_visualization_group(entry, output_path)
                vtk_mode = self.vtk_mode_combo.currentText()
                self.visualization_stack.setCurrentIndex(0)
                if vtk_group == "geometry":
                    runtime_case, _runtime_path = self._resolve_runtime_outputs()
                    if runtime_case is None:
                        raise RuntimeError("Runtime case is not available for geometry preview")
                    if vtk_mode == "3D preview":
                        message = self.visualization_canvas.plot_runtime_geometry_3d(runtime_case)
                    else:
                        message = self.visualization_canvas.plot_runtime_geometry_section(runtime_case)
                elif vtk_group == "field":
                    plane_axis = self.plane_combo.currentText()
                    plane_value = float(self.z_plane_spin.value())
                    show_geometry = self.show_geometry_checkbox.isChecked()
                    show_trajectories = self.show_trajectories_checkbox.isChecked()
                    trajectories_count = int(self.trajectories_count_spin.value())
                    
                    if vtk_mode == "Line plot":
                        start_point = (
                            float(self.line_start_x_spin.value()),
                            float(self.line_start_y_spin.value()),
                            float(self.line_start_z_spin.value()),
                        )
                        direction = (
                            float(self.line_dir_x_spin.value()),
                            float(self.line_dir_y_spin.value()),
                            float(self.line_dir_z_spin.value()),
                        )
                        num_points = int(self.line_num_points_spin.value())
                        step_length = float(self.line_step_length_spin.value())
                        message = self.visualization_canvas.plot_vtk_line_plot(
                            output_path,
                            start_point,
                            direction,
                            num_points,
                            step_length,
                        )
                    else:
                        message = self.visualization_canvas.plot_vtk_structured_slice(
                            output_path,
                            plane_value,
                            plane_axis=plane_axis,
                            show_geometry=show_geometry,
                            show_trajectories=show_trajectories,
                            trajectories_count=trajectories_count,
                            scalar_min=float(self.scalar_min_spin.value()),
                            scalar_max=float(self.scalar_max_spin.value()),
                            scalar_display_mode=self.scalar_display_mode_combo.currentText(),
                            contour_lines=int(self.contour_lines_spin.value()),
                            trajectory_options=trajectory_options,
                        )
                else:
                    if vtk_mode == "3D trajectories":
                        message = self.visualization_canvas.plot_vtk_3d(
                            output_path,
                            int(self.max_trajectories_spin.value()),
                            trajectory_options,
                        )
                    elif vtk_mode == "Phase space x-x'":
                        message = self.visualization_canvas.plot_vtk_phase_space(
                            output_path,
                            float(self.z_plane_spin.value()),
                            transverse_axis="x",
                            trajectory_options=trajectory_options,
                        )
                    elif vtk_mode == "Phase space y-y'":
                        message = self.visualization_canvas.plot_vtk_phase_space(
                            output_path,
                            float(self.z_plane_spin.value()),
                            transverse_axis="y",
                            trajectory_options=trajectory_options,
                        )
                    else:
                        message = self.visualization_canvas.plot_vtk_cross_section(
                            output_path,
                            float(self.z_plane_spin.value()),
                            trajectory_options,
                        )
            elif artifact_type == "PNG plots":
                message = self._show_png(output_path)
            else:
                message = self._show_text_output(output_path)
        except Exception as exc:  # pragma: no cover - visualization depends on data quality
            QMessageBox.critical(self, "Visualization failed", str(exc))
            return

        self.visualization_status.setText(message)
        self.statusBar().showMessage(message)

    def _build_image_viewer(self):
        image_label = QLabel("Select a PNG plot to visualize.")
        image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        image_label.setWordWrap(True)
        image_label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        image_label.setMinimumSize(320, 240)
        container = QWidget()
        container_layout = QVBoxLayout(container)
        container_layout.setContentsMargins(0, 0, 0, 0)
        container_layout.addWidget(image_label, 1)
        preview = QPlainTextEdit()
        preview.setReadOnly(True)
        preview.setMaximumHeight(0)
        container_layout.addWidget(preview)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(container)
        self.visualization_image_label = image_label
        self.visualization_image_preview = preview
        return scroll

    def _resolve_runtime_outputs(self):
        if self.current_runtime_case is None:
            try:
                runtime_case, runtime_path = self.build_runtime_case_preview()
            except Exception:
                return None, None
        else:
            runtime_case = self.current_runtime_case
            runtime_path = self.current_runtime_path
        return runtime_case, runtime_path

    def _resolve_output_directory(self, output_key: str, default_directory: str) -> Path | None:
        runtime_case, runtime_path = self._resolve_runtime_outputs()
        if runtime_case is None or runtime_path is None:
            return None

        outputs = runtime_case.get("outputs", {})
        output_config = outputs.get(output_key, {})
        directory_name = output_config.get("directory", default_directory) if isinstance(output_config, dict) else default_directory
        root_directory = outputs.get("rootDirectory") or str(runtime_path.parent)

        path = Path(directory_name)
        if not path.is_absolute():
            path = Path(root_directory) / path
        return path

    def _resolve_visualization_catalog(self) -> list[dict[str, object]]:
        runtime_case, runtime_path = self._resolve_runtime_outputs()
        if runtime_case is None or runtime_path is None:
            return []
        return build_output_catalog(runtime_case, runtime_path)

    def _resolve_visualization_entries(self) -> list[dict[str, object]]:
        category = self.output_category_combo.currentText()
        artifact_type = self.artifact_type_combo.currentText()
        return [
            entry
            for entry in self._resolve_visualization_catalog()
            if entry.get("category") == category and entry.get("artifactType") == artifact_type
        ]

    def _available_artifact_types_for_category(self, category: str) -> list[str]:
        available = {
            str(entry.get("artifactType"))
            for entry in self._resolve_visualization_catalog()
            if entry.get("category") == category
        }
        if not available:
            return list(ARTIFACT_TYPES)
        return [artifact_type for artifact_type in ARTIFACT_TYPES if artifact_type in available]

    def _refresh_artifact_type_options(self) -> None:
        current_text = self.artifact_type_combo.currentText()
        available = self._available_artifact_types_for_category(self.output_category_combo.currentText())
        self.artifact_type_combo.blockSignals(True)
        self.artifact_type_combo.clear()
        self.artifact_type_combo.addItems(available)
        if current_text in available:
            self.artifact_type_combo.setCurrentText(current_text)
        elif available:
            self.artifact_type_combo.setCurrentIndex(0)
        self.artifact_type_combo.blockSignals(False)

    def _empty_visualization_message(self) -> str:
        category = self.output_category_combo.currentText()
        artifact_type = self.artifact_type_combo.currentText()
        return f"No {artifact_type.lower()} available in {category.lower()} yet."

    def _on_visualization_category_changed(self) -> None:
        self._refresh_artifact_type_options()
        self._on_visualization_controls_changed()

    def _on_visualization_artifact_changed(self) -> None:
        self._on_visualization_controls_changed()

    def _on_visualization_controls_changed(self) -> None:
        self._refresh_visualization_options()
        self.refresh_vtk_files()

    def _refresh_visualization_options(self) -> None:
        entry = self._selected_visualization_entry()
        self._sync_vtk_mode_options(entry)
        artifact_type = self.artifact_type_combo.currentText()
        is_vtk = artifact_type == "VTK"
        is_text = artifact_type == "Text outputs"
        vtk_group = self._resolve_vtk_visualization_group(entry)
        vtk_mode = self.vtk_mode_combo.currentText()
        is_field_vtk = is_vtk and vtk_group == "field"
        is_line_plot = is_field_vtk and vtk_mode == "Line plot"
        is_scalar_field_vtk = is_field_vtk and not is_line_plot
        needs_trajectory_render_controls = is_vtk and (
            vtk_group == "trajectory" or (is_scalar_field_vtk and self.show_trajectories_checkbox.isChecked())
        )
        needs_z_plane = is_vtk and (
            (is_field_vtk and not is_line_plot)
            or vtk_mode in {"Cross-section", "Phase space x-x'", "Phase space y-y'"}
        )
        needs_max_trajectories = is_vtk and vtk_group == "trajectory" and vtk_mode == "3D trajectories"
        needs_contour_lines = (
            is_scalar_field_vtk and self.scalar_display_mode_combo.currentText() == "Contour lines"
        )
        plane_value_label = f"{self.plane_combo.currentText()}-plane [m]" if is_field_vtk and not is_line_plot else "z-plane [m]"
        
        self.vtk_mode_combo.setVisible(is_vtk)
        self.visualization_artifact_form.labelForField(self.vtk_mode_combo).setVisible(is_vtk)
        self.quantity_combo.setVisible(is_text)
        self.visualization_artifact_form.labelForField(self.quantity_combo).setVisible(is_text)
        self._set_form_row_label(self.visualization_options_form, self.z_plane_spin, plane_value_label)
        self._set_form_row_visible(self.visualization_options_form, self.z_plane_spin, needs_z_plane)
        self._set_form_row_visible(
            self.visualization_options_form,
            self.max_trajectories_spin,
            needs_max_trajectories,
        )
        
        # Field VTK controls visibility
        self._set_form_row_visible(self.visualization_options_form, self.plane_combo, is_field_vtk and not is_line_plot)
        self._set_form_row_visible(self.visualization_options_form, self.scalar_display_mode_combo, is_scalar_field_vtk)
        self._set_form_row_visible(self.visualization_options_form, self.contour_lines_spin, needs_contour_lines)
        self._set_form_row_visible(self.visualization_options_form, self.scalar_min_spin, is_scalar_field_vtk)
        self._set_form_row_visible(self.visualization_options_form, self.scalar_max_spin, is_scalar_field_vtk)
        self._set_form_row_visible(self.visualization_options_form, self.scalar_range_reset_button, is_scalar_field_vtk)
        self.show_geometry_checkbox.setVisible(is_scalar_field_vtk)
        self.show_trajectories_checkbox.setVisible(is_scalar_field_vtk)
        self._set_form_row_visible(
            self.visualization_options_form,
            self.trajectories_count_spin,
            is_scalar_field_vtk and self.show_trajectories_checkbox.isChecked()
        )
        self._set_form_row_visible(
            self.visualization_options_form,
            self.trajectory_color_mode_combo,
            needs_trajectory_render_controls,
        )
        self._set_form_row_visible(
            self.visualization_options_form,
            self.trajectory_energy_min_spin,
            needs_trajectory_render_controls,
        )
        self._set_form_row_visible(
            self.visualization_options_form,
            self.trajectory_energy_max_spin,
            needs_trajectory_render_controls,
        )
        self._set_form_row_visible(
            self.visualization_options_form,
            self.trajectory_energy_reset_button,
            needs_trajectory_render_controls,
        )
        self._set_form_row_visible(
            self.visualization_options_form,
            self.trajectory_species_widget,
            needs_trajectory_render_controls,
        )
        
        # Line plot controls visibility
        self._set_form_row_visible(self.visualization_options_form, self.line_start_widget, is_line_plot)
        self._set_form_row_visible(self.visualization_options_form, self.line_dir_widget, is_line_plot)
        self._set_form_row_visible(self.visualization_options_form, self.line_num_points_spin, is_line_plot)
        self._set_form_row_visible(self.visualization_options_form, self.line_step_length_spin, is_line_plot)

        self._refresh_scalar_slice_limits()
        
        if not is_text and self.visualization_stack.currentIndex() == 2:
            self.visualization_text_preview.clear()

    def _set_form_row_visible(self, form: QFormLayout, field: QWidget, visible: bool) -> None:
        field.setVisible(visible)
        label = form.labelForField(field)
        if label is not None:
            label.setVisible(visible)

    def _set_form_row_label(self, form: QFormLayout, field: QWidget, text: str) -> None:
        label = form.labelForField(field)
        if isinstance(label, QLabel):
            label.setText(text)

    def _sync_vtk_mode_options(self, entry: dict[str, object] | None) -> None:
        if self.artifact_type_combo.currentText() != ARTIFACT_TYPE_VTK:
            desired_items = TRAJECTORY_VTK_MODES
        else:
            vtk_group = self._resolve_vtk_visualization_group(entry)
            if vtk_group == "geometry":
                desired_items = GEOMETRY_VTK_MODES
            elif vtk_group == "field":
                desired_items = FIELD_VTK_MODES
            else:
                desired_items = TRAJECTORY_VTK_MODES
        current_items = [self.vtk_mode_combo.itemText(index) for index in range(self.vtk_mode_combo.count())]
        if current_items == desired_items:
            return

        previous_text = self.vtk_mode_combo.currentText()
        self.vtk_mode_combo.blockSignals(True)
        self.vtk_mode_combo.clear()
        self.vtk_mode_combo.addItems(desired_items)
        if previous_text in desired_items:
            self.vtk_mode_combo.setCurrentText(previous_text)
        self.vtk_mode_combo.blockSignals(False)

    def _selected_visualization_entry(self) -> dict[str, object] | None:
        selected_path = self.visualization_file_combo.currentData()
        if not selected_path:
            return None
        return self.visualization_entries_by_path.get(str(selected_path))

    def _scalar_field_visualization_selected(self) -> bool:
        entry = self._selected_visualization_entry()
        return (
            self.artifact_type_combo.currentText() == ARTIFACT_TYPE_VTK
            and self._resolve_vtk_visualization_group(entry) == "field"
            and self.vtk_mode_combo.currentText() != "Line plot"
        )

    def _set_scalar_slice_limit_values(self, scalar_min: float, scalar_max: float) -> None:
        self.scalar_min_spin.blockSignals(True)
        self.scalar_max_spin.blockSignals(True)
        self.scalar_min_spin.setValue(float(scalar_min))
        self.scalar_max_spin.setValue(float(scalar_max))
        self.scalar_min_spin.blockSignals(False)
        self.scalar_max_spin.blockSignals(False)

    def _set_trajectory_spinbox_limit(
        self,
        spin_box: QSpinBox,
        available_count: int | None,
        *,
        allow_all: bool,
    ) -> None:
        maximum = DEFAULT_TRAJECTORY_SPINBOX_MAX if available_count is None else max(1, int(available_count))
        minimum = -1 if allow_all else 1
        current_value = int(spin_box.value())

        spin_box.blockSignals(True)
        spin_box.setRange(minimum, maximum)
        if allow_all and current_value <= -1:
            spin_box.setValue(-1)
        else:
            spin_box.setValue(min(max(current_value, minimum), maximum))
        spin_box.blockSignals(False)

    def _refresh_trajectory_spin_limits(
        self,
        output_path: Path | None,
        entry: dict[str, object] | None,
    ) -> None:
        available_count: int | None = None
        if output_path is not None and self.artifact_type_combo.currentText() == ARTIFACT_TYPE_VTK:
            vtk_group = self._resolve_vtk_visualization_group(entry, output_path)
            try:
                if vtk_group == "field":
                    available_count = self.visualization_canvas.get_vtk_overlay_trajectory_count(output_path)
                elif vtk_group == "trajectory":
                    available_count = self.visualization_canvas.get_vtk_trajectory_count(output_path)
            except Exception:
                available_count = None

        self._set_trajectory_spinbox_limit(
            self.trajectories_count_spin,
            available_count,
            allow_all=True,
        )
        self._set_trajectory_spinbox_limit(
            self.max_trajectories_spin,
            available_count,
            allow_all=False,
        )

    def _refresh_scalar_slice_limits(self, reset_values: bool = False) -> None:
        selected_path = self.visualization_file_combo.currentData()
        if not selected_path or not self._scalar_field_visualization_selected():
            self._default_scalar_slice_limits = None
            self._default_scalar_slice_limits_path = None
            return

        output_path = Path(str(selected_path))
        path_key = str(output_path)
        if (
            not reset_values
            and self._default_scalar_slice_limits is not None
            and self._default_scalar_slice_limits_path == path_key
        ):
            return

        scalar_min, scalar_max = self.visualization_canvas.get_vtk_structured_scalar_range(output_path)
        self._default_scalar_slice_limits = (scalar_min, scalar_max)
        self._default_scalar_slice_limits_path = path_key
        self._set_scalar_slice_limit_values(scalar_min, scalar_max)

    def _reset_scalar_slice_limits(self) -> None:
        if not self._scalar_field_visualization_selected():
            return

        self._refresh_scalar_slice_limits(reset_values=True)
        self.plot_selected_output()

    def _resolve_runtime_particle_family(self) -> str:
        runtime_case, _runtime_path = self._resolve_runtime_outputs()
        if not isinstance(runtime_case, dict):
            return PARTICLE_FAMILIES[0]

        particles = runtime_case.get("particles")
        if not isinstance(particles, dict):
            return PARTICLE_FAMILIES[0]

        kinds: list[str] = []
        raw_types = particles.get("types") if isinstance(particles.get("types"), list) else []
        type_lookup: dict[str, str] = {}
        for raw_type in raw_types:
            if not isinstance(raw_type, dict):
                continue
            kind = str(raw_type.get("kind") or "").strip()
            type_id = str(raw_type.get("id") or "").strip()
            if kind:
                kinds.append(kind)
            if kind and type_id:
                type_lookup[type_id] = kind

        raw_sources = particles.get("sources") if isinstance(particles.get("sources"), list) else []
        for raw_source in raw_sources:
            if not isinstance(raw_source, dict):
                continue
            source_kind = str(raw_source.get("kind") or "").strip()
            if source_kind:
                kinds.append(source_kind)
                continue
            particle_type_id = str(raw_source.get("particleTypeId") or "").strip()
            mapped_kind = type_lookup.get(particle_type_id)
            if mapped_kind:
                kinds.append(mapped_kind)

        return detect_particle_family_from_kinds(kinds)

    def _clear_trajectory_species_checkboxes(self) -> None:
        for checkbox in self._trajectory_species_checkboxes.values():
            self.trajectory_species_list_layout.removeWidget(checkbox)
            checkbox.deleteLater()
        self._trajectory_species_checkboxes = {}

    def _selected_trajectory_particle_kinds(self) -> tuple[int, ...] | None:
        if not self._trajectory_species_checkboxes:
            return None
        return tuple(
            sorted(
                kind_id
                for kind_id, checkbox in self._trajectory_species_checkboxes.items()
                if checkbox.isChecked()
            )
        )

    def _populate_trajectory_species_checkboxes(
        self,
        particle_kind_ids: tuple[int, ...],
        family: str,
        selected_kind_ids: set[int],
    ) -> None:
        self._clear_trajectory_species_checkboxes()

        ordered_kind_ids = tuple(sorted(int(kind_id) for kind_id in particle_kind_ids))
        for kind_id in ordered_kind_ids:
            checkbox = QCheckBox(particle_label_from_export_id(kind_id, family))
            checkbox.setChecked(kind_id in selected_kind_ids)
            self.trajectory_species_list_layout.insertWidget(
                self.trajectory_species_list_layout.count() - 1,
                checkbox,
            )
            self._trajectory_species_checkboxes[kind_id] = checkbox

        has_species = bool(ordered_kind_ids)
        self.trajectory_species_empty_label.setVisible(not has_species)
        if not has_species:
            self.trajectory_species_empty_label.setText("No trajectory metadata available.")
        self.trajectory_species_all_button.setEnabled(has_species)
        self.trajectory_species_none_button.setEnabled(has_species)

    def _set_all_trajectory_species(self, checked: bool) -> None:
        for checkbox in self._trajectory_species_checkboxes.values():
            checkbox.blockSignals(True)
            checkbox.setChecked(checked)
            checkbox.blockSignals(False)

    def _set_trajectory_energy_limit_values(self, energy_min_ev: float, energy_max_ev: float) -> None:
        self.trajectory_energy_min_spin.blockSignals(True)
        self.trajectory_energy_max_spin.blockSignals(True)
        self.trajectory_energy_min_spin.setValue(float(energy_min_ev))
        self.trajectory_energy_max_spin.setValue(float(energy_max_ev))
        self.trajectory_energy_min_spin.blockSignals(False)
        self.trajectory_energy_max_spin.blockSignals(False)

    def _refresh_trajectory_filter_controls(self, reset_values: bool = False) -> None:
        selected_path = self.visualization_file_combo.currentData()
        if not selected_path or self.artifact_type_combo.currentText() != ARTIFACT_TYPE_VTK:
            self._trajectory_metadata_source_path = None
            self._default_trajectory_energy_limits = None
            self._default_trajectory_energy_limits_path = None
            self._populate_trajectory_species_checkboxes((), PARTICLE_FAMILIES[0], set())
            return

        output_path = Path(str(selected_path))
        entry = self.visualization_entries_by_path.get(str(output_path))
        vtk_group = self._resolve_vtk_visualization_group(entry, output_path)
        if vtk_group not in {"trajectory", "field"}:
            self._trajectory_metadata_source_path = None
            self._default_trajectory_energy_limits = None
            self._default_trajectory_energy_limits_path = None
            self._populate_trajectory_species_checkboxes((), PARTICLE_FAMILIES[0], set())
            return

        family = self._resolve_runtime_particle_family()
        metadata = self.visualization_canvas.get_vtk_trajectory_metadata(output_path, family)
        if metadata is None:
            self._trajectory_metadata_source_path = None
            self._trajectory_metadata_family = family
            self._default_trajectory_energy_limits = None
            self._default_trajectory_energy_limits_path = None
            self._populate_trajectory_species_checkboxes((), family, set())
            return

        source_path_key = str(metadata.source_path)
        preserve_selection = source_path_key == self._trajectory_metadata_source_path and not reset_values
        current_selection = set(self._selected_trajectory_particle_kinds() or ()) if preserve_selection else set()
        available_kind_ids = set(int(kind_id) for kind_id in metadata.available_particle_kinds)
        if preserve_selection:
            selected_kind_ids = {kind_id for kind_id in current_selection if kind_id in available_kind_ids}
        else:
            selected_kind_ids = set(available_kind_ids)

        self._populate_trajectory_species_checkboxes(metadata.available_particle_kinds, family, selected_kind_ids)
        self._trajectory_metadata_source_path = source_path_key
        self._trajectory_metadata_family = family

        if metadata.energy_min_ev is not None and metadata.energy_max_ev is not None:
            self._default_trajectory_energy_limits = (metadata.energy_min_ev, metadata.energy_max_ev)
            self._default_trajectory_energy_limits_path = source_path_key
            if not preserve_selection or reset_values:
                self._set_trajectory_energy_limit_values(metadata.energy_min_ev, metadata.energy_max_ev)
        else:
            self._default_trajectory_energy_limits = None
            self._default_trajectory_energy_limits_path = None

    def _reset_trajectory_energy_limits(self) -> None:
        if self._default_trajectory_energy_limits is None:
            return
        self._set_trajectory_energy_limit_values(*self._default_trajectory_energy_limits)
        self.plot_selected_output()

    def _build_trajectory_render_options(
        self,
        output_path: Path,
        entry: dict[str, object] | None,
    ) -> TrajectoryRenderOptions | None:
        if self.artifact_type_combo.currentText() != ARTIFACT_TYPE_VTK:
            return None
        vtk_group = self._resolve_vtk_visualization_group(entry, output_path)
        if vtk_group not in {"trajectory", "field"}:
            return None

        energy_min_ev: float | None = None
        energy_max_ev: float | None = None
        if self._default_trajectory_energy_limits is not None:
            energy_min_ev = float(self.trajectory_energy_min_spin.value())
            energy_max_ev = float(self.trajectory_energy_max_spin.value())
            if energy_max_ev < energy_min_ev:
                raise ValueError("Trajectory energy maximum must be greater than or equal to energy minimum")

        return TrajectoryRenderOptions(
            color_mode=self.trajectory_color_mode_combo.currentText(),
            selected_particle_kinds=self._selected_trajectory_particle_kinds(),
            energy_min_ev=energy_min_ev,
            energy_max_ev=energy_max_ev,
            particle_family=self._trajectory_metadata_family,
        )

    def _resolve_vtk_visualization_group(
        self,
        entry: dict[str, object] | None,
        output_path: Path | None = None,
    ) -> str:
        if not isinstance(entry, dict):
            return "trajectory"

        kind = str(entry.get("kind") or "").lower()
        name = str(entry.get("name") or (output_path.name if output_path is not None else "")).lower()

        if kind == "geometry-vtk":
            return "geometry"
        if kind == "trajectories-vtk":
            return "trajectory"
        if kind in {"potential-vtk", "space-charge-vtk", "simulation-state-vtk"}:
            return "field"
        if kind == "iteration-artifact":
            if name.endswith("_trajectories.vtk"):
                return "trajectory"
            if (
                name.endswith("_potential.vtk")
                or name.endswith("_scharge.vtk")
                or name.endswith("_simulation_state.vtk")
            ):
                return "field"
            return "trajectory"
        if kind == "vtk":
            if (
                name.endswith("_potential.vtk")
                or name.endswith("_scharge.vtk")
                or name.endswith("_simulation_state.vtk")
            ):
                return "field"
            return "trajectory"
        return "trajectory"

    def _quantity_options_for_entry(self, entry: dict[str, object] | None) -> list[str]:
        if entry is None:
            return []
        kind = str(entry.get("kind") or "raw")
        if kind in {"iteration-artifact", "diagnostic-summary"}:
            return list(_BEAM_SUMMARY_QUANTITIES)
        if kind == "grid-power-summary":
            return list(_GRID_POWER_QUANTITIES)
        if kind == "grid-power-breakdown":
            return list(_GRID_POWER_BREAKDOWN_QUANTITIES)
        if kind == "emitter-map":
            return list(_EMITTER_QUANTITIES)
        return list(_RAW_TEXT_QUANTITIES)

    def _refresh_quantity_options(self, entry: dict[str, object] | None = None) -> None:
        if self.artifact_type_combo.currentText() != ARTIFACT_TYPE_TEXT:
            self.quantity_combo.blockSignals(True)
            self.quantity_combo.clear()
            self.quantity_combo.setEnabled(False)
            self.quantity_combo.blockSignals(False)
            return

        current_text = self.quantity_combo.currentText()
        options = self._quantity_options_for_entry(entry)
        self.quantity_combo.blockSignals(True)
        self.quantity_combo.clear()
        self.quantity_combo.addItems(options)
        if current_text in options:
            self.quantity_combo.setCurrentText(current_text)
        elif options:
            self.quantity_combo.setCurrentIndex(0)
        self.quantity_combo.setEnabled(bool(options))
        self.quantity_combo.blockSignals(False)

    def _on_visualization_file_changed(self) -> None:
        selected_path = self.visualization_file_combo.currentData()
        if not selected_path:
            self._refresh_quantity_options(None)
            self._refresh_trajectory_spin_limits(None, None)
            self._refresh_trajectory_filter_controls(reset_values=True)
            self._refresh_visualization_options()
            self.visualization_metadata.setText(self._empty_visualization_message())
            return
        output_path = Path(str(selected_path))
        entry = self.visualization_entries_by_path.get(str(output_path))
        self._refresh_quantity_options(entry)
        self._refresh_trajectory_spin_limits(output_path, entry)
        self._refresh_trajectory_filter_controls(reset_values=True)
        self._refresh_visualization_options()
        self._refresh_scalar_slice_limits(reset_values=True)
        if entry is None:
            self.visualization_metadata.setText(str(output_path))
            return

        lines = [
            f"Category: {entry.get('category', 'Unknown')}",
            f"Artifact type: {entry.get('artifactType', 'Unknown')}",
            f"Kind: {entry.get('kind', 'raw')}",
        ]
        if entry.get("iteration") is not None:
            lines.append(f"Iteration: {entry['iteration']}")
        lines.append("")
        lines.append(str(output_path))
        self.visualization_metadata.setText("\n".join(lines))

    def _show_png(self, png_path: Path) -> str:
        from PySide6.QtGui import QPixmap

        pixmap = QPixmap(str(png_path))
        if pixmap.isNull():
            raise RuntimeError(f"Could not load image: {png_path}")
        self.visualization_stack.setCurrentIndex(1)
        target_size = self.visualization_image_label.size()
        scaled = pixmap.scaled(
            target_size if target_size.width() > 1 and target_size.height() > 1 else pixmap.size(),
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        self.visualization_image_label.setPixmap(scaled)
        self.visualization_image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        return f"Displayed {png_path.name}"

    def _show_text_output(self, text_path: Path) -> str:
        entry = self.visualization_entries_by_path.get(str(text_path))
        kind = str(entry.get("kind") or "raw") if entry is not None else "raw"
        quantity = self.quantity_combo.currentText() or "Raw text preview"
        metric = TEXT_QUANTITY_MAP.get(quantity, "raw")

        if metric == "raw":
            return self._show_raw_text(text_path)

        if kind == "iteration-artifact" and metric in {
            "beam_current",
            "centroid_x",
            "centroid_y",
            "divergence_x",
            "divergence_y",
            "average_voltage",
            "centroids_grouped",
            "divergences_grouped",
            "deviations_grouped",
            "convergence",
        }:
            self.visualization_stack.setCurrentIndex(0)
            return self.visualization_canvas.plot_iteration_summary(text_path, metric)

        if kind == "grid-power-summary" and metric in {"power", "current", "particles"}:
            self.visualization_stack.setCurrentIndex(0)
            return self.visualization_canvas.plot_grid_power_summary(text_path, metric)

        if kind == "emitter-map" and metric == "emitter":
            self.visualization_stack.setCurrentIndex(0)
            return self.visualization_canvas.plot_emitter_map(text_path)

        if kind == "diagnostic-summary" and metric in {
            "beam_current",
            "centroid_x",
            "centroid_y",
            "divergence_x",
            "divergence_y",
            "average_voltage",
            "centroids_grouped",
            "divergences_grouped",
            "deviations_grouped",
        }:
            self.visualization_stack.setCurrentIndex(0)
            return self.visualization_canvas.plot_diagnostic_summary(text_path, metric)

        return self._show_raw_text(text_path)

    def _show_raw_text(self, text_path: Path) -> str:
        self.visualization_stack.setCurrentIndex(2)
        self.visualization_text_preview.setPlainText(text_path.read_text(encoding="utf-8", errors="replace"))
        return f"Loaded raw text preview for {text_path.name}"