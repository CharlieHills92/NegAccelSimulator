"""Magnetic-field workspace section for the NegAccel GUI."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)

from ..common import ParameterizedEditor, REPO_ROOT, WorkflowError, nested_get, parse_number_list
from ..visualization import configure_matplotlib_canvas, configure_matplotlib_toolbar


def _merge_directory_file_path(directory: str, file_path: str) -> str:
    normalized_path = file_path.strip()
    normalized_directory = directory.strip()
    if not normalized_path or not normalized_directory:
        return normalized_path

    candidate = Path(normalized_path)
    if candidate.is_absolute():
        return candidate.as_posix()

    directory_path = Path(normalized_directory)
    try:
        candidate.relative_to(directory_path)
        return candidate.as_posix()
    except ValueError:
        return (directory_path / candidate).as_posix()


def _double_spin(minimum: float, maximum: float, decimals: int, step: float) -> QDoubleSpinBox:
    spin = QDoubleSpinBox()
    spin.setRange(minimum, maximum)
    spin.setDecimals(decimals)
    spin.setSingleStep(step)
    return spin


class MagneticProfileCanvas(FigureCanvasQTAgg):
    def __init__(self) -> None:
        self.figure = Figure(figsize=(6.8, 5.2), tight_layout=True)
        super().__init__(self.figure)
        configure_matplotlib_canvas(self)
        self.plot_placeholder("Configure magnetic sources and profile settings to preview B components.")

    def plot_placeholder(self, message: str) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.axis("off")
        axis.text(0.5, 0.5, message, ha="center", va="center", fontsize=11)
        self.draw_idle()

    def plot_components(self, s_values: np.ndarray, components: dict[str, np.ndarray], selected: list[str]) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        for name in selected:
            axis.plot(s_values, components[name], label=name, linewidth=1.8)
        axis.set_xlabel("Path length s [m]")
        axis.set_ylabel("B component [T]")
        axis.set_title("Magnetic profile along selected line")
        axis.grid(True, alpha=0.3)
        axis.legend()
        self.draw_idle()


class _FieldPayloadWidget(QWidget):
    def __init__(self, window, source_type: str = "constant", value: str = "", changed_callback=None) -> None:
        super().__init__()
        self._window = window
        self._changed_callback = changed_callback

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        self._edit = QLineEdit(value)
        self._browse = QPushButton("Browse")
        layout.addWidget(self._edit, 1)
        layout.addWidget(self._browse)

        self._edit.textChanged.connect(self._notify_change)
        self._browse.clicked.connect(self._browse_file)
        self.set_source_type(source_type)

    def _notify_change(self) -> None:
        if self._changed_callback is not None:
            self._changed_callback()

    def _browse_file(self) -> None:
        current_value = self._edit.text().strip()
        start_dir = ""
        if current_value:
            current_path = Path(current_value).expanduser()
            if current_path.is_file():
                start_dir = str(current_path.parent)
            else:
                start_dir = str(current_path.parent)
        elif self._window.authoring_path is not None:
            start_dir = str(self._window.authoring_path.parent)
        else:
            start_dir = str(REPO_ROOT)

        selected, _ = QFileDialog.getOpenFileName(
            self,
            "Select magnetic field file",
            start_dir,
            "Magnetic field files (*.fld)",
        )
        if selected:
            self._edit.setText(selected)

    def set_source_type(self, source_type: str) -> None:
        normalized = source_type.strip().lower()
        is_file = normalized == "file"
        self._browse.setVisible(is_file)
        self._browse.setEnabled(is_file)
        if is_file:
            self._edit.setPlaceholderText("file path (*.fld)")
        else:
            self._edit.setPlaceholderText("bx, by, bz")

    def text(self) -> str:
        return self._edit.text()

    def setText(self, value: str) -> None:
        self._edit.setText(value)


class _MagneticWorkspaceWidget(QWidget):
    def __init__(self, window) -> None:
        super().__init__()
        self._window = window
        self._field_cache: dict[str, dict[str, Any]] = {}
        self._build_ui()

    def _row_for_widget(self, widget: QWidget, column: int) -> int:
        for row in range(self._table.rowCount()):
            cell_widget = self._table.cellWidget(row, column)
            if cell_widget is widget:
                return row
            if isinstance(cell_widget, ParameterizedEditor) and cell_widget.editor_widget() is widget:
                return row
        return 0

    def _unwrap_editor(self, widget: QWidget | None, expected_type):
        if isinstance(widget, expected_type):
            return widget
        if isinstance(widget, ParameterizedEditor):
            editor = widget.editor_widget()
            if isinstance(editor, expected_type):
                return editor
        return None

    def _build_ui(self) -> None:
        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        root.addWidget(splitter)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(6, 6, 6, 6)

        fields_box = QGroupBox("Magnetic field sources")
        fields_layout = QVBoxLayout(fields_box)

        self._enabled = QCheckBox("Enable external magnetic field")

        header_form = QFormLayout()
        header_form.addRow(self._enabled)
        fields_layout.addLayout(header_form)

        self._table = QTableWidget(0, 4)
        self._table.setHorizontalHeaderLabels(["Name", "Type", "Scale", "Value or file"])  # constant uses bx,by,bz
        self._table.verticalHeader().setVisible(False)
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.setMinimumHeight(280)
        fields_layout.addWidget(self._table)

        buttons = QHBoxLayout()
        self._add_button = QPushButton("+ Field")
        self._remove_button = QPushButton("- Field")
        buttons.addWidget(self._add_button)
        buttons.addWidget(self._remove_button)
        buttons.addStretch(1)
        fields_layout.addLayout(buttons)

        left_layout.addWidget(fields_box)
        left_layout.addStretch(1)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(6, 6, 6, 6)

        plot_box = QGroupBox("Magnetic profile")
        plot_layout = QVBoxLayout(plot_box)

        controls = QFormLayout()
        components_row = QHBoxLayout()
        self._bx = QCheckBox("Bx")
        self._by = QCheckBox("By")
        self._bz = QCheckBox("Bz")
        self._bx.setChecked(True)
        self._by.setChecked(True)
        self._bz.setChecked(True)
        components_row.addWidget(self._bx)
        components_row.addWidget(self._by)
        components_row.addWidget(self._bz)
        components_row.addStretch(1)

        components_widget = QWidget()
        components_widget.setLayout(components_row)
        self._start = QLineEdit("0, 0, 0")
        self._direction = QLineEdit("0, 0, 1")
        self._points = QSpinBox()
        self._points.setRange(2, 10000)
        self._points.setValue(300)
        self._step = _double_spin(1.0e-9, 10.0, 8, 0.001)
        self._step.setValue(0.001)

        controls.addRow("Components", components_widget)
        controls.addRow("Start point x,y,z", self._start)
        controls.addRow("Direction dirx,diry,dirz", self._direction)
        controls.addRow(
            "Number of points",
            self._window._build_parameterized_editor(
                lambda: "magneticField.plotProfile.points",
                "Number of points",
                self._points,
            ),
        )
        controls.addRow(
            "Step [m]",
            self._window._build_parameterized_editor(
                lambda: "magneticField.plotProfile.step",
                "Step [m]",
                self._step,
            ),
        )
        plot_layout.addLayout(controls)

        self._plot_error = QLabel("")
        self._plot_error.setStyleSheet("color: #8f1d1d;")
        self._plot_error.setWordWrap(True)
        plot_layout.addWidget(self._plot_error)

        self._canvas = MagneticProfileCanvas()
        self._toolbar = NavigationToolbar2QT(self._canvas, plot_box)
        configure_matplotlib_toolbar(self._toolbar)
        plot_layout.addWidget(self._toolbar)
        plot_layout.addWidget(self._canvas, 1)

        right_layout.addWidget(plot_box)

        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 2)
        splitter.setSizes([560, 920])

        self._add_button.clicked.connect(self._add_row)
        self._remove_button.clicked.connect(self._remove_row)

        for widget in (
            self._enabled,
            self._start,
            self._direction,
            self._points,
            self._step,
            self._bx,
            self._by,
            self._bz,
        ):
            if isinstance(widget, QLineEdit):
                widget.textChanged.connect(self._on_changed)
            elif isinstance(widget, (QSpinBox, QDoubleSpinBox)):
                widget.valueChanged.connect(self._on_changed)
            elif isinstance(widget, QCheckBox):
                widget.toggled.connect(self._on_changed)

    def _on_changed(self) -> None:
        self._window.schedule_preview_refresh()
        self.refresh_plot()

    def _field_name_widget(self, value: str = "") -> QLineEdit:
        widget = QLineEdit(value)
        widget.textChanged.connect(self._on_changed)
        return widget

    def _field_type_widget(self, value: str = "constant") -> QComboBox:
        widget = QComboBox()
        widget.addItems(["constant", "file"])
        normalized = value.strip().lower()
        if normalized not in {"constant", "file"}:
            normalized = "constant"
        widget.setCurrentText(normalized)
        return widget

    def _field_scale_widget(self, value: float = 1.0) -> ParameterizedEditor:
        widget = _double_spin(0.0, 1.0e9, 6, 0.05)
        widget.setValue(value)
        widget.valueChanged.connect(self._on_changed)
        return self._window._build_parameterized_editor(
            lambda w=widget: f"magneticField.fields[{self._row_for_widget(w, 2)}].scale",
            "Field scale",
            widget,
        )

    def _field_payload_widget(self, source_type: str, value: str = "") -> _FieldPayloadWidget:
        return _FieldPayloadWidget(
            self._window,
            source_type=source_type,
            value=value,
            changed_callback=self._on_changed,
        )

    def _on_row_source_type_changed(self, payload_widget: _FieldPayloadWidget, source_type: str) -> None:
        payload_widget.set_source_type(source_type)
        self._on_changed()

    def _add_row(self, row_data: dict[str, Any] | None = None) -> None:
        row = self._table.rowCount()
        self._table.insertRow(row)
        row_data = row_data or {}

        name = str(row_data.get("name", f"field_{row + 1}"))
        source_type = str(row_data.get("sourceType", "constant"))
        scale = float(row_data.get("scale", 1.0))

        payload = ""
        if source_type == "constant":
            constant = row_data.get("constantValue")
            if isinstance(constant, dict):
                payload = f"{float(constant.get('bx', 0.0)):g}, {float(constant.get('by', 0.0)):g}, {float(constant.get('bz', 0.0)):g}"
        else:
            payload = str(row_data.get("filePath", ""))

        self._table.setCellWidget(row, 0, self._field_name_widget(name))
        type_widget = self._field_type_widget(source_type)
        payload_widget = self._field_payload_widget(source_type, payload)
        type_widget.currentTextChanged.connect(
            lambda selected, widget=payload_widget: self._on_row_source_type_changed(widget, selected)
        )
        self._table.setCellWidget(row, 1, type_widget)
        self._table.setCellWidget(row, 2, self._field_scale_widget(scale))
        self._table.setCellWidget(row, 3, payload_widget)
        self._on_changed()

    def _remove_row(self) -> None:
        row = self._table.currentRow()
        if row < 0:
            row = self._table.rowCount() - 1
        if row >= 0:
            self._table.removeRow(row)
            self._on_changed()

    def _iter_fields(self) -> list[dict[str, Any]]:
        fields: list[dict[str, Any]] = []
        for row in range(self._table.rowCount()):
            name_widget = self._table.cellWidget(row, 0)
            type_widget = self._table.cellWidget(row, 1)
            scale_widget = self._table.cellWidget(row, 2)
            payload_widget = self._table.cellWidget(row, 3)
            if not isinstance(name_widget, QLineEdit):
                continue
            if not isinstance(type_widget, QComboBox):
                continue
            scale_widget = self._unwrap_editor(scale_widget, QDoubleSpinBox)
            if not isinstance(scale_widget, QDoubleSpinBox):
                continue
            if not isinstance(payload_widget, _FieldPayloadWidget):
                continue

            name = name_widget.text().strip()
            source_type = type_widget.currentText().strip().lower()
            scale = float(scale_widget.value())
            payload = payload_widget.text().strip()
            if not name:
                continue
            if source_type not in {"constant", "file"}:
                raise WorkflowError(f"magneticField.fields[{row}].sourceType must be 'constant' or 'file'")

            entry: dict[str, Any] = {
                "name": name,
                "sourceType": source_type,
                "scale": scale,
            }
            if source_type == "constant":
                values = parse_number_list(payload, f"magneticField.fields[{row}].constantValue")
                if len(values) != 3:
                    raise WorkflowError(f"magneticField.fields[{row}].constantValue must contain exactly 3 values")
                entry["constantValue"] = {
                    "bx": float(values[0]),
                    "by": float(values[1]),
                    "bz": float(values[2]),
                }
            else:
                if not payload:
                    raise WorkflowError(f"magneticField.fields[{row}].filePath is required")
                entry["filePath"] = payload
            fields.append(entry)
        return fields

    def collect(self) -> dict[str, Any]:
        profile_components: list[str] = []
        if self._bx.isChecked():
            profile_components.append("Bx")
        if self._by.isChecked():
            profile_components.append("By")
        if self._bz.isChecked():
            profile_components.append("Bz")

        magnetic: dict[str, Any] = {
            "enabled": self._enabled.isChecked(),
            "fields": self._iter_fields(),
            "plotProfile": {
                "components": profile_components,
                "startPoint": self._start.text().strip(),
                "direction": self._direction.text().strip(),
                "points": int(self._points.value()),
                "step": float(self._step.value()),
            },
        }
        return magnetic

    def populate(self, spec: dict[str, Any]) -> None:
        magnetic = nested_get(spec, "magneticField", default={})
        if not isinstance(magnetic, dict):
            magnetic = {}

        self._enabled.setChecked(bool(magnetic.get("enabled", False)))

        self._table.setRowCount(0)
        fields = magnetic.get("fields")
        directory = str(magnetic.get("directory", "")).strip()
        if isinstance(fields, list):
            for field in fields:
                if isinstance(field, dict):
                    normalized_field = dict(field)
                    if str(normalized_field.get("sourceType", "")).strip().lower() == "file":
                        normalized_field["filePath"] = _merge_directory_file_path(
                            directory,
                            str(normalized_field.get("filePath", "")),
                        )
                    self._add_row(normalized_field)

        # legacy fallback import
        if self._table.rowCount() == 0 and isinstance(magnetic, dict):
            source_mode = str(magnetic.get("sourceMode", "none"))
            legacy_scale = float(magnetic.get("scale", 1.0))
            if source_mode == "file" and magnetic.get("file"):
                self._add_row(
                    {
                        "name": "external",
                        "sourceType": "file",
                        "scale": legacy_scale,
                        "filePath": str(magnetic.get("file", "")),
                    }
                )
            elif source_mode == "directory" and magnetic.get("directory"):
                legacy_path = (Path(str(magnetic.get("directory", ""))) / "EXTfield.fld").as_posix()
                self._add_row(
                    {
                        "name": "external",
                        "sourceType": "file",
                        "scale": legacy_scale,
                        "filePath": legacy_path,
                    }
                )

        profile = magnetic.get("plotProfile") if isinstance(magnetic, dict) else None
        if isinstance(profile, dict):
            components = profile.get("components")
            selected = set(components) if isinstance(components, list) else {"Bx", "By", "Bz"}
            self._bx.setChecked("Bx" in selected)
            self._by.setChecked("By" in selected)
            self._bz.setChecked("Bz" in selected)
            self._start.setText(str(profile.get("startPoint", "0, 0, 0")))
            self._direction.setText(str(profile.get("direction", "0, 0, 1")))
            self._points.setValue(int(profile.get("points", 300)))
            self._step.setValue(float(profile.get("step", 0.001)))
        else:
            self._bx.setChecked(True)
            self._by.setChecked(True)
            self._bz.setChecked(True)
            self._start.setText("0, 0, 0")
            self._direction.setText("0, 0, 1")
            self._points.setValue(300)
            self._step.setValue(0.001)

        self.refresh_plot()

    def _resolve_field_path(self, file_path: str) -> Path:
        candidate = Path(file_path).expanduser()
        if candidate.is_absolute():
            return candidate

        if self._window.authoring_path is not None:
            return (self._window.authoring_path.parent / candidate).resolve()
        return (REPO_ROOT / candidate).resolve()

    def _load_field_grid(self, path: Path) -> dict[str, Any]:
        key = str(path)
        if key in self._field_cache:
            return self._field_cache[key]

        if not path.exists():
            raise WorkflowError(f"Magnetic field file not found: {path}")

        raw = np.loadtxt(path, comments="#")
        if raw.ndim == 1:
            raw = raw.reshape(1, -1)
        if raw.shape[1] < 6:
            raise WorkflowError(f"Magnetic field file must have at least 6 numeric columns: {path}")

        # Keep runtime parity with C++ loader (xscale = 1e-3).
        xyz = raw[:, :3] * 1.0e-3
        bvals = raw[:, 3:6]

        xs = np.unique(xyz[:, 0])
        ys = np.unique(xyz[:, 1])
        zs = np.unique(xyz[:, 2])

        values = np.zeros((len(xs), len(ys), len(zs), 3), dtype=float)
        index_x = {value: idx for idx, value in enumerate(xs)}
        index_y = {value: idx for idx, value in enumerate(ys)}
        index_z = {value: idx for idx, value in enumerate(zs)}
        for point, vector in zip(xyz, bvals):
            values[index_x[point[0]], index_y[point[1]], index_z[point[2]], :] = vector

        parsed = {
            "xs": xs,
            "ys": ys,
            "zs": zs,
            "values": values,
        }
        self._field_cache[key] = parsed
        return parsed

    @staticmethod
    def _axis_interp_indices(axis_values: np.ndarray, coordinate: float) -> tuple[int, int, float] | None:
        if axis_values.size == 0:
            return None
        if coordinate < axis_values[0] or coordinate > axis_values[-1]:
            return None
        if axis_values.size == 1:
            return 0, 0, 0.0

        upper = int(np.searchsorted(axis_values, coordinate, side="right"))
        if upper <= 0:
            return 0, 0, 0.0
        if upper >= axis_values.size:
            return axis_values.size - 1, axis_values.size - 1, 0.0
        lower = upper - 1
        span = float(axis_values[upper] - axis_values[lower])
        if span <= 0.0:
            return lower, upper, 0.0
        weight = float((coordinate - axis_values[lower]) / span)
        return lower, upper, max(0.0, min(1.0, weight))

    def _sample_grid(self, grid: dict[str, Any], point: np.ndarray) -> np.ndarray:
        ix = self._axis_interp_indices(grid["xs"], float(point[0]))
        iy = self._axis_interp_indices(grid["ys"], float(point[1]))
        iz = self._axis_interp_indices(grid["zs"], float(point[2]))
        if ix is None or iy is None or iz is None:
            return np.zeros(3, dtype=float)

        i0, i1, wx = ix
        j0, j1, wy = iy
        k0, k1, wz = iz
        values = grid["values"]

        def lerp(a: np.ndarray, b: np.ndarray, w: float) -> np.ndarray:
            return a * (1.0 - w) + b * w

        c000 = values[i0, j0, k0]
        c100 = values[i1, j0, k0]
        c010 = values[i0, j1, k0]
        c110 = values[i1, j1, k0]
        c001 = values[i0, j0, k1]
        c101 = values[i1, j0, k1]
        c011 = values[i0, j1, k1]
        c111 = values[i1, j1, k1]

        c00 = lerp(c000, c100, wx)
        c10 = lerp(c010, c110, wx)
        c01 = lerp(c001, c101, wx)
        c11 = lerp(c011, c111, wx)
        c0 = lerp(c00, c10, wy)
        c1 = lerp(c01, c11, wy)
        return lerp(c0, c1, wz)

    def refresh_plot(self) -> None:
        self._plot_error.setText("")
        try:
            magnetic = self.collect()
            components = magnetic.get("plotProfile", {}).get("components", [])
            selected = [name for name in ("Bx", "By", "Bz") if name in components]
            if not selected:
                self._canvas.plot_placeholder("Select at least one component (Bx, By, or Bz).")
                return

            if not bool(magnetic.get("enabled", False)):
                self._canvas.plot_placeholder("Magnetic field is disabled.")
                return

            fields = magnetic.get("fields", [])
            if not fields:
                self._canvas.plot_placeholder("No magnetic field sources defined.")
                return

            start_values = parse_number_list(str(magnetic["plotProfile"].get("startPoint", "")), "magneticField.plotProfile.startPoint")
            direction_values = parse_number_list(str(magnetic["plotProfile"].get("direction", "")), "magneticField.plotProfile.direction")
            if len(start_values) != 3:
                raise WorkflowError("magneticField.plotProfile.startPoint must contain exactly 3 values")
            if len(direction_values) != 3:
                raise WorkflowError("magneticField.plotProfile.direction must contain exactly 3 values")

            start = np.asarray(start_values, dtype=float)
            direction = np.asarray(direction_values, dtype=float)
            norm = float(np.linalg.norm(direction))
            if norm <= 0.0:
                raise WorkflowError("magneticField.plotProfile.direction must be a non-zero vector")
            direction = direction / norm

            points = int(magnetic["plotProfile"].get("points", 300))
            if points < 2:
                raise WorkflowError("magneticField.plotProfile.points must be >= 2")
            step = float(magnetic["plotProfile"].get("step", 0.001))
            if step <= 0.0:
                raise WorkflowError("magneticField.plotProfile.step must be > 0")

            s_values = np.arange(points, dtype=float) * step
            sample_points = start + np.outer(s_values, direction)

            totals = {
                "Bx": np.zeros(points, dtype=float),
                "By": np.zeros(points, dtype=float),
                "Bz": np.zeros(points, dtype=float),
            }
            for field in fields:
                scale = float(field.get("scale", 1.0))
                if field.get("sourceType") == "constant":
                    constant = field.get("constantValue", {})
                    vector = np.asarray(
                        [
                            float(constant.get("bx", 0.0)),
                            float(constant.get("by", 0.0)),
                            float(constant.get("bz", 0.0)),
                        ],
                        dtype=float,
                    ) * scale
                    totals["Bx"] += vector[0]
                    totals["By"] += vector[1]
                    totals["Bz"] += vector[2]
                else:
                    resolved_path = self._resolve_field_path(str(field.get("filePath", "")))
                    grid = self._load_field_grid(resolved_path)
                    for idx, point in enumerate(sample_points):
                        sample = self._sample_grid(grid, point) * scale
                        totals["Bx"][idx] += sample[0]
                        totals["By"][idx] += sample[1]
                        totals["Bz"][idx] += sample[2]

            self._canvas.plot_components(s_values, totals, selected)
        except Exception as exc:
            self._plot_error.setText(str(exc))
            self._canvas.plot_placeholder("Plot unavailable due to invalid profile or field definition.")


def build_workspace(window) -> QWidget:
    workspace = _MagneticWorkspaceWidget(window)
    window.widgets["magnetic.workspace"] = workspace
    window.widgets["magnetic.enabled"] = workspace._enabled
    return workspace


def build_form(window) -> QFormLayout:
    # The Magnetic section uses a dedicated workspace with left controls and right plot.
    layout = QFormLayout()
    layout.addRow(QLabel("Magnetic section is rendered as a dedicated workspace."))
    return layout


def populate(window, spec: dict[str, object]) -> None:
    workspace = window.widgets.get("magnetic.workspace")
    if isinstance(workspace, _MagneticWorkspaceWidget):
        workspace.populate(spec)


def collect(window, spec: dict[str, object]) -> None:
    workspace = window.widgets.get("magnetic.workspace")
    if not isinstance(workspace, _MagneticWorkspaceWidget):
        return
    spec["magneticField"] = workspace.collect()
