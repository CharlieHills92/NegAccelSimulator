"""Boundary-condition form section for the NegAccel GUI."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHeaderView,
    QLabel,
    QLineEdit,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..common import ParameterizedEditor, QCheckBox, WorkflowError, nested_get
from ...workflow.domains.boundary import resolve_boundary_value_expressions


_DOMAIN_BOUNDARIES = [
    {"boundaryId": 1, "name": "x-min", "conditionType": "neumann", "value": 0.0},
    {"boundaryId": 2, "name": "x-max", "conditionType": "neumann", "value": 0.0},
    {"boundaryId": 3, "name": "y-min", "conditionType": "neumann", "value": 0.0},
    {"boundaryId": 4, "name": "y-max", "conditionType": "neumann", "value": 0.0},
    {"boundaryId": 5, "name": "z-min", "conditionType": "dirichlet", "value": 0.0},
    {"boundaryId": 6, "name": "z-max", "conditionType": "neumann", "value": 0.0},
]

_LEGACY_STAGE_NAMES = [
    {"GG", "G1", "AG1"},
    {"REP", "G2", "AG2"},
    {"G3", "AG3"},
    {"G4", "AG4"},
    {"G5", "AG5"},
]


class _BoundaryEditorWidget(QWidget):
    def __init__(self, window, parent=None):
        super().__init__(parent)
        self._window = window
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(4)
        self._table = QTableWidget(0, 5)
        self._table.setHorizontalHeaderLabels(["ID", "Name", "Type", "Expression", "Value"])
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(28)
        self._table.setShowGrid(False)
        self._table.setStyleSheet(
            "QComboBox, QLineEdit { padding: 2px 4px; min-height: 22px; }"
        )
        header = self._table.horizontalHeader()
        header.setFixedHeight(26)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(4, QHeaderView.ResizeMode.Stretch)
        root.addWidget(self._table)

    def _format_boundary_value(self, value: object) -> str:
        if isinstance(value, bool):
            return str(value)
        if isinstance(value, (int, float)):
            return format(float(value), ".12g")
        return str(value)

    def _parse_boundary_value(self, raw_value: str, context: str) -> object:
        value = raw_value.strip()
        if not value:
            raise WorkflowError(f"{context} value must not be empty")
        try:
            return float(value)
        except ValueError:
            return value

    def _expression_widget(self, row: int) -> QLineEdit | None:
        widget = self._table.cellWidget(row, 3)
        if isinstance(widget, QLineEdit):
            return widget
        if isinstance(widget, ParameterizedEditor):
            editor = widget.editor_widget()
            if isinstance(editor, QLineEdit):
                return editor
        return None

    def _set_resolved_value(self, row: int, text: str, tooltip: str | None = None) -> None:
        item = self._table.item(row, 4)
        if item is None:
            item = QTableWidgetItem()
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self._table.setItem(row, 4, item)
        item.setText(text)
        item.setToolTip(tooltip or text)

    def _refresh_value_hints(self) -> None:
        boundaries: list[dict[str, object]] = []
        row_by_id: dict[int, int] = {}
        for row in range(self._table.rowCount()):
            boundary_id_item = self._table.item(row, 0)
            name_item = self._table.item(row, 1)
            condition_widget = self._table.cellWidget(row, 2)
            expression_widget = self._expression_widget(row)
            if boundary_id_item is None or name_item is None:
                continue
            if not isinstance(condition_widget, QComboBox) or not isinstance(expression_widget, QLineEdit):
                continue
            boundary_id = int(boundary_id_item.text())
            boundaries.append(
                {
                    "boundaryId": boundary_id,
                    "name": name_item.text(),
                    "conditionType": condition_widget.currentText(),
                    "value": self._parse_boundary_value(
                        expression_widget.text(),
                        f"Boundary {boundary_id} ({name_item.text()})",
                    ),
                }
            )
            row_by_id[boundary_id] = row

        try:
            resolved = resolve_boundary_value_expressions(boundaries, context_prefix="boundary editor")
        except WorkflowError as exc:
            for row in row_by_id.values():
                self._set_resolved_value(row, "Invalid", str(exc))
            return

        for boundary in resolved:
            row = row_by_id.get(int(boundary["boundaryId"]))
            if row is not None:
                self._set_resolved_value(row, self._format_boundary_value(boundary["value"]))

    def inventory_rows(self) -> list[dict[str, object]]:
        rows = [dict(boundary) for boundary in _DOMAIN_BOUNDARIES]
        solids_editor = self._window.widgets.get("geometry.solidsEditor")

        solids = solids_editor.get_solids() if hasattr(solids_editor, "get_solids") else []
        geometry_path = getattr(solids_editor, "_geometry_path", None)

        if geometry_path is not None or solids:
            used_ids: set[int] = set()
            for index, solid in enumerate(solids):
                boundary_id = solid.get("boundaryId", 7 + index)
                if not isinstance(boundary_id, int) or boundary_id < 7 or boundary_id in used_ids:
                    continue
                used_ids.add(boundary_id)
                rows.append(
                    {
                        "boundaryId": boundary_id,
                        "name": str(solid.get("name") or f"Solid {boundary_id}"),
                        "conditionType": "dirichlet",
                        "value": 0.0,
                    }
                )
            return rows
        return rows

    def _current_boundaries_by_id(self) -> dict[int, dict[str, object]]:
        boundaries: dict[int, dict[str, object]] = {}
        for row in range(self._table.rowCount()):
            boundary_id_item = self._table.item(row, 0)
            name_item = self._table.item(row, 1)
            condition_widget = self._table.cellWidget(row, 2)
            value_widget = self._expression_widget(row)
            if boundary_id_item is None or name_item is None:
                continue
            if not isinstance(condition_widget, QComboBox) or not isinstance(value_widget, QLineEdit):
                continue
            context = f"Boundary {boundary_id_item.text()} ({name_item.text()})"
            boundaries[int(boundary_id_item.text())] = {
                "boundaryId": int(boundary_id_item.text()),
                "name": name_item.text(),
                "conditionType": condition_widget.currentText(),
                "value": self._parse_boundary_value(value_widget.text(), context),
            }
        return boundaries

    def sync_from_geometry(self, boundaries: list[dict[str, object]] | None = None) -> None:
        current = self._current_boundaries_by_id()
        if boundaries is not None:
            current = {
                boundary["boundaryId"]: dict(boundary)
                for boundary in boundaries
                if isinstance(boundary, dict) and isinstance(boundary.get("boundaryId"), int)
            }

        merged_rows: list[dict[str, object]] = []
        for boundary in self.inventory_rows():
            current_row = current.get(boundary["boundaryId"], {})
            merged_rows.append(
                {
                    "boundaryId": boundary["boundaryId"],
                    "name": boundary["name"],
                    "conditionType": str(current_row.get("conditionType", boundary["conditionType"])).lower(),
                    "value": current_row.get("value", boundary["value"]),
                }
            )

        self._table.setRowCount(0)
        for row, boundary in enumerate(merged_rows):
            self._table.insertRow(row)

            boundary_id_item = QTableWidgetItem(str(boundary["boundaryId"]))
            boundary_id_item.setFlags(boundary_id_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self._table.setItem(row, 0, boundary_id_item)

            name_item = QTableWidgetItem(str(boundary["name"]))
            name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self._table.setItem(row, 1, name_item)

            condition_widget = QComboBox()
            condition_widget.addItems(["dirichlet", "neumann"])
            condition_widget.setCurrentText(str(boundary["conditionType"]).lower())
            condition_widget.currentTextChanged.connect(self._window.schedule_preview_refresh)
            condition_widget.currentTextChanged.connect(self._refresh_value_hints)
            self._table.setCellWidget(row, 2, condition_widget)

            value_widget = QLineEdit(self._format_boundary_value(boundary["value"]))
            value_widget.setPlaceholderText("e.g. 7000 or 5 * EG")
            value_widget.textChanged.connect(self._window.schedule_preview_refresh)
            value_widget.textChanged.connect(self._refresh_value_hints)
            parameterized_value_widget = self._window._build_parameterized_editor(
                lambda row_index=row: f"boundaryConditions.boundaries[{row_index}].value",
                f"Boundary {boundary['boundaryId']} ({boundary['name']}) expression",
                value_widget,
                parameter_type="string",
                value_reader=lambda widget=value_widget: widget.text().strip(),
            )
            self._table.setCellWidget(row, 3, parameterized_value_widget)

            self._set_resolved_value(row, self._format_boundary_value(boundary["value"]))

        self._refresh_value_hints()

    def get_boundaries(self) -> list[dict[str, object]]:
        boundaries: list[dict[str, object]] = []
        current = self._current_boundaries_by_id()
        for boundary_id in sorted(current):
            boundaries.append(current[boundary_id])
        return boundaries


def _legacy_boundaries(editor: _BoundaryEditorWidget, spec: dict[str, object]) -> list[dict[str, object]]:
    boundaries = {boundary["boundaryId"]: dict(boundary) for boundary in editor.inventory_rows()}
    extraction_grid = nested_get(spec, "boundaryConditions", "gridVoltagesVolts", "extractionGrid", default=None)
    if isinstance(extraction_grid, (int, float)):
        for boundary in boundaries.values():
            if str(boundary.get("name", "")).upper() == "EG":
                boundary["value"] = float(extraction_grid)
                break

    accelerator_stages = nested_get(
        spec,
        "boundaryConditions",
        "gridVoltagesVolts",
        "acceleratorStages",
        default=[],
    )
    if isinstance(accelerator_stages, list):
        for index, voltage in enumerate(accelerator_stages[: len(_LEGACY_STAGE_NAMES)]):
            if not isinstance(voltage, (int, float)):
                continue
            for boundary in boundaries.values():
                if str(boundary.get("name", "")).upper() in _LEGACY_STAGE_NAMES[index]:
                    boundary["value"] = float(voltage)
                    break

    return [boundaries[boundary_id] for boundary_id in sorted(boundaries)]


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    container = QWidget()
    root = QVBoxLayout(container)
    root.setContentsMargins(0, 0, 0, 0)
    root.setSpacing(10)

    intro = QLabel(
        "Configure simulation boundary IDs directly. IDs 1-6 are the domain faces; IDs above 6 are solid objects. Values may be numbers or expressions such as 5 * 700 or 5 * EG."
    )
    intro.setWordWrap(True)
    root.addWidget(intro)

    editor = _BoundaryEditorWidget(window)
    window.widgets["boundary.editor"] = editor
    window.sync_boundary_editor_from_geometry = editor.sync_from_geometry
    editor.sync_from_geometry()
    root.addWidget(editor)

    periodic_box = QGroupBox("Periodic boundaries")
    periodic_layout = QFormLayout(periodic_box)
    window.widgets["boundary.periodic.enabled"] = QCheckBox("Enable periodic boundaries")
    window.widgets["boundary.periodic.xMin"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["boundary.periodic.xMax"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["boundary.periodic.yMin"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["boundary.periodic.yMax"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    periodic_layout.addRow(window.widgets["boundary.periodic.enabled"])
    periodic_layout.addRow("Periodic x min [m]", window.widgets["boundary.periodic.xMin"])
    periodic_layout.addRow("Periodic x max [m]", window.widgets["boundary.periodic.xMax"])
    periodic_layout.addRow("Periodic y min [m]", window.widgets["boundary.periodic.yMin"])
    periodic_layout.addRow("Periodic y max [m]", window.widgets["boundary.periodic.yMax"])
    root.addWidget(periodic_box)

    layout.addRow(container)
    return layout


def populate(window, spec: dict[str, object]) -> None:
    boundaries = nested_get(spec, "boundaryConditions", "boundaries", default=None)
    editor = window.widgets["boundary.editor"]
    if isinstance(boundaries, list) and boundaries:
        editor.sync_from_geometry(boundaries)
    else:
        editor.sync_from_geometry(_legacy_boundaries(editor, spec))

    window.widgets["boundary.periodic.enabled"].setChecked(
        bool(nested_get(spec, "boundaryConditions", "periodicBoundaries", "enabled", default=False))
    )
    window.widgets["boundary.periodic.xMin"].setValue(
        float(nested_get(spec, "boundaryConditions", "periodicBoundaries", "xMinMeters", default=0.0))
    )
    window.widgets["boundary.periodic.xMax"].setValue(
        float(nested_get(spec, "boundaryConditions", "periodicBoundaries", "xMaxMeters", default=0.0))
    )
    window.widgets["boundary.periodic.yMin"].setValue(
        float(nested_get(spec, "boundaryConditions", "periodicBoundaries", "yMinMeters", default=0.0))
    )
    window.widgets["boundary.periodic.yMax"].setValue(
        float(nested_get(spec, "boundaryConditions", "periodicBoundaries", "yMaxMeters", default=0.0))
    )


def collect(window, spec: dict[str, object]) -> None:
    boundary = spec.setdefault("boundaryConditions", {})
    boundary["boundaries"] = window.widgets["boundary.editor"].get_boundaries()
    if "gridVoltagesVolts" in boundary:
        boundary.pop("gridVoltagesVolts")

    seen_ids: set[int] = set()
    for boundary_definition in boundary["boundaries"]:
        boundary_id = boundary_definition["boundaryId"]
        if boundary_id in seen_ids:
            raise WorkflowError(f"Duplicate boundaryConditions.boundaries boundaryId: {boundary_id}")
        seen_ids.add(boundary_id)

    if window.widgets["boundary.periodic.enabled"].isChecked():
        boundary["periodicBoundaries"] = {
            "enabled": True,
            "xMinMeters": float(window.widgets["boundary.periodic.xMin"].value()),
            "xMaxMeters": float(window.widgets["boundary.periodic.xMax"].value()),
            "yMinMeters": float(window.widgets["boundary.periodic.yMin"].value()),
            "yMaxMeters": float(window.widgets["boundary.periodic.yMax"].value()),
        }
    elif "periodicBoundaries" in boundary:
        boundary.pop("periodicBoundaries")
