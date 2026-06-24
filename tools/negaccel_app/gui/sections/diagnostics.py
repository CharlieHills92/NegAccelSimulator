"""Diagnostics configuration section for the NegAccel GUI."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QAbstractItemView, QHBoxLayout, QHeaderView, QLabel, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget
from ...workflow.domains.geometry_validation import derive_diagnostic_planes

from ..common import (
    QCheckBox,
    QFormLayout,
    QLineEdit,
    WorkflowError,
    format_number_list,
    nested_get,
    parse_number_list,
)


_DOMAIN_BOUNDARY_NAMES = {
    1: "x-min",
    2: "x-max",
    3: "y-min",
    4: "y-max",
    5: "z-min",
    6: "z-max",
}


class _GridPowerIdTable(QWidget):
    def __init__(self, window, parent=None):
        super().__init__(parent)
        self._window = window
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(4)

        intro = QLabel(
            "ID 0 means particles not ending on a boundary. IDs 1-6 are the domain faces; IDs above 6 are geometry solids. Range is derived from the current geometry."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

        self._table = QTableWidget(0, 4)
        self._table.setHorizontalHeaderLabels(["ID", "Name", "Derived Range [m]", "Include In Total"])
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(28)
        self._table.setShowGrid(False)
        header = self._table.horizontalHeader()
        header.setFixedHeight(26)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        root.addWidget(self._table)

    def _domain_z_bounds(self) -> tuple[float, float]:
        z_start_widget = self._window.widgets.get("geometry.domain.zStart")
        z_size_widget = self._window.widgets.get("geometry.domain.z")
        z_start = float(z_start_widget.value()) if hasattr(z_start_widget, "value") else 0.0
        z_size = float(z_size_widget.value()) if hasattr(z_size_widget, "value") else 0.0
        return z_start, z_start + z_size

    def _solid_map(self) -> dict[int, dict[str, object]]:
        solids_editor = self._window.widgets.get("geometry.solidsEditor")
        solids = solids_editor.get_solids() if hasattr(solids_editor, "get_solids") else []
        solid_map: dict[int, dict[str, object]] = {}
        for index, solid in enumerate(solids):
            if not isinstance(solid, dict):
                continue
            boundary_id = solid.get("boundaryId", 7 + index)
            if isinstance(boundary_id, int) and boundary_id >= 7:
                solid_map[boundary_id] = solid
        return solid_map

    def _derive_extent(self, row_id: int) -> str:
        z_start, z_end = self._domain_z_bounds()
        if row_id == 0:
            return f"inside domain, z: {z_start:g} .. {z_end:g}"
        if row_id in {1, 2, 3, 4}:
            return f"z: {z_start:g} .. {z_end:g}"
        if row_id == 5:
            return f"z = {z_start:g}"
        if row_id == 6:
            return f"z = {z_end:g}"

        solid = self._solid_map().get(row_id)
        z_profile = solid.get("zProfileMeters") if isinstance(solid, dict) else None
        if isinstance(z_profile, list) and z_profile:
            try:
                z_values = [float(value) for value in z_profile]
            except (TypeError, ValueError):
                return "derived at runtime"
            return f"z: {min(z_values):g} .. {max(z_values):g}"
        return "derived at runtime"

    def _inventory_rows(self) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = [
            {"id": 0, "name": "In volume", "extent": self._derive_extent(0)}
        ]

        boundary_editor = self._window.widgets.get("boundary.editor")
        boundary_rows = boundary_editor.get_boundaries() if hasattr(boundary_editor, "get_boundaries") else []
        if boundary_rows:
            for boundary in boundary_rows:
                boundary_id = boundary.get("boundaryId")
                if not isinstance(boundary_id, int) or boundary_id < 1:
                    continue
                rows.append(
                    {
                        "id": boundary_id,
                        "name": str(boundary.get("name") or f"Boundary {boundary_id}"),
                        "extent": self._derive_extent(boundary_id),
                    }
                )
            return rows

        solid_map = self._solid_map()
        for boundary_id in range(1, 7):
            rows.append(
                {
                    "id": boundary_id,
                    "name": _DOMAIN_BOUNDARY_NAMES[boundary_id],
                    "extent": self._derive_extent(boundary_id),
                }
            )
        for boundary_id in sorted(solid_map):
            solid = solid_map[boundary_id]
            rows.append(
                {
                    "id": boundary_id,
                    "name": str(solid.get("name") or f"Solid {boundary_id}"),
                    "extent": self._derive_extent(boundary_id),
                }
            )
        return rows

    def _current_rows_by_id(self) -> dict[int, dict[str, object]]:
        rows: dict[int, dict[str, object]] = {}
        for row in range(self._table.rowCount()):
            id_item = self._table.item(row, 0)
            include_widget = self._table.cellWidget(row, 3)
            if id_item is None or not isinstance(include_widget, QCheckBox):
                continue
            rows[int(id_item.text())] = {
                "id": int(id_item.text()),
                "includeInTotal": include_widget.isChecked(),
            }
        return rows

    def _default_include(self, row_id: int, row_name: str) -> bool:
        if row_id <= 6:
            return False
        return row_name.strip().upper() not in {"PG", "PLASMA GRID"}

    def _set_read_only_item(self, row: int, column: int, text: str) -> None:
        item = QTableWidgetItem(text)
        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
        self._table.setItem(row, column, item)

    def sync_from_geometry(self, ranges: list[dict[str, object]] | None = None) -> None:
        current = self._current_rows_by_id()
        if ranges is not None:
            current = {
                int(range_definition["id"]): dict(range_definition)
                for range_definition in ranges
                if isinstance(range_definition, dict) and isinstance(range_definition.get("id"), int)
            }

        self._table.setRowCount(0)
        for row, inventory_row in enumerate(self._inventory_rows()):
            row_id = int(inventory_row["id"])
            row_name = str(inventory_row["name"])
            include_in_total = bool(
                current.get(row_id, {}).get("includeInTotal", self._default_include(row_id, row_name))
            )

            self._table.insertRow(row)
            self._set_read_only_item(row, 0, str(row_id))
            self._set_read_only_item(row, 1, row_name)
            self._set_read_only_item(row, 2, str(inventory_row["extent"]))

            include_widget = QCheckBox()
            include_widget.setChecked(include_in_total)
            include_widget.stateChanged.connect(self._window.schedule_preview_refresh)
            self._table.setCellWidget(row, 3, include_widget)

    def get_ranges(self) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for row in range(self._table.rowCount()):
            id_item = self._table.item(row, 0)
            include_widget = self._table.cellWidget(row, 3)
            if id_item is None or not isinstance(include_widget, QCheckBox):
                continue
            rows.append(
                {
                    "id": int(id_item.text()),
                    "includeInTotal": include_widget.isChecked(),
                }
            )
        return rows


def _current_geometry_document(window) -> dict[str, object]:
    solids_editor = window.widgets.get("geometry.solidsEditor")
    if solids_editor is None or not hasattr(solids_editor, "get_solids"):
        return {
            "domain": {"zStartMeters": 0.0, "zSizeMeters": 0.0},
            "solids": [],
        }

    z_start_widget = window.widgets.get("geometry.domain.zStart")
    z_size_widget = window.widgets.get("geometry.domain.z")
    z_start = float(z_start_widget.value()) if hasattr(z_start_widget, "value") else 0.0
    z_size = float(z_size_widget.value()) if hasattr(z_size_widget, "value") else 0.0
    return {
        "domain": {
            "zStartMeters": z_start,
            "zSizeMeters": z_size,
        },
        "solids": solids_editor.get_solids(),
    }


def _set_diagnostic_plane_values(
    window,
    *,
    sample_positions: list[float] | None = None,
    summary_z: float | None = None,
    emitter_z: float | None = None,
    transmission_z: float | None = None,
) -> None:
    sample_widget = window.widgets.get("diagnostics.sampleZPositionsMeters")
    summary_widget = window.widgets.get("diagnostics.summaryZPositionMeters")
    emitter_widget = window.widgets.get("diagnostics.emitterExportZPositionMeters")
    transmission_widget = window.widgets.get("diagnostics.transmissionPlaneZPositionMeters")

    if isinstance(sample_widget, QLineEdit) and sample_positions is not None:
        sample_widget.blockSignals(True)
        sample_widget.setText(format_number_list(sample_positions))
        sample_widget.blockSignals(False)
    if hasattr(summary_widget, "setValue") and summary_z is not None:
        summary_widget.blockSignals(True)
        summary_widget.setValue(summary_z)
        summary_widget.blockSignals(False)
    if hasattr(emitter_widget, "setValue") and emitter_z is not None:
        emitter_widget.blockSignals(True)
        emitter_widget.setValue(emitter_z)
        emitter_widget.blockSignals(False)
    if hasattr(transmission_widget, "setValue") and transmission_z is not None:
        transmission_widget.blockSignals(True)
        transmission_widget.setValue(transmission_z)
        transmission_widget.blockSignals(False)


def _derived_planes_from_geometry(window) -> dict[str, object] | None:
    try:
        return derive_diagnostic_planes(_current_geometry_document(window))
    except WorkflowError:
        return None


def _sync_derived_planes_from_geometry(window) -> None:
    derived = _derived_planes_from_geometry(window)
    if derived is None:
        return

    _set_diagnostic_plane_values(
        window,
        sample_positions=list(derived["sampleZPositionsMeters"]),
        summary_z=float(derived["summaryZPositionMeters"]),
        emitter_z=float(derived["emitterExportZPositionMeters"]),
        transmission_z=float(derived["transmissionPlaneZPositionMeters"]),
    )


def _restore_derived_planes(window) -> None:
    _sync_derived_planes_from_geometry(window)
    window.schedule_preview_refresh()


def build_form(window) -> QFormLayout:
    layout = QFormLayout()

    window.widgets["diagnostics.sampleZPositionsMeters"] = QLineEdit()
    window.widgets["diagnostics.summaryZPositionMeters"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["diagnostics.emitterExportZPositionMeters"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["diagnostics.transmissionPlaneZPositionMeters"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["diagnostics.apertureRadiusMeters"] = window._double_spin(0.0, 10.0, 6, 0.001)

    window.widgets["diagnostics.sampleZPositionsMeters"].setToolTip(
        "Defaults to the geometry-derived solid boundary z extents. You can edit the list manually or restore the geometry defaults."
    )
    window.widgets["diagnostics.summaryZPositionMeters"].setToolTip(
        "Defaults to the maximum z coordinate of the simulation domain. You can edit it manually or restore the geometry defaults."
    )
    window.widgets["diagnostics.emitterExportZPositionMeters"].setToolTip(
        "Defaults to the maximum z coordinate of the simulation domain. You can edit it manually or restore the geometry defaults."
    )

    window.widgets["diagnostics.transmissionPlaneZPositionMeters"].setToolTip(
        "Defaults to the maximum z coordinate of solid boundary ID 7 and falls back to domain max z when ID 7 is absent. You can edit it manually or restore the geometry defaults."
    )

    window.widgets["diagnostics.writePerSpeciesDiagnostics"] = QCheckBox("Write per-species diagnostic files")
    window.widgets["diagnostics.writePerSpeciesGridPower"] = QCheckBox("Write per-species grid power summaries")
    window.widgets["diagnostics.writePerSpeciesPlots"] = QCheckBox("Render per-species plots")
    window.widgets["diagnostics.writeNegativeIonSummary"] = QCheckBox("Write dedicated negative-ion summary")

    restore_plane_defaults_button = QPushButton("Restore to defaults")
    restore_plane_defaults_button.clicked.connect(lambda: _restore_derived_planes(window))
    restore_plane_defaults_row = QHBoxLayout()
    restore_plane_defaults_row.addStretch()
    restore_plane_defaults_row.addWidget(restore_plane_defaults_button)
    restore_plane_defaults_widget = QWidget()
    restore_plane_defaults_widget.setLayout(restore_plane_defaults_row)

    grid_ranges = _GridPowerIdTable(window)
    window.widgets["diagnostics.gridRanges"] = grid_ranges
    window.sync_diagnostics_grid_power_from_geometry = grid_ranges.sync_from_geometry
    grid_ranges.sync_from_geometry()
    _sync_derived_planes_from_geometry(window)

    layout.addRow("Along-Z sample planes [m]", window.widgets["diagnostics.sampleZPositionsMeters"])
    layout.addRow("Summary plane z [m]", window.widgets["diagnostics.summaryZPositionMeters"])
    layout.addRow("Emitter export plane z [m]", window.widgets["diagnostics.emitterExportZPositionMeters"])
    layout.addRow("Extracted-current evaluation plane z [m]", window.widgets["diagnostics.transmissionPlaneZPositionMeters"])
    layout.addRow("", restore_plane_defaults_widget)
    layout.addRow("Aperture radius [m]", window.widgets["diagnostics.apertureRadiusMeters"])
    layout.addRow(window.widgets["diagnostics.writePerSpeciesDiagnostics"])
    layout.addRow(window.widgets["diagnostics.writePerSpeciesGridPower"])
    layout.addRow(window.widgets["diagnostics.writePerSpeciesPlots"])
    layout.addRow(window.widgets["diagnostics.writeNegativeIonSummary"])
    layout.addRow("Grid-power ranges", window.widgets["diagnostics.gridRanges"])
    return layout


def populate(window, spec: dict[str, object]) -> None:
    window.widgets["diagnostics.apertureRadiusMeters"].setValue(
        float(nested_get(spec, "diagnostics", "summary", "apertureRadiusMeters", default=0.007))
    )

    window.widgets["diagnostics.writePerSpeciesDiagnostics"].setChecked(
        bool(nested_get(spec, "diagnostics", "species", "writePerSpeciesDiagnostics", default=True))
    )
    window.widgets["diagnostics.writePerSpeciesGridPower"].setChecked(
        bool(nested_get(spec, "diagnostics", "species", "writePerSpeciesGridPower", default=True))
    )
    window.widgets["diagnostics.writePerSpeciesPlots"].setChecked(
        bool(nested_get(spec, "diagnostics", "species", "writePerSpeciesPlots", default=True))
    )
    window.widgets["diagnostics.writeNegativeIonSummary"].setChecked(
        bool(nested_get(spec, "diagnostics", "species", "writeNegativeIonSummary", default=True))
    )

    window.widgets["diagnostics.gridRanges"].sync_from_geometry(
        list(nested_get(spec, "diagnostics", "gridPower", "ranges", default=[]))
    )
    derived = _derived_planes_from_geometry(window)
    sample_positions = nested_get(spec, "diagnostics", "planes", "sampleZPositionsMeters", default=None)
    if sample_positions is None and derived is not None:
        sample_positions = list(derived["sampleZPositionsMeters"])
    summary_z = nested_get(spec, "diagnostics", "planes", "summaryZPositionMeters", default=None)
    if summary_z is None and derived is not None:
        summary_z = float(derived["summaryZPositionMeters"])
    emitter_z = nested_get(spec, "diagnostics", "planes", "emitterExportZPositionMeters", default=None)
    if emitter_z is None and derived is not None:
        emitter_z = float(derived["emitterExportZPositionMeters"])
    transmission_z = nested_get(
        spec,
        "diagnostics",
        "summary",
        "transmissionPlaneZPositionMeters",
        default=None,
    )
    if transmission_z is None and derived is not None:
        transmission_z = float(derived["transmissionPlaneZPositionMeters"])

    _set_diagnostic_plane_values(
        window,
        sample_positions=list(sample_positions) if sample_positions is not None else None,
        summary_z=float(summary_z) if summary_z is not None else None,
        emitter_z=float(emitter_z) if emitter_z is not None else None,
        transmission_z=float(transmission_z) if transmission_z is not None else None,
    )


def collect(window, spec: dict[str, object]) -> None:
    diagnostics = spec.setdefault("diagnostics", {})

    planes = diagnostics.setdefault("planes", {})
    sample_positions = parse_number_list(
        window.widgets["diagnostics.sampleZPositionsMeters"].text(),
        "diagnostics.planes.sampleZPositionsMeters",
    )
    if not sample_positions:
        raise WorkflowError("diagnostics.planes.sampleZPositionsMeters must contain at least one value")
    planes["sampleZPositionsMeters"] = sample_positions
    planes["summaryZPositionMeters"] = float(window.widgets["diagnostics.summaryZPositionMeters"].value())
    planes["emitterExportZPositionMeters"] = float(window.widgets["diagnostics.emitterExportZPositionMeters"].value())

    species = diagnostics.setdefault("species", {})
    species["writePerSpeciesDiagnostics"] = window.widgets[
        "diagnostics.writePerSpeciesDiagnostics"
    ].isChecked()
    species["writePerSpeciesGridPower"] = window.widgets[
        "diagnostics.writePerSpeciesGridPower"
    ].isChecked()
    species["writePerSpeciesPlots"] = window.widgets["diagnostics.writePerSpeciesPlots"].isChecked()
    species["writeNegativeIonSummary"] = window.widgets[
        "diagnostics.writeNegativeIonSummary"
    ].isChecked()

    grid_power = diagnostics.setdefault("gridPower", {})
    grid_power["ranges"] = window.widgets["diagnostics.gridRanges"].get_ranges()
    if not grid_power["ranges"]:
        raise WorkflowError("diagnostics.gridPower.ranges must contain at least one row")

    summary = diagnostics.setdefault("summary", {})
    summary["apertureRadiusMeters"] = float(window.widgets["diagnostics.apertureRadiusMeters"].value())
    summary["transmissionPlaneZPositionMeters"] = float(
        window.widgets["diagnostics.transmissionPlaneZPositionMeters"].value()
    )