"""Geometry form section for the NegAccel GUI – structured solid editor."""

from __future__ import annotations

import copy
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QMenu,
    QFileDialog,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from ..common import WorkflowError, nested_get, serialize_geometry_file_path
from ..visualization import GeometryCanvas, NavigationToolbar2QT, configure_matplotlib_toolbar
from ...workflow import TEMPLATES_DIR, load_geometry_file, load_geometry_reference, write_geometry_file


# ── widget helpers ────────────────────────────────────────────────────────────

def _dspin(min_: float, max_: float, decimals: int, step: float) -> QDoubleSpinBox:
    w = QDoubleSpinBox()
    w.setRange(min_, max_)
    w.setDecimals(decimals)
    w.setSingleStep(step)
    return w


def _ispin(min_: int, max_: int, value: int = 0) -> QSpinBox:
    w = QSpinBox()
    w.setRange(min_, max_)
    w.setValue(value)
    return w


def _default_profile() -> dict:
    return {
        "zProfileMeters": [0.0, 0.001],
        "rProfileMeters": [0.01, 0.01],
        "roundingRadiiMeters": [0.0, 0.0],
    }


def _mtf_aperture_pattern(
    layout: str = "rectangular-grid",
    *,
    x_offset: float = 0.0,
    y_offset: float = 0.0,
    row_shift: float = 0.0,
) -> dict:
    pattern = {
        "layout": layout,
        "countX": 5,
        "countY": 5,
        "pitchXMeters": 0.019,
        "pitchYMeters": 0.021,
        "outsidePatternIsSolid": True,
    }
    if abs(x_offset) > 1.0e-12:
        pattern["xOffsetMeters"] = x_offset
    if abs(y_offset) > 1.0e-12:
        pattern["yOffsetMeters"] = y_offset
    if layout == "staggered-grid" and abs(row_shift) > 1.0e-12:
        pattern["rowShiftXMeters"] = row_shift
    return pattern


def _mtf_solid(
    name: str,
    boundary_id: int,
    z_profile: list[float],
    r_profile: list[float],
    aperture_pattern: dict,
) -> dict:
    return {
        "name": name,
        "kind": "solid",
        "boundaryId": boundary_id,
        "zProfileMeters": z_profile,
        "rProfileMeters": r_profile,
        "roundingRadiiMeters": [0.0] * len(z_profile),
        "aperturePattern": copy.deepcopy(aperture_pattern),
    }


def _normalize_solid(solid: dict) -> dict:
    normalized = copy.deepcopy(solid)
    normalized_kind = str(normalized.get("kind", "solid") or "solid").strip()
    normalized["kind"] = "diagnosticPlane" if normalized_kind == "diagnosticPlane" else "solid"
    if "zProfileMeters" not in normalized or "rProfileMeters" not in normalized:
        normalized.update(_default_profile())
        return normalized

    z_values = normalized.get("zProfileMeters") or []
    if "roundingRadiiMeters" not in normalized or len(normalized.get("roundingRadiiMeters") or []) != len(z_values):
        normalized["roundingRadiiMeters"] = [0.0] * len(z_values)
    return normalized


def _load_geometry_document(path: Path) -> dict:
    return load_geometry_file(path)


def _validate_solids(solids: list[dict]) -> None:
    seen_names: set[str] = set()
    for solid_index, solid in enumerate(solids):
        solid_name = str(solid.get("name", "")).strip()
        label = solid_name or f"Solid {solid_index + 1}"
        if not solid_name:
            raise WorkflowError(f"{label}: name is required")
        normalized_name = solid_name.upper()
        if normalized_name in seen_names:
            raise WorkflowError(f"Duplicate solid name: {solid_name}")
        seen_names.add(normalized_name)

        boundary_id = solid.get("boundaryId", 7 + solid_index)
        if not isinstance(boundary_id, int) or boundary_id < 7:
            raise WorkflowError(f"{solid_name}: boundaryId must be an integer >= 7")

        z_values = solid.get("zProfileMeters")
        r_values = solid.get("rProfileMeters")
        rounding_values = solid.get("roundingRadiiMeters")
        if not isinstance(z_values, list) or not isinstance(r_values, list):
            raise WorkflowError(f"{solid_name}: profile points are required")
        if len(z_values) < 2 or len(r_values) < 2:
            raise WorkflowError(f"{solid_name}: add at least two profile points")
        if len(z_values) != len(r_values):
            raise WorkflowError(f"{solid_name}: z and r point counts must match")
        if rounding_values is not None:
            if not isinstance(rounding_values, list) or len(rounding_values) != len(z_values):
                raise WorkflowError(f"{solid_name}: rounding list must match the number of profile points")
        else:
            rounding_values = [0.0] * len(z_values)

        previous_z = None
        has_positive_span = False
        for point_index, (z_value, r_value, rounding_value) in enumerate(
            zip(z_values, r_values, rounding_values), start=1
        ):
            try:
                z_number = float(z_value)
                r_number = float(r_value)
                rounding_number = float(rounding_value)
            except (TypeError, ValueError) as exc:
                raise WorkflowError(f"{solid_name}: point {point_index} must be numeric") from exc
            if previous_z is not None:
                if z_number < previous_z:
                    raise WorkflowError(f"{solid_name}: z values must be non-decreasing")
                if z_number > previous_z:
                    has_positive_span = True
            if r_number < 0.0:
                raise WorkflowError(f"{solid_name}: radius values must be non-negative")
            previous_z = z_number
        if not has_positive_span:
            raise WorkflowError(f"{solid_name}: profile must span a non-zero z range")


# ── profile-section editor ────────────────────────────────────────────────────

class _ProfileEditorWidget(QWidget):
    """Edits one solid profile: z/r/rounding table plus optional aperture pattern."""

    def __init__(self, change_callback=None, parent=None):
        super().__init__(parent)
        self._cb = change_callback
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(4, 4, 4, 4)
        root.setSpacing(6)

        # z / r profile table
        profile_box = QGroupBox(
            "Profile points  (z may repeat for steps; rounding row i applies to segment i-1 -> i)"
        )
        pl = QVBoxLayout(profile_box)
        pt_row = QHBoxLayout()
        add_pt = QPushButton("+ Point")
        add_pt.setFixedWidth(80)
        add_pt.clicked.connect(self._add_point)
        rem_pt = QPushButton("− Point")
        rem_pt.setFixedWidth(80)
        rem_pt.clicked.connect(self._remove_point)
        pt_row.addWidget(add_pt)
        pt_row.addWidget(rem_pt)
        pt_row.addStretch()
        pl.addLayout(pt_row)
        self._table = QTableWidget(0, 3)
        self._table.setHorizontalHeaderLabels(["z [m]", "r [m]", "segment rounding [m]"])
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._table.setAlternatingRowColors(True)
        self._table.setMinimumHeight(220)
        self._table.itemChanged.connect(self._fire)
        pl.addWidget(self._table)
        rounding_note = QLabel(
            "Row 0 rounding is ignored. Positive values keep the arc closer to the aperture axis; "
            "negative values keep the farther arc."
        )
        rounding_note.setWordWrap(True)
        pl.addWidget(rounding_note)
        root.addWidget(profile_box)

        # Aperture pattern (collapsible via checkable group box)
        self._ap_box = QGroupBox("Aperture pattern")
        self._ap_box.setCheckable(True)
        self._ap_box.setChecked(False)
        self._ap_box.toggled.connect(self._on_ap_toggled)
        self._ap_box.toggled.connect(self._fire)
        self._ap_content = QWidget()
        af = QFormLayout(self._ap_content)
        self._ap_layout_cb = QComboBox()
        self._ap_layout_cb.addItems(["single", "rectangular-grid", "staggered-grid"])
        self._ap_layout_cb.currentTextChanged.connect(self._on_ap_layout_changed)
        self._ap_layout_cb.currentTextChanged.connect(self._fire)
        self._ap_cx = _ispin(1, 500, 1)
        self._ap_cy = _ispin(1, 500, 1)
        self._ap_px = _dspin(1e-5, 1.0, 6, 1e-3)
        self._ap_px.setValue(0.014)
        self._ap_py = _dspin(1e-5, 1.0, 6, 1e-3)
        self._ap_py.setValue(0.014)
        self._ap_margin = _dspin(0.0, 1.0, 6, 1e-4)
        self._ap_xoff = _dspin(-1.0, 1.0, 6, 1e-4)
        self._ap_yoff = _dspin(-1.0, 1.0, 6, 1e-4)
        self._ap_row_shift = _dspin(-1.0, 1.0, 6, 1e-4)
        self._ap_outside = QCheckBox("Outside pattern is solid")
        self._ap_outside.setChecked(True)
        self._ap_row_shift_lbl = QLabel("Row shift X [m]")
        for w in (self._ap_cx, self._ap_cy, self._ap_px, self._ap_py,
                  self._ap_margin, self._ap_xoff, self._ap_yoff, self._ap_row_shift):
            w.valueChanged.connect(self._fire)
        self._ap_outside.toggled.connect(self._fire)
        af.addRow("Layout", self._ap_layout_cb)
        af.addRow("Count X", self._ap_cx)
        af.addRow("Count Y", self._ap_cy)
        af.addRow("Pitch X [m]", self._ap_px)
        af.addRow("Pitch Y [m]", self._ap_py)
        af.addRow("Margin [m]", self._ap_margin)
        af.addRow("X offset [m]", self._ap_xoff)
        af.addRow("Y offset [m]", self._ap_yoff)
        af.addRow(self._ap_row_shift_lbl, self._ap_row_shift)
        af.addRow("", self._ap_outside)
        ap_outer = QVBoxLayout(self._ap_box)
        ap_outer.setContentsMargins(4, 4, 4, 4)
        ap_outer.addWidget(self._ap_content)
        self._ap_content.setVisible(False)
        root.addWidget(self._ap_box)
        root.addStretch()
        self._on_ap_layout_changed("single")

    def _fire(self, *_) -> None:
        if self._cb:
            self._cb()

    def _on_ap_toggled(self, checked: bool) -> None:
        self._ap_content.setVisible(checked)

    def _on_ap_layout_changed(self, layout: str) -> None:
        is_grid = layout != "single"
        is_staggered = layout == "staggered-grid"
        for w in (self._ap_cx, self._ap_cy, self._ap_px, self._ap_py, self._ap_margin):
            w.setEnabled(is_grid)
        self._ap_row_shift.setEnabled(is_staggered)
        self._ap_row_shift_lbl.setEnabled(is_staggered)

    def _add_point(self) -> None:
        self._table.blockSignals(True)
        row = self._table.rowCount()
        self._table.insertRow(row)
        self._table.setItem(row, 0, QTableWidgetItem("0.0"))
        self._table.setItem(row, 1, QTableWidgetItem("0.01"))
        self._table.setItem(row, 2, QTableWidgetItem("0.0"))
        self._table.blockSignals(False)
        self._fire()

    def _remove_point(self) -> None:
        row = self._table.currentRow()
        if row >= 0:
            self._table.removeRow(row)
            self._fire()

    def set_profile(self, data: dict) -> None:
        self._table.blockSignals(True)
        self._table.setRowCount(0)
        rounding_radii = data.get("roundingRadiiMeters") or []
        if len(rounding_radii) != len(data.get("zProfileMeters", [])):
            rounding_radii = [0.0] * len(data.get("zProfileMeters", []))
        for z, r, rounding in zip(data.get("zProfileMeters", []), data.get("rProfileMeters", []), rounding_radii):
            row = self._table.rowCount()
            self._table.insertRow(row)
            self._table.setItem(row, 0, QTableWidgetItem(str(z)))
            self._table.setItem(row, 1, QTableWidgetItem(str(r)))
            self._table.setItem(row, 2, QTableWidgetItem(str(rounding)))
        self._table.blockSignals(False)
        ap = data.get("aperturePattern")
        self._ap_box.blockSignals(True)
        self._ap_box.setChecked(ap is not None)
        self._ap_box.blockSignals(False)
        self._ap_content.setVisible(ap is not None)
        if ap:
            self._ap_layout_cb.setCurrentText(ap.get("layout", "single"))
            self._ap_cx.setValue(int(ap.get("countX", 1)))
            self._ap_cy.setValue(int(ap.get("countY", 1)))
            self._ap_px.setValue(float(ap.get("pitchXMeters", 0.014)))
            self._ap_py.setValue(float(ap.get("pitchYMeters", 0.014)))
            self._ap_margin.setValue(float(ap.get("marginMeters", 0.0)))
            self._ap_xoff.setValue(float(ap.get("xOffsetMeters", 0.0)))
            self._ap_yoff.setValue(float(ap.get("yOffsetMeters", 0.0)))
            self._ap_row_shift.setValue(float(ap.get("rowShiftXMeters", 0.0)))
            self._ap_outside.setChecked(bool(ap.get("outsidePatternIsSolid", True)))
        else:
            self._ap_layout_cb.setCurrentText("single")
            self._ap_cx.setValue(1)
            self._ap_cy.setValue(1)
            self._ap_px.setValue(0.014)
            self._ap_py.setValue(0.014)
            self._ap_margin.setValue(0.0)
            self._ap_xoff.setValue(0.0)
            self._ap_yoff.setValue(0.0)
            self._ap_row_shift.setValue(0.0)
            self._ap_outside.setChecked(True)

    def get_profile(self) -> dict:
        z_vals: list = []
        r_vals: list = []
        rounding_vals: list = []
        for row in range(self._table.rowCount()):
            zi = self._table.item(row, 0)
            ri = self._table.item(row, 1)
            rounding_item = self._table.item(row, 2)
            try:
                z_vals.append(float(zi.text() if zi else "0"))
                r_vals.append(float(ri.text() if ri else "0"))
                rounding_vals.append(float(rounding_item.text() if rounding_item else "0"))
            except ValueError:
                z_vals.append(0.0)
                r_vals.append(0.0)
                rounding_vals.append(0.0)
        sec: dict = {
            "zProfileMeters": z_vals,
            "rProfileMeters": r_vals,
            "roundingRadiiMeters": rounding_vals,
        }
        if self._ap_box.isChecked():
            layout_name = self._ap_layout_cb.currentText()
            ap: dict = {"layout": layout_name, "outsidePatternIsSolid": self._ap_outside.isChecked()}
            if layout_name != "single":
                ap["countX"] = self._ap_cx.value()
                ap["countY"] = self._ap_cy.value()
                ap["pitchXMeters"] = self._ap_px.value()
                ap["pitchYMeters"] = self._ap_py.value()
                if self._ap_margin.value() > 0.0:
                    ap["marginMeters"] = self._ap_margin.value()
                if abs(self._ap_xoff.value()) > 1e-12:
                    ap["xOffsetMeters"] = self._ap_xoff.value()
                if abs(self._ap_yoff.value()) > 1e-12:
                    ap["yOffsetMeters"] = self._ap_yoff.value()
                if layout_name == "staggered-grid" and abs(self._ap_row_shift.value()) > 1e-12:
                    ap["rowShiftXMeters"] = self._ap_row_shift.value()
            sec["aperturePattern"] = ap
        return sec


# ── solid detail editor ───────────────────────────────────────────────────────

class _SolidDetailWidget(QWidget):
    """Edits one solid: metadata, one profile table, and aperture settings."""

    _KINDS = ["solid", "diagnosticPlane"]

    def __init__(self, change_callback=None, parent=None):
        super().__init__(parent)
        self._cb = change_callback
        self._list_item = None
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(8, 4, 4, 4)
        root.setSpacing(8)

        # Metadata group
        meta_box = QGroupBox("Solid properties")
        mf = QFormLayout(meta_box)
        self._name = QLineEdit()
        self._name.setPlaceholderText("e.g. PG, EG, AG1, AG2 …")
        self._name.textChanged.connect(self._on_name_changed)
        self._kind = QComboBox()
        self._kind.addItems(self._KINDS)
        self._kind.currentTextChanged.connect(self._fire)
        self._boundary_en = QCheckBox("Override  ID =")
        self._boundary_val = _ispin(7, 100000, 7)
        self._boundary_val.setEnabled(False)
        self._boundary_en.toggled.connect(self._boundary_val.setEnabled)
        self._boundary_en.toggled.connect(self._fire)
        self._boundary_val.valueChanged.connect(self._fire)
        boundary_w = QWidget()
        bl = QHBoxLayout(boundary_w)
        bl.setContentsMargins(0, 0, 0, 0)
        bl.addWidget(self._boundary_en)
        bl.addWidget(self._boundary_val)
        bl.addStretch()
        mf.addRow("Name", self._name)
        mf.addRow("Kind", self._kind)
        mf.addRow("Boundary ID", boundary_w)
        root.addWidget(meta_box)

        # Profile
        sec_box = QGroupBox("Profile")
        sl = QVBoxLayout(sec_box)
        self._profile = _ProfileEditorWidget(change_callback=self._cb)
        sl.addWidget(self._profile)
        root.addWidget(sec_box)

    def _fire(self, *_) -> None:
        if self._cb:
            self._cb()

    def _on_name_changed(self, text: str) -> None:
        if self._list_item is not None:
            self._list_item.setText(text.strip() or "unnamed")
        self._fire()

    def set_solid(self, solid: dict, list_item=None) -> None:
        solid = _normalize_solid(solid)
        self._list_item = list_item
        self._name.blockSignals(True)
        self._name.setText(solid.get("name", ""))
        self._name.blockSignals(False)
        self._kind.setCurrentText(solid.get("kind", "solid"))
        boundary_id = solid.get("boundaryId")
        has_boundary_id = isinstance(boundary_id, int) and boundary_id >= 7
        self._boundary_en.setChecked(has_boundary_id)
        self._boundary_val.setValue(boundary_id if has_boundary_id else 7)
        self._profile.set_profile(solid)

    def get_solid(self) -> dict:
        solid: dict = {
            "name": self._name.text().strip() or "unnamed",
            "kind": self._kind.currentText(),
        }
        if self._boundary_en.isChecked():
            solid["boundaryId"] = self._boundary_val.value()
        solid.update(self._profile.get_profile())
        return solid


class _GeometryPreviewWidget(QWidget):
    """Plots the current structured geometry as a profile or simple 3D envelope."""

    _MODES = ["Selected solid profile", "Selected solid 3D", "Full geometry 3D"]

    def __init__(self, parent=None):
        super().__init__(parent)
        self._source_mode = "Built-in template"
        self._template_name = "MTF"
        self._solids: list[dict] = []
        self._selected_solid: dict | None = None
        self._domain: dict[str, float] = {}
        self._view_direction: str | None = None
        self._view_buttons: dict[str, QToolButton] = {}
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(6)

        controls = QHBoxLayout()
        controls.addWidget(QLabel("Preview mode"))
        self._mode = QComboBox()
        self._mode.addItems(self._MODES)
        self._mode.currentTextChanged.connect(self.refresh_preview)
        controls.addWidget(self._mode)
        self._show_solid = QCheckBox("Show solid")
        self._show_solid.toggled.connect(self.refresh_preview)
        controls.addWidget(self._show_solid)
        for direction in ("+x", "-x", "+y", "-y", "+z", "-z"):
            button = QToolButton()
            button.setText(direction.upper())
            button.setCheckable(True)
            button.setFixedWidth(34)
            button.clicked.connect(lambda _checked=False, d=direction: self._set_view_direction(d))
            button.setVisible(False)
            self._view_buttons[direction] = button
            controls.addWidget(button)
        self._full_solids_button = QToolButton()
        self._full_solids_button.setText("Solids")
        self._full_solids_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self._full_solids_menu = QMenu(self._full_solids_button)
        self._full_solids_button.setMenu(self._full_solids_menu)
        self._full_solids_button.setVisible(False)
        controls.addWidget(self._full_solids_button)
        self._full_solids_all = QPushButton("All")
        self._full_solids_all.setFixedWidth(56)
        self._full_solids_all.clicked.connect(self._select_all_full_solids)
        self._full_solids_all.setVisible(False)
        controls.addWidget(self._full_solids_all)
        self._full_solids_none = QPushButton("None")
        self._full_solids_none.setFixedWidth(56)
        self._full_solids_none.clicked.connect(self._clear_all_full_solids)
        self._full_solids_none.setVisible(False)
        controls.addWidget(self._full_solids_none)
        controls.addStretch(1)
        root.addLayout(controls)

        self._canvas = GeometryCanvas()
        self._toolbar = NavigationToolbar2QT(self._canvas, self)
        configure_matplotlib_toolbar(self._toolbar)
        self._status = QLabel("Geometry preview is ready.")
        self._status.setWordWrap(True)
        root.addWidget(self._toolbar)
        root.addWidget(self._canvas, 1)
        root.addWidget(self._status)

    def set_geometry_data(
        self,
        source_mode: str,
        template_name: str,
        solids: list[dict],
        selected_solid: dict | None,
        domain: dict[str, float] | None = None,
    ) -> None:
        self._source_mode = source_mode
        self._template_name = template_name
        self._solids = solids
        self._selected_solid = selected_solid
        self._domain = dict(domain or {})
        self._sync_full_solids_selector()
        self.refresh_preview()

    def _sync_full_solids_selector(self) -> None:
        preserved_states: dict[int, bool] = {}
        for action in self._full_solids_menu.actions():
            index = action.data()
            if isinstance(index, int):
                preserved_states[index] = action.isChecked()

        self._full_solids_menu.clear()
        for index, solid in enumerate(self._solids):
            solid_name = str(solid.get("name", f"Solid {index + 1}"))
            action = QAction(solid_name, self._full_solids_menu)
            action.setData(index)
            action.setCheckable(True)
            action.setChecked(preserved_states.get(index, True))
            action.toggled.connect(self._on_full_solids_action_toggled)
            self._full_solids_menu.addAction(action)
        self._update_full_solids_controls()

    def _update_full_solids_controls(self) -> None:
        actions = self._full_solids_menu.actions()
        action_count = len(actions)
        selected_count = sum(1 for action in actions if action.isChecked())
        is_full_geometry = self._mode.currentText() == "Full geometry 3D" and action_count > 0
        self._full_solids_button.setVisible(is_full_geometry)
        self._full_solids_all.setVisible(is_full_geometry)
        self._full_solids_none.setVisible(is_full_geometry)
        self._full_solids_button.setEnabled(action_count > 0)
        self._full_solids_all.setEnabled(action_count > 0)
        self._full_solids_none.setEnabled(action_count > 0)
        if action_count > 0:
            self._full_solids_button.setText(f"Solids ({selected_count}/{action_count})")
        else:
            self._full_solids_button.setText("Solids")

    def _set_all_full_solids(self, checked: bool) -> None:
        for action in self._full_solids_menu.actions():
            action.blockSignals(True)
            action.setChecked(checked)
            action.blockSignals(False)
        self._update_full_solids_controls()
        self.refresh_preview()

    def _select_all_full_solids(self) -> None:
        self._set_all_full_solids(True)

    def _clear_all_full_solids(self) -> None:
        self._set_all_full_solids(False)

    def _on_full_solids_action_toggled(self, _checked: bool) -> None:
        self._update_full_solids_controls()
        self.refresh_preview()

    def _selected_full_geometry_solids(self) -> list[dict]:
        selected_solids: list[dict] = []
        for action in self._full_solids_menu.actions():
            if not action.isChecked():
                continue
            index = action.data()
            if isinstance(index, int) and 0 <= index < len(self._solids):
                selected_solids.append(self._solids[index])
        return selected_solids

    def _set_view_direction(self, direction: str) -> None:
        self._view_direction = direction
        self._update_3d_view_controls()
        self.refresh_preview()

    def _update_3d_view_controls(self) -> None:
        is_3d_mode = self._mode.currentText() != "Selected solid profile"
        for direction, button in self._view_buttons.items():
            button.setVisible(is_3d_mode)
            button.setEnabled(is_3d_mode)
            button.blockSignals(True)
            button.setChecked(self._view_direction == direction)
            button.blockSignals(False)

    def refresh_preview(self, *_args) -> None:
        self._update_full_solids_controls()
        self._update_3d_view_controls()
        selected_solid = self._selected_solid or (self._solids[0] if self._solids else None)
        if not self._solids:
            if self._source_mode == "Built-in template":
                message = (
                    f"{self._template_name} is using the built-in geometry generator. "
                    "Switch Source to Structured profile solids to preview authored solids here."
                )
            else:
                message = "Add at least one solid to preview the structured geometry."
            self._canvas.plot_placeholder(message)
            self._status.setText(message)
            return

        mode = self._mode.currentText()
        self._show_solid.setEnabled(mode != "Selected solid profile")
        if mode == "Selected solid profile":
            if selected_solid is None:
                message = "Select a solid to preview its radial profile."
                self._canvas.plot_placeholder(message)
                self._status.setText(message)
                return
            status = self._canvas.plot_solid_profile(selected_solid)
            self._status.setText(status)
            return

        if mode == "Selected solid 3D":
            if selected_solid is None:
                message = "Select a solid to preview it in 3D."
                self._canvas.plot_placeholder(message)
                self._status.setText(message)
                return
            status = self._canvas.plot_solid_3d(
                selected_solid,
                domain=self._domain,
                show_solid=self._show_solid.isChecked(),
                view_direction=self._view_direction,
            )
            self._status.setText(status)
            return

        full_geometry_solids = self._selected_full_geometry_solids()
        if not full_geometry_solids:
            message = "Select at least one solid to preview in Full geometry 3D mode."
            self._canvas.plot_placeholder(message)
            self._status.setText(message)
            return

        status = self._canvas.plot_geometry_3d(
            full_geometry_solids,
            domain=self._domain,
            show_solid=self._show_solid.isChecked(),
            view_direction=self._view_direction,
        )
        self._status.setText(f"{status} ({len(full_geometry_solids)} of {len(self._solids)} solids selected)")


class _GeometryWorkspaceWidget(QWidget):
    """Full geometry tab workspace: inputs + solids list, detail editor, and preview."""

    def __init__(self, window, parent=None):
        super().__init__(parent)
        self._window = window
        self._cb = window.schedule_preview_refresh
        self._prev_row = -1
        self._geometry_path: Path | None = None
        self._geometry_name = "Geometry"
        self._build_ui()
        self._sync_geometry_mode()
        self._refresh_preview()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        root.addWidget(splitter, 1)

        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(8)

        inputs_box = QGroupBox("Geometry inputs")
        inputs_form = QFormLayout(inputs_box)
        self._load_geometry_btn = QPushButton("Open")
        self._load_geometry_btn.clicked.connect(self._choose_geometry_file)
        self._geometry_path_edit = QLineEdit()
        self._geometry_path_edit.setReadOnly(True)
        button_row = QHBoxLayout()
        button_row.setContentsMargins(0, 0, 0, 0)
        button_row.setSpacing(6)
        self._save_geometry_btn = QPushButton("Save")
        self._save_geometry_btn.clicked.connect(self._save_geometry_file)
        button_row.addWidget(self._load_geometry_btn)
        button_row.addWidget(self._save_geometry_btn)
        button_row.addStretch(1)
        button_wrapper = QWidget()
        button_wrapper.setLayout(button_row)
        self._geometry_name_edit = QLineEdit()
        self._geometry_name_edit.setText(self._geometry_name)
        self._mesh_size = self._window._double_spin(0.0, 1.0, 6, 0.0005)
        self._domain_x = self._window._double_spin(0.0, 10.0, 6, 0.001)
        self._domain_y = self._window._double_spin(0.0, 10.0, 6, 0.001)
        self._domain_z = self._window._double_spin(0.0, 10.0, 6, 0.001)
        self._domain_z_start = self._window._double_spin(-10.0, 10.0, 6, 0.001)
        self._export_vtk = QCheckBox("Prefer geometry export in VTK")
        inputs_form.addRow("Geometry file", self._geometry_path_edit)
        inputs_form.addRow("", button_wrapper)
        inputs_form.addRow("Name", self._geometry_name_edit)
        inputs_form.addRow("Mesh size [m]", self._mesh_size)
        inputs_form.addRow("Domain x [m]", self._domain_x)
        inputs_form.addRow("Domain y [m]", self._domain_y)
        inputs_form.addRow("Domain z [m]", self._domain_z)
        inputs_form.addRow("Domain z start [m]", self._domain_z_start)
        inputs_form.addRow(self._export_vtk)
        left_layout.addWidget(inputs_box)

        self._solids_box = QGroupBox("Solids")
        solids_layout = QVBoxLayout(self._solids_box)
        solids_buttons = QHBoxLayout()
        add_btn = QPushButton("+ Add")
        add_btn.clicked.connect(self._add_solid)
        remove_btn = QPushButton("− Remove")
        remove_btn.clicked.connect(self._remove_solid)
        solids_buttons.addWidget(add_btn)
        solids_buttons.addWidget(remove_btn)
        solids_layout.addLayout(solids_buttons)
        self._list = QListWidget()
        self._list.setMinimumWidth(140)
        self._list.currentRowChanged.connect(self._on_row_changed)
        solids_layout.addWidget(self._list, 1)
        left_layout.addWidget(self._solids_box, 1)

        self._detail = _SolidDetailWidget(change_callback=self._on_detail_changed)
        self._detail_scroll = QScrollArea()
        self._detail_scroll.setWidget(self._detail)
        self._detail_scroll.setWidgetResizable(True)
        self._detail_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._detail.setEnabled(False)

        self._preview = _GeometryPreviewWidget()

        splitter.addWidget(left_panel)
        splitter.addWidget(self._detail_scroll)
        splitter.addWidget(self._preview)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 1)
        splitter.setSizes([300, 520, 520])

        self._window.widgets["geometry.path"] = self._geometry_path_edit
        self._window.widgets["geometry.name"] = self._geometry_name_edit
        self._window.widgets["geometry.meshSizeMeters"] = self._mesh_size
        self._window.widgets["geometry.domain.x"] = self._domain_x
        self._window.widgets["geometry.domain.y"] = self._domain_y
        self._window.widgets["geometry.domain.z"] = self._domain_z
        self._window.widgets["geometry.domain.zStart"] = self._domain_z_start
        self._window.widgets["geometry.exportGeometryVtk"] = self._export_vtk
        self._window.widgets["geometry.solidsEditor"] = self

        self._geometry_name_edit.textChanged.connect(self._on_geometry_name_changed)
        for widget in (self._mesh_size, self._domain_x, self._domain_y, self._domain_z, self._domain_z_start):
            widget.valueChanged.connect(self._on_geometry_controls_changed)
        self._export_vtk.toggled.connect(self._on_geometry_controls_changed)

    def _fire(self) -> None:
        if self._cb:
            self._cb()
        sync_boundary_editor = getattr(self._window, "sync_boundary_editor_from_geometry", None)
        if callable(sync_boundary_editor):
            sync_boundary_editor()
        sync_diagnostics_grid_power = getattr(self._window, "sync_diagnostics_grid_power_from_geometry", None)
        if callable(sync_diagnostics_grid_power):
            sync_diagnostics_grid_power()
        sync_diagnostics_planes = getattr(self._window, "sync_diagnostics_planes_from_geometry", None)
        if callable(sync_diagnostics_planes):
            sync_diagnostics_planes()

    def _on_geometry_controls_changed(self, *_args) -> None:
        self._fire()
        self._refresh_preview()

    def _on_geometry_name_changed(self, text: str) -> None:
        stripped = text.strip()
        self._geometry_name = stripped or "Geometry"
        self._fire()
        self._refresh_preview()

    def _sync_geometry_mode(self, *_args) -> None:
        self._solids_box.setEnabled(True)
        self._detail_scroll.setEnabled(self._list.currentRow() >= 0)
        self._refresh_preview()

    def _apply_geometry_document(self, geometry_document: dict, source_path: Path | None) -> None:
        self._geometry_path = source_path.resolve() if source_path is not None else None
        self._geometry_name = str(
            geometry_document.get("name") or (source_path.stem if source_path is not None else "Geometry")
        )
        self._geometry_path_edit.setText(
            serialize_geometry_file_path(self._geometry_path) if self._geometry_path is not None else ""
        )
        self._geometry_name_edit.blockSignals(True)
        self._geometry_name_edit.setText(self._geometry_name)
        self._geometry_name_edit.blockSignals(False)
        self._mesh_size.setValue(float(geometry_document.get("meshSizeMeters", 0.002)))
        self._export_vtk.setChecked(bool(geometry_document.get("exportGeometryVtk", True)))
        domain = geometry_document.get("domain", {})
        self._domain_x.setValue(float(domain.get("xSizeMeters", 0.08)))
        self._domain_y.setValue(float(domain.get("ySizeMeters", 0.08)))
        self._domain_z.setValue(float(domain.get("zSizeMeters", 0.567)))
        self._domain_z_start.setValue(float(domain.get("zStartMeters", 0.0)))
        self.set_solids(geometry_document.get("solids", []))
        self._sync_geometry_mode()

    def _current_geometry_document(self) -> dict:
        solids = self.get_solids()
        if not solids:
            solids = []
        geometry_document = {
            "name": self._geometry_name,
            "meshSizeMeters": float(self._mesh_size.value()),
            "exportGeometryVtk": self._export_vtk.isChecked(),
            "domain": {
                "xSizeMeters": float(self._domain_x.value()),
                "ySizeMeters": float(self._domain_y.value()),
                "zSizeMeters": float(self._domain_z.value()),
                "zStartMeters": float(self._domain_z_start.value()),
            },
            "solids": solids,
        }
        return geometry_document

    def current_geometry_document(self) -> dict:
        return copy.deepcopy(self._current_geometry_document())

    def load_geometry_path(self, path: Path) -> None:
        geometry_document = _load_geometry_document(path)
        self._apply_geometry_document(geometry_document, path)
        self._fire()
        self._refresh_preview()

    def save_geometry_path(self, path: Path) -> Path:
        geometry_document = self._current_geometry_document()
        resolved_path = write_geometry_file(path.resolve(), geometry_document)
        self._geometry_path = resolved_path
        self._geometry_path_edit.setText(serialize_geometry_file_path(resolved_path))
        self._fire()
        self._refresh_preview()
        return resolved_path

    def _default_geometry_directory(self) -> Path:
        if self._geometry_path is not None:
            return self._geometry_path.parent
        if TEMPLATES_DIR.is_dir():
            return TEMPLATES_DIR
        return Path.cwd()

    def _choose_geometry_file(self) -> None:
        selected, _ = QFileDialog.getOpenFileName(
            self,
            "Open geometry JSON",
            str(self._default_geometry_directory()),
            "JSON files (*.json)",
        )
        if not selected:
            return
        self.load_geometry_path(Path(selected).resolve())

    def _save_geometry_file(self) -> None:
        suggested_name = self._geometry_path or (self._default_geometry_directory() / f"{self._geometry_name or 'Geometry'}.json")
        selected, _ = QFileDialog.getSaveFileName(
            self,
            "Save geometry JSON",
            str(suggested_name),
            "JSON files (*.json)",
        )
        if not selected:
            return
        self.save_geometry_path(Path(selected).resolve())

    def _flush_current(self) -> None:
        if 0 <= self._prev_row < self._list.count():
            item = self._list.item(self._prev_row)
            updated = self._detail.get_solid()
            item.setData(Qt.ItemDataRole.UserRole, updated)
            item.setText(updated.get("name", "unnamed"))

    def _refresh_preview(self) -> None:
        self._preview.set_geometry_data(
            "Structured profile solids",
            self._geometry_name,
            self.get_solids(),
            self.current_solid(),
            {
                "xSizeMeters": float(self._domain_x.value()),
                "ySizeMeters": float(self._domain_y.value()),
                "zSizeMeters": float(self._domain_z.value()),
                "zStartMeters": float(self._domain_z_start.value()),
            },
        )

    def _on_row_changed(self, row: int) -> None:
        self._flush_current()
        self._prev_row = row
        if row < 0:
            self._detail.setEnabled(False)
            self._detail_scroll.setEnabled(False)
            self._refresh_preview()
            return
        item = self._list.item(row)
        self._detail.setEnabled(True)
        self._detail_scroll.setEnabled(True)
        self._detail.set_solid(item.data(Qt.ItemDataRole.UserRole) or {}, item)
        self._refresh_preview()

    def _on_detail_changed(self, *_args) -> None:
        self._flush_current()
        self._fire()
        self._refresh_preview()

    def _add_solid(self) -> None:
        new_solid = {
            "name": f"SOLID{self._list.count() + 1}",
            "kind": "solid",
            **_default_profile(),
        }
        item = QListWidgetItem(new_solid["name"])
        item.setData(Qt.ItemDataRole.UserRole, new_solid)
        self._list.addItem(item)
        self._list.setCurrentItem(item)
        self._fire()
        self._refresh_preview()

    def _remove_solid(self) -> None:
        row = self._list.currentRow()
        if row < 0:
            return
        self._list.takeItem(row)
        self._prev_row = -1
        if self._list.count() > 0:
            self._list.setCurrentRow(min(row, self._list.count() - 1))
        else:
            self._detail.setEnabled(False)
            self._detail_scroll.setEnabled(False)
            self._refresh_preview()
        self._fire()

    def current_solid(self) -> dict | None:
        self._flush_current()
        item = self._list.currentItem()
        if item is None:
            return None
        data = item.data(Qt.ItemDataRole.UserRole)
        return data if isinstance(data, dict) else None

    def set_solids(self, solids: list) -> None:
        self._list.blockSignals(True)
        self._list.clear()
        for solid in solids:
            normalized = _normalize_solid(solid)
            item = QListWidgetItem(normalized.get("name", "unnamed"))
            item.setData(Qt.ItemDataRole.UserRole, normalized)
            self._list.addItem(item)
        self._list.blockSignals(False)
        self._prev_row = -1
        if self._list.count() > 0:
            self._list.setCurrentRow(0)
        else:
            self._detail.setEnabled(False)
            self._detail_scroll.setEnabled(False)
        self._refresh_preview()

    def get_solids(self) -> list:
        self._flush_current()
        solids: list = []
        for index in range(self._list.count()):
            item = self._list.item(index)
            data = item.data(Qt.ItemDataRole.UserRole)
            if data:
                solids.append(data)
        return solids


# ── master-detail solid list ──────────────────────────────────────────────────

class _SolidsEditorWidget(QWidget):
    """Left: solid list; right: solid detail form."""

    def __init__(self, change_callback=None, parent=None):
        super().__init__(parent)
        self._cb = change_callback
        self._prev_row = -1
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        root.addWidget(splitter)

        # Left: list + buttons
        left = QWidget()
        ll = QVBoxLayout(left)
        ll.setContentsMargins(0, 0, 4, 0)
        ll.setSpacing(4)
        btn_row = QHBoxLayout()
        add_btn = QPushButton("+ Add")
        add_btn.clicked.connect(self._add_solid)
        rem_btn = QPushButton("− Remove")
        rem_btn.clicked.connect(self._remove_solid)
        btn_row.addWidget(add_btn)
        btn_row.addWidget(rem_btn)
        ll.addLayout(btn_row)
        self._list = QListWidget()
        self._list.setMinimumWidth(100)
        self._list.setMaximumWidth(220)
        self._list.currentRowChanged.connect(self._on_row_changed)
        ll.addWidget(self._list)

        # Right: detail widget inside a scroll area
        self._detail = _SolidDetailWidget(change_callback=self._cb)
        detail_scroll = QScrollArea()
        detail_scroll.setWidget(self._detail)
        detail_scroll.setWidgetResizable(True)
        detail_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        splitter.addWidget(left)
        splitter.addWidget(detail_scroll)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([130, 500])
        self.setMinimumHeight(520)

    def _fire(self) -> None:
        if self._cb:
            self._cb()

    def _flush_current(self) -> None:
        """Save the detail form's state back into the currently selected list item."""
        if 0 <= self._prev_row < self._list.count():
            item = self._list.item(self._prev_row)
            updated = self._detail.get_solid()
            item.setData(Qt.ItemDataRole.UserRole, updated)
            item.setText(updated.get("name", "unnamed"))

    def _on_row_changed(self, row: int) -> None:
        self._flush_current()
        self._prev_row = row
        if row < 0:
            return
        item = self._list.item(row)
        self._detail.set_solid(item.data(Qt.ItemDataRole.UserRole) or {}, item)

    def _add_solid(self) -> None:
        new_solid = {
            "name": f"SOLID{self._list.count() + 1}",
            "kind": "solid",
            **_default_profile(),
        }
        item = QListWidgetItem(new_solid["name"])
        item.setData(Qt.ItemDataRole.UserRole, new_solid)
        self._list.addItem(item)
        self._list.setCurrentItem(item)
        self._fire()

    def _remove_solid(self) -> None:
        row = self._list.currentRow()
        if row < 0:
            return
        self._list.takeItem(row)
        self._prev_row = -1
        self._fire()

    def set_solids(self, solids: list) -> None:
        self._list.blockSignals(True)
        self._list.clear()
        for s in solids:
            item = QListWidgetItem(s.get("name", "unnamed"))
            item.setData(Qt.ItemDataRole.UserRole, s)
            self._list.addItem(item)
        self._list.blockSignals(False)
        self._prev_row = -1
        if self._list.count() > 0:
            self._list.setCurrentRow(0)

    def get_solids(self) -> list:
        self._flush_current()
        result = []
        for i in range(self._list.count()):
            item = self._list.item(i)
            data = item.data(Qt.ItemDataRole.UserRole)
            if data:
                result.append(data)
        return result


# ── section API ───────────────────────────────────────────────────────────────

def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    layout.addRow(build_workspace(window))
    return layout


def build_workspace(window) -> QWidget:
    return _GeometryWorkspaceWidget(window)


def populate(window, spec: dict) -> None:
    geometry = nested_get(spec, "geometry", default={})
    if not isinstance(geometry, dict):
        geometry = {}
    geometry_path_value = geometry.get("path")
    workspace = window.widgets["geometry.solidsEditor"]
    if isinstance(geometry_path_value, str) and geometry_path_value.strip():
        base_path = window.authoring_path if getattr(window, "authoring_path", None) is not None else None
        path, _geometry_document = load_geometry_reference(geometry_path_value, base_path)
        workspace.load_geometry_path(path)
        return

    solids = geometry.get("solids", [])
    normalized_solids = [_normalize_solid(solid) for solid in solids] if isinstance(solids, list) else []
    geometry_name = str(nested_get(spec, "geometry", "name", default="Geometry"))
    window.widgets["geometry.meshSizeMeters"].setValue(float(nested_get(spec, "geometry", "meshSizeMeters", default=0.002)))
    window.widgets["geometry.domain.x"].setValue(float(nested_get(spec, "geometry", "domain", "xSizeMeters", default=0.08)))
    window.widgets["geometry.domain.y"].setValue(float(nested_get(spec, "geometry", "domain", "ySizeMeters", default=0.08)))
    window.widgets["geometry.domain.z"].setValue(float(nested_get(spec, "geometry", "domain", "zSizeMeters", default=0.567)))
    window.widgets["geometry.domain.zStart"].setValue(float(nested_get(spec, "geometry", "domain", "zStartMeters", default=0.0)))
    window.widgets["geometry.exportGeometryVtk"].setChecked(bool(nested_get(spec, "geometry", "exportGeometryVtk", default=True)))
    workspace._geometry_path = None
    workspace._geometry_name = geometry_name.strip() or "Geometry"
    window.widgets["geometry.path"].setText("")
    window.widgets["geometry.name"].setText(workspace._geometry_name)
    workspace.set_solids(normalized_solids)


def collect(window, spec: dict) -> None:
    geometry = spec.setdefault("geometry", {})
    workspace = window.widgets["geometry.solidsEditor"]
    if workspace._geometry_path is None:
        raise WorkflowError("Open a geometry JSON, then save it to a geometry file before saving authoring JSON")
    geometry.clear()
    geometry["path"] = serialize_geometry_file_path(workspace._geometry_path)
