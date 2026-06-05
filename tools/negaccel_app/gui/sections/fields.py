"""Magnetic-field form section for the NegAccel GUI."""

from __future__ import annotations

from ..common import MAGNETIC_SOURCE_MODES, QCheckBox, QComboBox, QFormLayout, QLineEdit, QSpinBox, nested_get, set_combo_value


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    window.widgets["magnetic.enabled"] = QCheckBox("Enable external magnetic field")
    source_mode = QComboBox()
    source_mode.addItems(MAGNETIC_SOURCE_MODES)
    window.widgets["magnetic.sourceMode"] = source_mode
    window.widgets["magnetic.case"] = QSpinBox()
    window.widgets["magnetic.case"].setRange(0, 9999)
    window.widgets["magnetic.scale"] = window._double_spin(0.0, 1.0e6, 6, 0.1)
    window.widgets["magnetic.directory"] = QLineEdit()
    window.widgets["magnetic.file"] = QLineEdit()
    window.widgets["magnetic.filePattern"] = QLineEdit()
    layout.addRow(window.widgets["magnetic.enabled"])
    layout.addRow("Source mode", window.widgets["magnetic.sourceMode"])
    layout.addRow("Field case index", window.widgets["magnetic.case"])
    layout.addRow("Field scale", window.widgets["magnetic.scale"])
    layout.addRow("Directory", window.widgets["magnetic.directory"])
    layout.addRow("File", window.widgets["magnetic.file"])
    layout.addRow("File pattern", window.widgets["magnetic.filePattern"])
    return layout


def populate(window, spec: dict[str, object]) -> None:
    window.widgets["magnetic.enabled"].setChecked(bool(nested_get(spec, "magneticField", "enabled", default=False)))
    set_combo_value(
        window.widgets["magnetic.sourceMode"],
        str(nested_get(spec, "magneticField", "sourceMode", default="auto-by-accelerator")),
    )
    window.widgets["magnetic.case"].setValue(int(nested_get(spec, "magneticField", "case", default=1)))
    window.widgets["magnetic.scale"].setValue(float(nested_get(spec, "magneticField", "scale", default=1.0)))
    window.widgets["magnetic.directory"].setText(str(nested_get(spec, "magneticField", "directory", default="")))
    window.widgets["magnetic.file"].setText(str(nested_get(spec, "magneticField", "file", default="")))
    window.widgets["magnetic.filePattern"].setText(str(nested_get(spec, "magneticField", "filePattern", default="")))


def collect(window, spec: dict[str, object]) -> None:
    magnetic = spec.setdefault("magneticField", {})
    magnetic["enabled"] = window.widgets["magnetic.enabled"].isChecked()
    magnetic["sourceMode"] = window.widgets["magnetic.sourceMode"].currentText()
    magnetic["case"] = int(window.widgets["magnetic.case"].value())
    magnetic["scale"] = float(window.widgets["magnetic.scale"].value())
    window._set_optional_text(magnetic, "directory", window.widgets["magnetic.directory"].text())
    window._set_optional_text(magnetic, "file", window.widgets["magnetic.file"].text())
    window._set_optional_text(magnetic, "filePattern", window.widgets["magnetic.filePattern"].text())
