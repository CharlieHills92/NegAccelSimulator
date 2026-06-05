"""Geometry form section for the NegAccel GUI."""

from __future__ import annotations

from ..common import GEOMETRY_TEMPLATES, QCheckBox, QComboBox, QFormLayout, nested_get, set_combo_value


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    template = QComboBox()
    template.addItems(GEOMETRY_TEMPLATES)
    window.widgets["geometry.template"] = template
    window.widgets["geometry.meshSizeMeters"] = window._double_spin(0.0, 1.0, 6, 0.0005)
    window.widgets["geometry.domain.x"] = window._double_spin(0.0, 10.0, 6, 0.001)
    window.widgets["geometry.domain.y"] = window._double_spin(0.0, 10.0, 6, 0.001)
    window.widgets["geometry.domain.z"] = window._double_spin(0.0, 10.0, 6, 0.001)
    window.widgets["geometry.gaps.accelerationGapMeters"] = window._double_spin(0.0, 2.0, 6, 0.001)
    window.widgets["geometry.exportGeometryVtk"] = QCheckBox("Prefer geometry export in VTK")
    layout.addRow("Template", window.widgets["geometry.template"])
    layout.addRow("Mesh size [m]", window.widgets["geometry.meshSizeMeters"])
    layout.addRow("Domain x [m]", window.widgets["geometry.domain.x"])
    layout.addRow("Domain y [m]", window.widgets["geometry.domain.y"])
    layout.addRow("Domain z [m]", window.widgets["geometry.domain.z"])
    layout.addRow("Acceleration gap [m]", window.widgets["geometry.gaps.accelerationGapMeters"])
    layout.addRow(window.widgets["geometry.exportGeometryVtk"])
    return layout


def populate(window, spec: dict[str, object]) -> None:
    set_combo_value(window.widgets["geometry.template"], str(nested_get(spec, "geometry", "template", default="MTF")))
    window.widgets["geometry.meshSizeMeters"].setValue(float(nested_get(spec, "geometry", "meshSizeMeters", default=0.002)))
    window.widgets["geometry.domain.x"].setValue(float(nested_get(spec, "geometry", "domain", "xSizeMeters", default=0.08)))
    window.widgets["geometry.domain.y"].setValue(float(nested_get(spec, "geometry", "domain", "ySizeMeters", default=0.08)))
    window.widgets["geometry.domain.z"].setValue(float(nested_get(spec, "geometry", "domain", "zSizeMeters", default=0.567)))
    window.widgets["geometry.gaps.accelerationGapMeters"].setValue(
        float(nested_get(spec, "geometry", "gaps", "accelerationGapMeters", default=0.0))
    )
    window.widgets["geometry.exportGeometryVtk"].setChecked(
        bool(nested_get(spec, "geometry", "exportGeometryVtk", default=True))
    )


def collect(window, spec: dict[str, object]) -> None:
    geometry = spec.setdefault("geometry", {})
    geometry["template"] = window.widgets["geometry.template"].currentText()
    geometry["meshSizeMeters"] = float(window.widgets["geometry.meshSizeMeters"].value())
    geometry["exportGeometryVtk"] = window.widgets["geometry.exportGeometryVtk"].isChecked()
    domain = geometry.setdefault("domain", {})
    domain["xSizeMeters"] = float(window.widgets["geometry.domain.x"].value())
    domain["ySizeMeters"] = float(window.widgets["geometry.domain.y"].value())
    domain["zSizeMeters"] = float(window.widgets["geometry.domain.z"].value())
    gaps = geometry.setdefault("gaps", {})
    gaps["accelerationGapMeters"] = float(window.widgets["geometry.gaps.accelerationGapMeters"].value())
