"""Boundary-condition form section for the NegAccel GUI."""

from __future__ import annotations

from ..common import QCheckBox, QFormLayout, QLineEdit, WorkflowError, format_number_list, nested_get, parse_number_list


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    window.widgets["boundary.extractionGrid"] = window._double_spin(-1.0e7, 1.0e7, 3, 100.0)
    window.widgets["boundary.acceleratorStages"] = QLineEdit()
    window.widgets["boundary.acceleratorStages"].setPlaceholderText(
        "Comma-separated voltages, e.g. 165600, 324000, 482400"
    )
    window.widgets["boundary.periodic.enabled"] = QCheckBox("Enable periodic boundaries")
    window.widgets["boundary.periodic.xMin"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["boundary.periodic.xMax"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["boundary.periodic.yMin"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["boundary.periodic.yMax"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    layout.addRow("Extraction grid [V]", window.widgets["boundary.extractionGrid"])
    layout.addRow("Accelerator stages [V]", window.widgets["boundary.acceleratorStages"])
    layout.addRow(window.widgets["boundary.periodic.enabled"])
    layout.addRow("Periodic x min [m]", window.widgets["boundary.periodic.xMin"])
    layout.addRow("Periodic x max [m]", window.widgets["boundary.periodic.xMax"])
    layout.addRow("Periodic y min [m]", window.widgets["boundary.periodic.yMin"])
    layout.addRow("Periodic y max [m]", window.widgets["boundary.periodic.yMax"])
    return layout


def populate(window, spec: dict[str, object]) -> None:
    window.widgets["boundary.extractionGrid"].setValue(
        float(nested_get(spec, "boundaryConditions", "gridVoltagesVolts", "extractionGrid", default=0.0))
    )
    window.widgets["boundary.acceleratorStages"].setText(
        format_number_list(
            list(nested_get(spec, "boundaryConditions", "gridVoltagesVolts", "acceleratorStages", default=[]))
        )
    )
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
    voltages = boundary.setdefault("gridVoltagesVolts", {})
    voltages["extractionGrid"] = float(window.widgets["boundary.extractionGrid"].value())
    stages = parse_number_list(
        window.widgets["boundary.acceleratorStages"].text(),
        "boundaryConditions.gridVoltagesVolts.acceleratorStages",
    )
    if not stages:
        raise WorkflowError("boundaryConditions.gridVoltagesVolts.acceleratorStages must not be empty")
    voltages["acceleratorStages"] = stages
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
