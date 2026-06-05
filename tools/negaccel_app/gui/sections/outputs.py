"""Output and logging form section for the NegAccel GUI."""

from __future__ import annotations

from ..common import LOG_LEVELS, QCheckBox, QComboBox, QFormLayout, QLineEdit, nested_get, set_combo_value


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    window.widgets["outputs.summary.enabled"] = QCheckBox("Write summary outputs")
    window.widgets["outputs.summary.directory"] = QLineEdit()
    window.widgets["outputs.plots.enabled"] = QCheckBox("Write plots")
    window.widgets["outputs.plots.directory"] = QLineEdit()
    window.widgets["outputs.data.enabled"] = QCheckBox("Write .dat state files")
    window.widgets["outputs.data.directory"] = QLineEdit()
    window.widgets["outputs.vtk.enabled"] = QCheckBox("Write VTK outputs")
    window.widgets["outputs.vtk.directory"] = QLineEdit()
    window.widgets["outputs.vtk.exportGeometry"] = QCheckBox("Export geometry")
    window.widgets["outputs.vtk.exportSimulationState"] = QCheckBox("Export simulation state")
    window.widgets["outputs.vtk.exportTracedParticles"] = QCheckBox("Export traced particles")
    console_level = QComboBox()
    console_level.addItems(LOG_LEVELS)
    window.widgets["outputs.logging.consoleLevel"] = console_level
    file_level = QComboBox()
    file_level.addItems(LOG_LEVELS)
    window.widgets["outputs.logging.fileLevel"] = file_level
    window.widgets["outputs.logging.captureStdout"] = QCheckBox(
        "Capture stdout/stderr into the configured log file"
    )
    window.widgets["outputs.logging.writeDebugArtifacts"] = QCheckBox("Enable debug-only artifacts")
    window.widgets["outputs.logging.structuredLogFile"] = QLineEdit()
    layout.addRow(window.widgets["outputs.summary.enabled"])
    layout.addRow("Summary directory", window.widgets["outputs.summary.directory"])
    layout.addRow(window.widgets["outputs.plots.enabled"])
    layout.addRow("Plots directory", window.widgets["outputs.plots.directory"])
    layout.addRow(window.widgets["outputs.data.enabled"])
    layout.addRow("Data directory", window.widgets["outputs.data.directory"])
    layout.addRow(window.widgets["outputs.vtk.enabled"])
    layout.addRow("VTK directory", window.widgets["outputs.vtk.directory"])
    layout.addRow(window.widgets["outputs.vtk.exportGeometry"])
    layout.addRow(window.widgets["outputs.vtk.exportSimulationState"])
    layout.addRow(window.widgets["outputs.vtk.exportTracedParticles"])
    layout.addRow("Console log level", window.widgets["outputs.logging.consoleLevel"])
    layout.addRow("File log level", window.widgets["outputs.logging.fileLevel"])
    layout.addRow(window.widgets["outputs.logging.captureStdout"])
    layout.addRow(window.widgets["outputs.logging.writeDebugArtifacts"])
    layout.addRow("Structured log file", window.widgets["outputs.logging.structuredLogFile"])
    return layout


def populate(window, spec: dict[str, object]) -> None:
    window.widgets["outputs.summary.enabled"].setChecked(bool(nested_get(spec, "outputs", "summary", "enabled", default=True)))
    window.widgets["outputs.summary.directory"].setText(
        str(nested_get(spec, "outputs", "summary", "directory", default="Summary"))
    )
    window.widgets["outputs.plots.enabled"].setChecked(bool(nested_get(spec, "outputs", "plots", "enabled", default=True)))
    window.widgets["outputs.plots.directory"].setText(
        str(nested_get(spec, "outputs", "plots", "directory", default="Plots"))
    )
    window.widgets["outputs.data.enabled"].setChecked(bool(nested_get(spec, "outputs", "data", "enabled", default=True)))
    window.widgets["outputs.data.directory"].setText(
        str(nested_get(spec, "outputs", "data", "directory", default="Data"))
    )
    window.widgets["outputs.vtk.enabled"].setChecked(bool(nested_get(spec, "outputs", "vtk", "enabled", default=True)))
    window.widgets["outputs.vtk.directory"].setText(
        str(nested_get(spec, "outputs", "vtk", "directory", default="VTK"))
    )
    window.widgets["outputs.vtk.exportGeometry"].setChecked(
        bool(nested_get(spec, "outputs", "vtk", "exportGeometry", default=True))
    )
    window.widgets["outputs.vtk.exportSimulationState"].setChecked(
        bool(nested_get(spec, "outputs", "vtk", "exportSimulationState", default=True))
    )
    window.widgets["outputs.vtk.exportTracedParticles"].setChecked(
        bool(nested_get(spec, "outputs", "vtk", "exportTracedParticles", default=True))
    )
    set_combo_value(
        window.widgets["outputs.logging.consoleLevel"],
        str(nested_get(spec, "outputs", "logging", "consoleLevel", default="info")),
    )
    set_combo_value(
        window.widgets["outputs.logging.fileLevel"],
        str(nested_get(spec, "outputs", "logging", "fileLevel", default="debug")),
    )
    window.widgets["outputs.logging.captureStdout"].setChecked(
        bool(nested_get(spec, "outputs", "logging", "captureStdout", default=True))
    )
    window.widgets["outputs.logging.writeDebugArtifacts"].setChecked(
        bool(nested_get(spec, "outputs", "logging", "writeDebugArtifacts", default=False))
    )
    window.widgets["outputs.logging.structuredLogFile"].setText(
        str(nested_get(spec, "outputs", "logging", "structuredLogFile", default="run.log"))
    )


def collect(window, spec: dict[str, object]) -> None:
    outputs = spec.setdefault("outputs", {})
    summary = outputs.setdefault("summary", {})
    summary["enabled"] = window.widgets["outputs.summary.enabled"].isChecked()
    summary["directory"] = window.widgets["outputs.summary.directory"].text().strip() or "Summary"
    plots = outputs.setdefault("plots", {})
    plots["enabled"] = window.widgets["outputs.plots.enabled"].isChecked()
    plots["directory"] = window.widgets["outputs.plots.directory"].text().strip() or "Plots"
    data = outputs.setdefault("data", {})
    data["enabled"] = window.widgets["outputs.data.enabled"].isChecked()
    data["directory"] = window.widgets["outputs.data.directory"].text().strip() or "Data"
    vtk = outputs.setdefault("vtk", {})
    vtk["enabled"] = window.widgets["outputs.vtk.enabled"].isChecked()
    vtk["directory"] = window.widgets["outputs.vtk.directory"].text().strip() or "VTK"
    vtk["exportGeometry"] = window.widgets["outputs.vtk.exportGeometry"].isChecked()
    vtk["exportSimulationState"] = window.widgets["outputs.vtk.exportSimulationState"].isChecked()
    vtk["exportTracedParticles"] = window.widgets["outputs.vtk.exportTracedParticles"].isChecked()
    logging = outputs.setdefault("logging", {})
    logging["consoleLevel"] = window.widgets["outputs.logging.consoleLevel"].currentText()
    logging["fileLevel"] = window.widgets["outputs.logging.fileLevel"].currentText()
    logging["captureStdout"] = window.widgets["outputs.logging.captureStdout"].isChecked()
    logging["writeDebugArtifacts"] = window.widgets["outputs.logging.writeDebugArtifacts"].isChecked()
    logging["structuredLogFile"] = window.widgets["outputs.logging.structuredLogFile"].text().strip() or "run.log"