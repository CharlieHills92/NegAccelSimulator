"""Output and logging form section for the NegAccel GUI."""

from __future__ import annotations

from ..common import (
    LOG_LEVELS,
    QCheckBox,
    QComboBox,
    QFormLayout,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QWidget,
    format_number_list,
    nested_get,
    parse_number_list,
    set_combo_value,
)


def _copy_iteration_planes_from_diagnostics(window) -> None:
    diagnostics_widget = window.widgets.get("diagnostics.sampleZPositionsMeters")
    iteration_widget = window.widgets.get("outputs.iteration.planeZPositionsMeters")
    if not isinstance(diagnostics_widget, QLineEdit) or not isinstance(iteration_widget, QLineEdit):
        return
    iteration_widget.setText(diagnostics_widget.text().strip())
    window.schedule_preview_refresh()


def _clear_iteration_plane_override(window) -> None:
    iteration_widget = window.widgets.get("outputs.iteration.planeZPositionsMeters")
    if not isinstance(iteration_widget, QLineEdit):
        return
    iteration_widget.clear()
    window.schedule_preview_refresh()


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
    window.widgets["outputs.iteration.enabled"] = QCheckBox("Write per-iteration outputs")
    iteration_stride = QSpinBox()
    iteration_stride.setRange(1, 100000)
    iteration_stride.setValue(1)
    window.widgets["outputs.iteration.everyNIterations"] = iteration_stride
    window.widgets["outputs.iteration.exportPlaneDiagnostics"] = QCheckBox(
        "Write per-iteration plane diagnostics"
    )
    window.widgets["outputs.iteration.exportSimulationState"] = QCheckBox(
        "Export per-iteration fields to VTK"
    )
    window.widgets["outputs.iteration.exportTracedParticles"] = QCheckBox(
        "Export per-iteration trajectories to VTK"
    )
    window.widgets["outputs.iteration.planeZPositionsMeters"] = QLineEdit()
    window.widgets["outputs.iteration.planeZPositionsMeters"].setToolTip(
        "Optional comma-separated z positions in meters. Leave blank to keep the legacy EG/OUT iteration summaries."
    )
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
    window.widgets["diagnostics.meniscusEnabled"] = QCheckBox("Render meniscus zoom plots")
    window.widgets["diagnostics.meniscusZMinMeters"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["diagnostics.meniscusZMaxMeters"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["diagnostics.meniscusTransverseMinMeters"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    window.widgets["diagnostics.meniscusTransverseMaxMeters"] = window._double_spin(-10.0, 10.0, 6, 0.001)
    copy_planes_button = QPushButton("Use Diagnostics Sample Planes")
    copy_planes_button.clicked.connect(lambda: _copy_iteration_planes_from_diagnostics(window))
    clear_planes_button = QPushButton("Use Legacy EG/OUT Defaults")
    clear_planes_button.clicked.connect(lambda: _clear_iteration_plane_override(window))
    iteration_plane_button_row = QHBoxLayout()
    iteration_plane_button_row.addWidget(copy_planes_button)
    iteration_plane_button_row.addWidget(clear_planes_button)
    iteration_plane_button_row.addStretch(1)
    iteration_plane_button_widget = QWidget()
    iteration_plane_button_widget.setLayout(iteration_plane_button_row)
    layout.addRow(window.widgets["outputs.summary.enabled"])
    layout.addRow("Summary directory", window.widgets["outputs.summary.directory"])
    layout.addRow(window.widgets["outputs.plots.enabled"])
    layout.addRow("Plots directory", window.widgets["outputs.plots.directory"])
    layout.addRow(window.widgets["diagnostics.meniscusEnabled"])
    layout.addRow("Meniscus z min [m]", window.widgets["diagnostics.meniscusZMinMeters"])
    layout.addRow("Meniscus z max [m]", window.widgets["diagnostics.meniscusZMaxMeters"])
    layout.addRow("Meniscus transverse min [m]", window.widgets["diagnostics.meniscusTransverseMinMeters"])
    layout.addRow("Meniscus transverse max [m]", window.widgets["diagnostics.meniscusTransverseMaxMeters"])
    layout.addRow(window.widgets["outputs.data.enabled"])
    layout.addRow("Data directory", window.widgets["outputs.data.directory"])
    layout.addRow(window.widgets["outputs.vtk.enabled"])
    layout.addRow("VTK directory", window.widgets["outputs.vtk.directory"])
    layout.addRow(window.widgets["outputs.vtk.exportGeometry"])
    layout.addRow(window.widgets["outputs.vtk.exportSimulationState"])
    layout.addRow(window.widgets["outputs.vtk.exportTracedParticles"])
    layout.addRow(window.widgets["outputs.iteration.enabled"])
    layout.addRow("Iteration output stride", window.widgets["outputs.iteration.everyNIterations"])
    layout.addRow(window.widgets["outputs.iteration.exportPlaneDiagnostics"])
    layout.addRow(window.widgets["outputs.iteration.exportSimulationState"])
    layout.addRow(window.widgets["outputs.iteration.exportTracedParticles"])
    layout.addRow("Iteration plane overrides [m]", window.widgets["outputs.iteration.planeZPositionsMeters"])
    layout.addRow("", iteration_plane_button_widget)
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
    window.widgets["diagnostics.meniscusEnabled"].setChecked(
        bool(nested_get(spec, "diagnostics", "plots", "meniscus", "enabled", default=True))
    )
    window.widgets["diagnostics.meniscusZMinMeters"].setValue(
        float(nested_get(spec, "diagnostics", "plots", "meniscus", "zMinMeters", default=0.0))
    )
    window.widgets["diagnostics.meniscusZMaxMeters"].setValue(
        float(nested_get(spec, "diagnostics", "plots", "meniscus", "zMaxMeters", default=0.034))
    )
    window.widgets["diagnostics.meniscusTransverseMinMeters"].setValue(
        float(nested_get(spec, "diagnostics", "plots", "meniscus", "transverseMinMeters", default=-0.01))
    )
    window.widgets["diagnostics.meniscusTransverseMaxMeters"].setValue(
        float(nested_get(spec, "diagnostics", "plots", "meniscus", "transverseMaxMeters", default=0.01))
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
    window.widgets["outputs.iteration.enabled"].setChecked(
        bool(nested_get(spec, "outputs", "iteration", "enabled", default=True))
    )
    window.widgets["outputs.iteration.everyNIterations"].setValue(
        int(nested_get(spec, "outputs", "iteration", "everyNIterations", default=1))
    )
    window.widgets["outputs.iteration.exportPlaneDiagnostics"].setChecked(
        bool(nested_get(spec, "outputs", "iteration", "exportPlaneDiagnostics", default=True))
    )
    window.widgets["outputs.iteration.exportSimulationState"].setChecked(
        bool(nested_get(spec, "outputs", "iteration", "exportSimulationState", default=False))
    )
    window.widgets["outputs.iteration.exportTracedParticles"].setChecked(
        bool(nested_get(spec, "outputs", "iteration", "exportTracedParticles", default=False))
    )
    iteration_planes = nested_get(spec, "outputs", "iteration", "planeZPositionsMeters", default=[])
    if isinstance(iteration_planes, list):
        window.widgets["outputs.iteration.planeZPositionsMeters"].setText(
            format_number_list(iteration_planes)
        )
    else:
        window.widgets["outputs.iteration.planeZPositionsMeters"].clear()
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
    diagnostics = spec.setdefault("diagnostics", {})
    diagnostic_plots = diagnostics.setdefault("plots", {})
    meniscus = diagnostic_plots.setdefault("meniscus", {})
    meniscus["enabled"] = window.widgets["diagnostics.meniscusEnabled"].isChecked()
    meniscus["zMinMeters"] = float(window.widgets["diagnostics.meniscusZMinMeters"].value())
    meniscus["zMaxMeters"] = float(window.widgets["diagnostics.meniscusZMaxMeters"].value())
    meniscus["transverseMinMeters"] = float(
        window.widgets["diagnostics.meniscusTransverseMinMeters"].value()
    )
    meniscus["transverseMaxMeters"] = float(
        window.widgets["diagnostics.meniscusTransverseMaxMeters"].value()
    )
    data = outputs.setdefault("data", {})
    data["enabled"] = window.widgets["outputs.data.enabled"].isChecked()
    data["directory"] = window.widgets["outputs.data.directory"].text().strip() or "Data"
    vtk = outputs.setdefault("vtk", {})
    vtk["enabled"] = window.widgets["outputs.vtk.enabled"].isChecked()
    vtk["directory"] = window.widgets["outputs.vtk.directory"].text().strip() or "VTK"
    vtk["exportGeometry"] = window.widgets["outputs.vtk.exportGeometry"].isChecked()
    vtk["exportSimulationState"] = window.widgets["outputs.vtk.exportSimulationState"].isChecked()
    vtk["exportTracedParticles"] = window.widgets["outputs.vtk.exportTracedParticles"].isChecked()
    iteration = outputs.setdefault("iteration", {})
    iteration["enabled"] = window.widgets["outputs.iteration.enabled"].isChecked()
    iteration["everyNIterations"] = int(window.widgets["outputs.iteration.everyNIterations"].value())
    iteration["exportPlaneDiagnostics"] = window.widgets[
        "outputs.iteration.exportPlaneDiagnostics"
    ].isChecked()
    iteration["exportSimulationState"] = window.widgets[
        "outputs.iteration.exportSimulationState"
    ].isChecked()
    iteration["exportTracedParticles"] = window.widgets[
        "outputs.iteration.exportTracedParticles"
    ].isChecked()
    iteration["planeZPositionsMeters"] = parse_number_list(
        window.widgets["outputs.iteration.planeZPositionsMeters"].text(),
        "outputs.iteration.planeZPositionsMeters",
    )
    logging = outputs.setdefault("logging", {})
    logging["consoleLevel"] = window.widgets["outputs.logging.consoleLevel"].currentText()
    logging["fileLevel"] = window.widgets["outputs.logging.fileLevel"].currentText()
    logging["captureStdout"] = window.widgets["outputs.logging.captureStdout"].isChecked()
    logging["writeDebugArtifacts"] = window.widgets["outputs.logging.writeDebugArtifacts"].isChecked()
    logging["structuredLogFile"] = window.widgets["outputs.logging.structuredLogFile"].text().strip() or "run.log"