"""Process execution helpers for the NegAccel GUI."""

from __future__ import annotations

import json
import shlex
from pathlib import Path

from .common import (
    DEFAULT_SETUP_SCRIPT,
    DEFAULT_SIMULATOR_PATH,
    QMessageBox,
    QProcess,
    QTextCursor,
    REPO_ROOT,
    nested_get,
    write_json,
)


class ExecutionMixin:
    def materialize_runtime_case(self) -> Path | None:
        try:
            authoring_spec = self.build_authoring_spec()
            runtime_case, runtime_path = self.build_runtime_case_preview(authoring_spec)
            write_json(runtime_path, runtime_case)
        except Exception as exc:  # pragma: no cover - UI error path
            QMessageBox.critical(self, "Materialization failed", str(exc))
            return None

        self.current_runtime_case = runtime_case
        self.current_runtime_path = runtime_path
        self.authoring_preview.setPlainText(json.dumps(authoring_spec, indent=2))
        self.runtime_preview.setPlainText(json.dumps(runtime_case, indent=2))
        self.append_log(f"Materialized runtime case: {runtime_path}")
        self.refresh_vtk_files()
        self.statusBar().showMessage(f"Materialized runtime case to {runtime_path}")
        return runtime_path

    def run_simulation(self) -> None:
        if self.process.state() != QProcess.ProcessState.NotRunning:
            QMessageBox.information(
                self, "Simulation already running", "Stop the active process before launching another run."
            )
            return

        runtime_path = self.materialize_runtime_case()
        if runtime_path is None:
            return

        if self.current_runtime_case and self.load_existing_checkbox.isChecked():
            data_enabled = bool(nested_get(self.current_runtime_case, "outputs", "data", "enabled", default=True))
            if not data_enabled:
                QMessageBox.critical(self, "Invalid run mode", "load_existing requires outputs.data.enabled to be true.")
                return

        simulator_path = DEFAULT_SIMULATOR_PATH.resolve()
        if not simulator_path.exists():
            QMessageBox.critical(self, "Missing simulator", f"Simulator executable not found: {simulator_path}")
            return

        command_parts = [shlex.quote(str(simulator_path)), shlex.quote(str(runtime_path.resolve()))]
        if self.load_existing_checkbox.isChecked():
            command_parts.append("1")

        command = f"source {shlex.quote(str(DEFAULT_SETUP_SCRIPT))} && {' '.join(command_parts)}"
        self.process_output.appendPlainText(f"\n$ {command}\n")
        self.process.setWorkingDirectory(str(REPO_ROOT))
        self.process.start("/bin/bash", ["-lc", command])
        if self.output_tabs is not None:
            self.output_tabs.setCurrentWidget(self.process_output)

    def stop_simulation(self) -> None:
        if self.process.state() == QProcess.ProcessState.NotRunning:
            return
        self.process.terminate()
        if not self.process.waitForFinished(2000):
            self.process.kill()

    def on_process_started(self) -> None:
        self.run_button.setEnabled(False)
        self.materialize_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.statusBar().showMessage("Simulation running")

    def on_process_finished(self, exit_code: int, _status: QProcess.ExitStatus) -> None:
        self.run_button.setEnabled(True)
        self.materialize_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.append_log(f"Simulation finished with exit code {exit_code}")
        self.refresh_vtk_files()
        self.statusBar().showMessage(f"Simulation finished with exit code {exit_code}")

    def append_process_output(self) -> None:
        data = bytes(self.process.readAllStandardOutput()).decode("utf-8", errors="replace")
        if not data:
            return
        cursor = self.process_output.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        self.process_output.setTextCursor(cursor)
        self.process_output.insertPlainText(data)
        self.process_output.ensureCursorVisible()

    def append_log(self, message: str) -> None:
        self.process_output.appendPlainText(message)