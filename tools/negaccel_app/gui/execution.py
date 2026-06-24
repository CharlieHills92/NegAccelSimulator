"""Process execution helpers for the NegAccel GUI."""

from __future__ import annotations

import json
import re
import shlex
import shutil
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
from negaccel_app.workflow.post_processing import write_runtime_setup_summary


_EXPLICIT_ITERATION_STATUS_PATTERN = re.compile(r"^Iteration\s+(\d+)\s*/\s*(\d+)\s*:\s*(.+?)\s*$")
_LEGACY_ITERATION_START_PATTERN = re.compile(r"Start iteration #\s*(\d+)\s*/\s*(\d+)", re.IGNORECASE)


class ExecutionMixin:
    def _prompt_existing_case_folder_action(self, case_folder: Path) -> str | None:
        message_box = QMessageBox(self)
        message_box.setIcon(QMessageBox.Icon.Warning)
        message_box.setWindowTitle("Simulation folder already exists")
        message_box.setText(
            "A folder with the corresponding case name was already found. What would you like to do?"
        )
        message_box.setInformativeText(str(case_folder))
        delete_button = message_box.addButton("Delete and recreate the folder", QMessageBox.ButtonRole.AcceptRole)
        overwrite_button = message_box.addButton("Overwrite the folder", QMessageBox.ButtonRole.DestructiveRole)
        stop_button = message_box.addButton("Stop the simulation", QMessageBox.ButtonRole.RejectRole)
        message_box.setDefaultButton(stop_button)
        message_box.exec()

        clicked = message_box.clickedButton()
        if clicked == delete_button:
            return "delete"
        if clicked == overwrite_button:
            return "overwrite"
        return None

    def _restore_authoring_file_if_needed(self, case_folder: Path) -> None:
        if self.authoring_path is None:
            return

        try:
            self.authoring_path.relative_to(case_folder)
        except ValueError:
            return

        payload = self.build_authoring_spec()
        write_json(self.authoring_path, payload)

    def _prepare_case_folder_for_run(self, runtime_path: Path) -> str | None:
        case_folder = runtime_path.parent
        if not case_folder.exists():
            return "unchanged"

        action = self._prompt_existing_case_folder_action(case_folder)
        if action == "overwrite":
            self.append_log(f"Using existing case folder without cleanup: {case_folder}")
            return "unchanged"
        if action is None:
            self.append_log(f"Simulation cancelled because the case folder already exists: {case_folder}")
            return None

        try:
            shutil.rmtree(case_folder)
            case_folder.mkdir(parents=True, exist_ok=True)
            self._restore_authoring_file_if_needed(case_folder)
        except Exception as exc:  # pragma: no cover - UI error path
            QMessageBox.critical(self, "Folder preparation failed", str(exc))
            return None

        self.append_log(f"Deleted and recreated case folder: {case_folder}")
        return "recreated"

    def _write_materialized_runtime_case(
        self,
        authoring_spec: dict[str, object],
        runtime_case: dict[str, object],
        runtime_path: Path,
    ) -> Path:
        write_json(runtime_path, runtime_case)
        setup_summary_path = write_runtime_setup_summary(runtime_case, runtime_path)

        self.current_runtime_case = runtime_case
        self.current_runtime_path = runtime_path
        self.authoring_preview.setPlainText(json.dumps(authoring_spec, indent=2))
        self.runtime_preview.setPlainText(json.dumps(runtime_case, indent=2))
        self.append_log(f"Materialized runtime case: {runtime_path}")
        if setup_summary_path is not None:
            self.append_log(f"Wrote setup summary: {setup_summary_path}")
        self.refresh_vtk_files()
        self.statusBar().showMessage(f"Materialized runtime case to {runtime_path}")
        return runtime_path

    def materialize_runtime_case(self) -> Path | None:
        try:
            authoring_spec = self.build_authoring_spec()
            runtime_case, runtime_path = self.build_runtime_case_preview(authoring_spec)
            return self._write_materialized_runtime_case(authoring_spec, runtime_case, runtime_path)
        except Exception as exc:  # pragma: no cover - UI error path
            QMessageBox.critical(self, "Materialization failed", str(exc))
            return None

    def run_simulation(self) -> None:
        if self.process.state() != QProcess.ProcessState.NotRunning:
            QMessageBox.information(
                self, "Simulation already running", "Stop the active process before launching another run."
            )
            return

        try:
            preview_authoring_spec = self.build_authoring_spec()
            preview_runtime_case, runtime_path = self.build_runtime_case_preview(preview_authoring_spec)
        except Exception as exc:  # pragma: no cover - UI error path
            QMessageBox.critical(self, "Materialization failed", str(exc))
            return

        simulator_path = DEFAULT_SIMULATOR_PATH.resolve()
        if not simulator_path.exists():
            QMessageBox.critical(self, "Missing simulator", f"Simulator executable not found: {simulator_path}")
            return

        try:
            program, arguments, display_command = self._build_launch_invocation(runtime_path, simulator_path)
        except ValueError as exc:
            QMessageBox.critical(self, "Invalid server command", str(exc))
            return

        folder_action = self._prepare_case_folder_for_run(runtime_path)
        if folder_action is None:
            return

        try:
            self.auto_save_authoring_file()
        except Exception as exc:  # pragma: no cover - UI error path
            QMessageBox.critical(self, "Save failed", str(exc))
            return

        try:
            authoring_spec = self.build_authoring_spec()
            runtime_case, runtime_path = self.build_runtime_case_preview(authoring_spec)
            self._write_materialized_runtime_case(authoring_spec, runtime_case, runtime_path)
        except Exception as exc:  # pragma: no cover - UI error path
            QMessageBox.critical(self, "Materialization failed", str(exc))
            return

        self.process_output.appendPlainText(f"\n$ {display_command}\n")
        self.process.setWorkingDirectory(str(REPO_ROOT))
        self._stop_requested = False
        self.process.start(program, arguments)
        if self.output_tabs is not None:
            self.output_tabs.setCurrentWidget(self.process_output)

    def stop_simulation(self) -> None:
        if self.process.state() == QProcess.ProcessState.NotRunning:
            return
        self._stop_requested = True
        self._set_simulation_status("Stopping simulation...", None)
        self.process.terminate()
        if not self.process.waitForFinished(2000):
            self.process.kill()

    def on_process_started(self) -> None:
        self.run_button.setEnabled(False)
        self.materialize_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self._process_output_buffer = ""
        self._simulation_current_iteration = 0
        self._simulation_total_iterations = self._simulation_iteration_total()
        self._set_simulation_status(
            "Initialization: mesh generation...",
            self._simulation_progress_value("mesh"),
        )

    def on_process_finished(self, exit_code: int, exit_status: QProcess.ExitStatus) -> None:
        self._flush_process_output_buffer()
        self.run_button.setEnabled(True)
        self.materialize_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.append_log(f"Simulation finished with exit code {exit_code}")
        self.refresh_vtk_files()

        if self._stop_requested:
            status_text = "Run stopped"
            self._set_simulation_status(status_text, None)
        elif exit_status == QProcess.ExitStatus.NormalExit and exit_code == 0:
            status_text = "Run completed!"
            self._set_simulation_status(status_text, self._simulation_progress_value("completed"))
        else:
            status_text = f"Run failed (exit code {exit_code})"
            self._set_simulation_status(status_text, None)

        self._stop_requested = False
        self.statusBar().showMessage(status_text)

    def append_process_output(self) -> None:
        data = bytes(self.process.readAllStandardOutput()).decode("utf-8", errors="replace")
        if not data:
            return
        cursor = self.process_output.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        self.process_output.setTextCursor(cursor)
        self.process_output.insertPlainText(data)
        self.process_output.ensureCursorVisible()
        self._consume_process_output(data)

    def append_log(self, message: str) -> None:
        self.process_output.appendPlainText(message)

    def _simulation_iteration_total(self) -> int:
        total = getattr(self, "_simulation_total_iterations", 0)
        if total:
            return int(total)

        runtime_iterations = None
        if isinstance(self.current_runtime_case, dict):
            runtime_iterations = nested_get(self.current_runtime_case, "simulation", "iterations", default=None)
        if runtime_iterations is None:
            iterations_widget = self.widgets.get("run.iterations") if isinstance(self.widgets, dict) else None
            if iterations_widget is not None:
                runtime_iterations = iterations_widget.value()

        try:
            return max(int(runtime_iterations), 0)
        except (TypeError, ValueError):
            return 0

    def _simulation_progress_max(self) -> int:
        total_iterations = self._simulation_iteration_total()
        if total_iterations > 0:
            return (2 * total_iterations) + 4
        return 4

    def _simulation_progress_value(self, stage: str, iteration: int | None = None) -> int:
        total_steps = self._simulation_progress_max()
        if stage == "mesh":
            return 1
        if stage == "boundary":
            return 2
        if stage == "solving":
            current_iteration = max(int(iteration or 1), 1)
            return min((2 * current_iteration) + 1, total_steps)
        if stage == "tracing":
            current_iteration = max(int(iteration or 1), 1)
            return min((2 * current_iteration) + 2, total_steps)
        if stage == "saving":
            return max(total_steps - 1, 0)
        if stage == "completed":
            return total_steps
        return 0

    def _set_simulation_status(self, text: str, value: int | None) -> None:
        total_steps = self._simulation_progress_max()
        self.simulation_status_bar.setRange(0, total_steps)
        if value is not None:
            self.simulation_status_bar.setValue(max(0, min(int(value), total_steps)))
        self.simulation_status_bar.setFormat(text)
        self.simulation_status_bar.setToolTip(text)
        self.statusBar().showMessage(text)

    def _consume_process_output(self, data: str) -> None:
        self._process_output_buffer += data
        lines = self._process_output_buffer.splitlines(keepends=True)
        remainder = ""
        if lines and not lines[-1].endswith(("\n", "\r")):
            remainder = lines.pop()

        for raw_line in lines:
            self._update_simulation_status_from_output_line(raw_line)

        self._process_output_buffer = remainder

    def _flush_process_output_buffer(self) -> None:
        if not self._process_output_buffer:
            return
        self._update_simulation_status_from_output_line(self._process_output_buffer)
        self._process_output_buffer = ""

    def _update_simulation_status_from_output_line(self, line: str) -> None:
        stripped = line.strip()
        if not stripped:
            return

        if stripped == "Initialization: mesh generation...":
            self._set_simulation_status(stripped, self._simulation_progress_value("mesh"))
            return

        if stripped == "Initialization: boundary definition...":
            self._set_simulation_status(stripped, self._simulation_progress_value("boundary"))
            return

        explicit_iteration = _EXPLICIT_ITERATION_STATUS_PATTERN.match(stripped)
        if explicit_iteration is not None:
            current_iteration = int(explicit_iteration.group(1))
            total_iterations = int(explicit_iteration.group(2))
            phase_text = explicit_iteration.group(3).strip().lower()
            self._simulation_current_iteration = current_iteration
            self._simulation_total_iterations = total_iterations

            if phase_text.startswith("solving poisson"):
                self._set_simulation_status(
                    f"Iteration {current_iteration}/{total_iterations}: solving Poisson...",
                    self._simulation_progress_value("solving", current_iteration),
                )
                return

            if phase_text.startswith("tracing particles"):
                self._set_simulation_status(
                    f"Iteration {current_iteration}/{total_iterations}: tracing particles...",
                    self._simulation_progress_value("tracing", current_iteration),
                )
                return

            self._set_simulation_status(stripped, None)
            return

        legacy_iteration = _LEGACY_ITERATION_START_PATTERN.search(stripped)
        if legacy_iteration is not None:
            current_iteration = int(legacy_iteration.group(1))
            total_iterations = int(legacy_iteration.group(2))
            self._simulation_current_iteration = current_iteration
            self._simulation_total_iterations = total_iterations
            self._set_simulation_status(
                f"Iteration {current_iteration}/{total_iterations}: solving Poisson...",
                self._simulation_progress_value("solving", current_iteration),
            )
            return

        if stripped.startswith("Starting particle trajectory calculation"):
            current_iteration = max(int(getattr(self, "_simulation_current_iteration", 0) or 1), 1)
            total_iterations = max(int(getattr(self, "_simulation_total_iterations", 0) or self._simulation_iteration_total()), 1)
            self._simulation_current_iteration = current_iteration
            self._simulation_total_iterations = total_iterations
            self._set_simulation_status(
                f"Iteration {current_iteration}/{total_iterations}: tracing particles...",
                self._simulation_progress_value("tracing", current_iteration),
            )
            return

        if stripped == "Simulation end: saving outputs...":
            self._set_simulation_status(stripped, self._simulation_progress_value("saving"))
            return

        if stripped.endswith(" saved") or stripped.startswith("Data output disabled by outputs.data settings"):
            self._set_simulation_status(
                "Simulation end: saving outputs...",
                self._simulation_progress_value("saving"),
            )
            return

        if stripped == "Simulation completed!":
            self._set_simulation_status("Run completed!", self._simulation_progress_value("completed"))

    def _current_server_command(self) -> str:
        if not isinstance(self.current_runtime_case, dict):
            return ""
        return str(nested_get(self.current_runtime_case, "execution", "serverCommand", default="")).strip()

    def _build_simulator_payload(self, runtime_path: Path, simulator_path: Path) -> str:
        command_parts = [shlex.quote(str(simulator_path)), shlex.quote(str(runtime_path.resolve()))]
        return f"source {shlex.quote(str(DEFAULT_SETUP_SCRIPT))} && {' '.join(command_parts)}"

    def _build_launch_invocation(self, runtime_path: Path, simulator_path: Path) -> tuple[str, list[str], str]:
        payload = self._build_simulator_payload(runtime_path, simulator_path)
        server_command = self._current_server_command()
        if not server_command:
            return "/bin/bash", ["-lc", payload], f"/bin/bash -lc {shlex.quote(payload)}"

        wrapper_parts = shlex.split(server_command)
        if not wrapper_parts:
            raise ValueError("Command to run on server is empty after parsing")

        program = wrapper_parts[0]
        arguments = [*wrapper_parts[1:], "-lc", payload]
        display_command = " ".join(shlex.quote(part) for part in [program, *arguments])
        return program, arguments, display_command