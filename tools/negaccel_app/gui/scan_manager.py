"""Scan Manager window for parameter scan orchestration and results visualization."""

from __future__ import annotations

import copy
import csv
from datetime import datetime
from pathlib import Path
from typing import Any

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import Qt, QTimer, QProcess
from PySide6.QtWidgets import (
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QLabel,
    QSpinBox,
    QComboBox,
    QTabWidget,
    QProgressBar,
    QTextEdit,
    QFileDialog,
    QMessageBox,
    QHeaderView,
)

from .common import REPO_ROOT, DEFAULT_SIMULATOR_PATH, WorkflowError, authored_to_runtime_case, ensure_case_metadata, load_json, nested_get, write_json
from ..workflow.common import set_path
from ..workflow.post_processing import aggregate_case_diagnostics
from .parameter_flagging import ScanProjectData, ScanCaseDefinition
from .visualization import configure_matplotlib_canvas, configure_matplotlib_toolbar


class ScanManagerWindow(QMainWindow):
    """Independent window for managing parameter scans."""
    
    def __init__(self, parent_window, authoring_path: Path | None = None, scan_project_path: Path | None = None):
        super().__init__()
        self.parent_window = parent_window
        self.authoring_path = authoring_path
        self.scan_project: ScanProjectData | None = None
        self.scan_project_path: Path | None = None
        self._running_processes: dict[int, QProcess] = {}
        self._process_buffer: dict[int, str] = {}
        self._pending_case_indices: list[int] = []
        self._worker_limit = 1
        self._completed_case_count = 0
        
        self.setWindowTitle("NegAccel Scan Manager")
        self.resize(1400, 900)
        
        self._build_ui()
        
        if scan_project_path and scan_project_path.exists():
            self.load_scan_project(scan_project_path)

    def seed_from_main_window(self) -> None:
        """Create an initial scan project from the currently flagged parameters in the main window."""
        project = self.parent_window.build_scan_project_seed()
        self.scan_project = project
        self.scan_project_path = self._default_scan_project_path()
        self.setWindowTitle(f"NegAccel Scan Manager - {project.scanProjectTag}")
        self._sync_table_from_project()
        if not project.parameters:
            self.status_label.setText("No parameters flagged yet. Tick 'P' next to numeric fields in the main window.")
        else:
            self.status_label.setText(f"Seeded {len(project.parameters)} parameter(s) from the main window")
    
    def _build_ui(self) -> None:
        """Build the Scan Manager UI."""
        central = QWidget()
        root_layout = QVBoxLayout(central)
        
        # Control bar
        control_bar = self._build_control_bar()
        root_layout.addWidget(control_bar)
        
        # Tab widget for Parameters and Results
        tabs = QTabWidget()
        self.parameters_tab = self._build_parameters_tab()
        self.results_tab = self._build_results_tab()
        
        tabs.addTab(self.parameters_tab, "Parameters")
        tabs.addTab(self.results_tab, "Results")
        root_layout.addWidget(tabs, 1)
        
        # Status/output bar
        status_bar = self._build_status_bar()
        root_layout.addWidget(status_bar)
        
        self.setCentralWidget(central)
    
    def _build_control_bar(self) -> QWidget:
        """Build the control bar with execution options."""
        bar = QWidget()
        layout = QHBoxLayout(bar)
        
        self.load_project_button = QPushButton("Load Project")
        self.load_project_button.clicked.connect(self.choose_scan_project)
        layout.addWidget(self.load_project_button)
        
        self.save_project_button = QPushButton("Save Project")
        self.save_project_button.clicked.connect(self.save_current_project)
        layout.addWidget(self.save_project_button)
        
        layout.addSpacing(20)
        
        # Execution mode and worker count
        layout.addWidget(QLabel("Execution:"))
        self.execution_mode_combo = QComboBox()
        self.execution_mode_combo.addItems(["Sequential", "Parallel"])
        layout.addWidget(self.execution_mode_combo)
        
        layout.addWidget(QLabel("Workers:"))
        self.worker_count_spin = QSpinBox()
        self.worker_count_spin.setRange(1, 16)
        self.worker_count_spin.setValue(4)
        layout.addWidget(self.worker_count_spin)
        
        layout.addSpacing(20)
        
        # Run and Stop buttons
        self.run_scan_button = QPushButton("Run Scan")
        self.run_scan_button.clicked.connect(self.run_scan)
        layout.addWidget(self.run_scan_button)
        
        self.stop_scan_button = QPushButton("Stop")
        self.stop_scan_button.setEnabled(False)
        self.stop_scan_button.clicked.connect(self.stop_scan)
        layout.addWidget(self.stop_scan_button)
        
        layout.addStretch()
        
        return bar
    
    def _build_parameters_tab(self) -> QWidget:
        """Build the parameters tab with table."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Parameter table
        self.parameters_table = QTableWidget()
        self.parameters_table.setColumnCount(0)
        self.parameters_table.setRowCount(0)
        self.parameters_table.horizontalHeader().setStretchLastSection(False)
        layout.addWidget(self.parameters_table, 1)
        
        # Table operations
        table_ops_layout = QHBoxLayout()
        self.add_row_button = QPushButton("Add Row")
        self.add_row_button.clicked.connect(self.add_table_row)
        table_ops_layout.addWidget(self.add_row_button)
        
        self.duplicate_row_button = QPushButton("Duplicate Selected")
        self.duplicate_row_button.clicked.connect(self.duplicate_table_row)
        table_ops_layout.addWidget(self.duplicate_row_button)
        
        self.delete_row_button = QPushButton("Delete Selected")
        self.delete_row_button.clicked.connect(self.delete_table_row)
        table_ops_layout.addWidget(self.delete_row_button)
        
        self.load_case_button = QPushButton("Load Case from JSON")
        self.load_case_button.clicked.connect(self.load_case_from_json)
        table_ops_layout.addWidget(self.load_case_button)
        
        table_ops_layout.addStretch()
        layout.addLayout(table_ops_layout)
        
        return tab
    
    def _build_results_tab(self) -> QWidget:
        """Build the results tab for visualization."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Placeholder for now; will be extended with matplotlib plotting
        layout.addWidget(QLabel("Results visualization will appear here after scan completion."))
        
        # Metrics selection (will be expanded)
        metrics_layout = QHBoxLayout()
        metrics_layout.addWidget(QLabel("X-Metric:"))
        self.x_metric_combo = QComboBox()
        metrics_layout.addWidget(self.x_metric_combo)
        
        metrics_layout.addWidget(QLabel("Y-Metric:"))
        self.y_metric_combo = QComboBox()
        metrics_layout.addWidget(self.y_metric_combo)
        
        self.plot_button = QPushButton("Plot")
        self.plot_button.clicked.connect(self.generate_plot)
        metrics_layout.addWidget(self.plot_button)
        
        self.export_csv_button = QPushButton("Export CSV")
        self.export_csv_button.clicked.connect(self.export_results_csv)
        metrics_layout.addWidget(self.export_csv_button)
        
        metrics_layout.addStretch()
        layout.addLayout(metrics_layout)
        
        self.results_figure = Figure(figsize=(7.2, 4.8), tight_layout=True)
        self.results_canvas = FigureCanvasQTAgg(self.results_figure)
        configure_matplotlib_canvas(self.results_canvas)
        self.results_toolbar = NavigationToolbar2QT(self.results_canvas, tab)
        configure_matplotlib_toolbar(self.results_toolbar)
        layout.addWidget(self.results_toolbar)
        layout.addWidget(self.results_canvas, 1)

        self.results_display = QTextEdit()
        self.results_display.setReadOnly(True)
        self.results_display.setMaximumHeight(180)
        layout.addWidget(self.results_display)

        self._plot_placeholder("Run or load scan cases to plot aggregated metrics.")
        
        return tab
    
    def _build_status_bar(self) -> QWidget:
        """Build the status bar with progress."""
        bar = QWidget()
        layout = QHBoxLayout(bar)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        layout.addWidget(self.progress_bar)
        
        self.status_label = QLabel("Ready")
        layout.addWidget(self.status_label)
        
        return bar

    def _default_scan_project_path(self) -> Path:
        if self.authoring_path is not None:
            return self.authoring_path.parent / f"negaccel-scan-project-{self.authoring_path.stem}.json"
        tag = self.scan_project.scanProjectTag if self.scan_project is not None else "untitled"
        return REPO_ROOT / f"negaccel-scan-project-{tag}.json"

    def _scan_output_root(self) -> Path:
        if self.scan_project_path is not None:
            root = self.scan_project_path.with_suffix("")
            root.mkdir(parents=True, exist_ok=True)
            return root
        root = self._default_scan_project_path().with_suffix("")
        root.mkdir(parents=True, exist_ok=True)
        return root

    def _plot_placeholder(self, message: str) -> None:
        self.results_figure.clear()
        axis = self.results_figure.add_subplot(111)
        axis.axis("off")
        axis.text(0.5, 0.5, message, ha="center", va="center", fontsize=11)
        self.results_canvas.draw_idle()

    def _collect_case_metric_rows(self) -> list[dict[str, Any]]:
        if self.scan_project is None:
            return []

        rows: list[dict[str, Any]] = []
        for case in self.scan_project.cases:
            row: dict[str, Any] = {
                "caseIndex": case.caseIndex,
                "caseLabel": case.caseLabel or self._case_file_stem(case.caseIndex),
                "status": case.status,
            }
            for parameter, value in zip(self.scan_project.parameters, case.parameterValues):
                row[parameter.label] = value

            if case.status == "completed" and case.runtimeJsonPath:
                try:
                    aggregated = aggregate_case_diagnostics(Path(case.runtimeJsonPath))
                    diagnostic_rows = aggregated.get("diagnosticSummary", {}).get("rows", [])
                    if diagnostic_rows:
                        final_row = diagnostic_rows[-1]
                        for key, value in final_row.items():
                            row[f"diagnostic:{key}"] = value
                    negative_rows = (aggregated.get("negativeIonSummary") or {}).get("rows", [])
                    if negative_rows:
                        final_negative_row = negative_rows[-1]
                        for key, value in final_negative_row.items():
                            row[f"negativeIon:{key}"] = value
                    for species_tag, grid_summary in aggregated.get("gridPowerSummaries", {}).items():
                        total_power = grid_summary.get("totalIncludedBeamPowerWatts")
                        if total_power is not None:
                            row[f"gridPower:{species_tag}:totalIncludedBeamPowerWatts"] = total_power
                except Exception as exc:
                    row["aggregationError"] = str(exc)

            rows.append(row)
        return rows

    def _refresh_metric_selectors(self) -> None:
        metric_rows = self._collect_case_metric_rows()
        metric_names: list[str] = []
        for row in metric_rows:
            for key, value in row.items():
                if key in {"caseIndex", "caseLabel", "status", "aggregationError"}:
                    continue
                if isinstance(value, (int, float)) and key not in metric_names:
                    metric_names.append(key)

        current_x = self.x_metric_combo.currentText()
        current_y = self.y_metric_combo.currentText()
        self.x_metric_combo.blockSignals(True)
        self.y_metric_combo.blockSignals(True)
        self.x_metric_combo.clear()
        self.y_metric_combo.clear()
        self.x_metric_combo.addItems(metric_names)
        self.y_metric_combo.addItems(metric_names)
        if current_x:
            index = self.x_metric_combo.findText(current_x)
            if index >= 0:
                self.x_metric_combo.setCurrentIndex(index)
        if current_y:
            index = self.y_metric_combo.findText(current_y)
            if index >= 0:
                self.y_metric_combo.setCurrentIndex(index)
        if not current_x and metric_names:
            self.x_metric_combo.setCurrentIndex(0)
        if not current_y and len(metric_names) > 1:
            self.y_metric_combo.setCurrentIndex(1)
        elif not current_y and metric_names:
            self.y_metric_combo.setCurrentIndex(0)
        self.x_metric_combo.blockSignals(False)
        self.y_metric_combo.blockSignals(False)
    
    def add_table_row(self) -> None:
        """Add a new empty row to the parameters table."""
        if not self.scan_project:
            QMessageBox.warning(self, "No Project", "Load or create a scan project first.")
            return
        
        row_index = len(self.scan_project.cases)
        new_case = ScanCaseDefinition(
            caseIndex=row_index,
            parameterValues=[None] * len(self.scan_project.parameters)
        )
        self.scan_project.cases.append(new_case)
        self._sync_table_from_project()
    
    def duplicate_table_row(self) -> None:
        """Duplicate the currently selected row."""
        if not self.scan_project:
            return
        
        current_row = self.parameters_table.currentRow()
        if current_row < 0 or current_row >= len(self.scan_project.cases):
            QMessageBox.warning(self, "No Selection", "Select a row to duplicate.")
            return
        
        source_case = self.scan_project.cases[current_row]
        new_index = len(self.scan_project.cases)
        new_case = ScanCaseDefinition(
            caseIndex=new_index,
            parameterValues=source_case.parameterValues.copy(),
            caseLabel=f"{source_case.caseLabel} (copy)" if source_case.caseLabel else ""
        )
        self.scan_project.cases.append(new_case)
        self._sync_table_from_project()
    
    def delete_table_row(self) -> None:
        """Delete the currently selected row."""
        if not self.scan_project:
            return
        
        current_row = self.parameters_table.currentRow()
        if current_row < 0 or current_row >= len(self.scan_project.cases):
            QMessageBox.warning(self, "No Selection", "Select a row to delete.")
            return
        
        self.scan_project.cases.pop(current_row)
        # Renumber remaining cases
        for i, case in enumerate(self.scan_project.cases):
            case.caseIndex = i
        self._sync_table_from_project()
    
    def load_case_from_json(self) -> None:
        """Load a case from an existing case_NNN.json file."""
        selected, _ = QFileDialog.getOpenFileName(
            self,
            "Select case JSON file",
            str(Path.home()),
            "Case JSON (case_*.json)",
        )
        if not selected:
            return
        
        try:
            case_data = load_json(Path(selected))
            # Extract parameter values based on current parameter definitions
            if not self.scan_project:
                QMessageBox.warning(self, "No Project", "Load a scan project first.")
                return
            
            param_values = []
            for param_def in self.scan_project.parameters:
                # Parse the path and get the value from case data
                path_parts = param_def.path.replace("[", ".").replace("]", "").split(".")
                value = nested_get(case_data, *path_parts)
                param_values.append(value)
            
            new_index = len(self.scan_project.cases)
            new_case = ScanCaseDefinition(
                caseIndex=new_index,
                parameterValues=param_values,
                caseLabel=f"Loaded from {Path(selected).name}"
            )
            self.scan_project.cases.append(new_case)
            self._sync_table_from_project()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load case: {e}")
    
    def _sync_table_from_project(self) -> None:
        """Sync the GUI table from the current project data."""
        if not self.scan_project:
            return
        
        # Clear and rebuild table columns
        self.parameters_table.clear()
        self.parameters_table.setColumnCount(len(self.scan_project.parameters) + 1)
        
        headers = [p.label for p in self.scan_project.parameters] + ["Status"]
        self.parameters_table.setHorizontalHeaderLabels(headers)
        
        # Populate rows
        self.parameters_table.setRowCount(len(self.scan_project.cases))
        for row_idx, case in enumerate(self.scan_project.cases):
            # Parameter columns
            for col_idx, param_val in enumerate(case.parameterValues):
                item = QTableWidgetItem(str(param_val) if param_val is not None else "")
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEditable)
                self.parameters_table.setItem(row_idx, col_idx, item)
            
            # Status column
            status_item = QTableWidgetItem(case.status)
            status_item.setFlags(status_item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.parameters_table.setItem(row_idx, len(self.scan_project.parameters), status_item)
        
        # Auto-resize columns
        self.parameters_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    
    def _sync_project_from_table(self) -> None:
        """Sync project data from the current GUI table state."""
        if not self.scan_project:
            return
        
        for row_idx in range(self.parameters_table.rowCount()):
            if row_idx < len(self.scan_project.cases):
                case = self.scan_project.cases[row_idx]
                # Update parameter values from table
                for col_idx in range(len(self.scan_project.parameters)):
                    item = self.parameters_table.item(row_idx, col_idx)
                    if item:
                        try:
                            param_def = self.scan_project.parameters[col_idx]
                            if param_def.parameterType == "number":
                                case.parameterValues[col_idx] = float(item.text())
                            elif param_def.parameterType == "integer":
                                case.parameterValues[col_idx] = int(item.text())
                            elif param_def.parameterType == "boolean":
                                case.parameterValues[col_idx] = item.text().lower() in ("true", "1", "yes")
                            else:
                                case.parameterValues[col_idx] = item.text()
                        except (ValueError, IndexError):
                            pass
    
    def choose_scan_project(self) -> None:
        """Choose a scan project file to load."""
        selected, _ = QFileDialog.getOpenFileName(
            self,
            "Select scan project JSON",
            str(Path.home()),
            "Scan Project (negaccel-scan-project-*.json)",
        )
        if selected:
            self.load_scan_project(Path(selected))
    
    def load_scan_project(self, project_path: Path) -> None:
        """Load a scan project from file."""
        try:
            self.scan_project = ScanProjectData.from_file(project_path)
            self.scan_project_path = project_path
            if self.authoring_path is None and self.scan_project.authoringCasePath:
                self.authoring_path = Path(self.scan_project.authoringCasePath)
            self.setWindowTitle(f"NegAccel Scan Manager - {self.scan_project.scanProjectTag}")
            self._sync_table_from_project()
            self._refresh_metric_selectors()
            self.status_label.setText(f"Loaded: {project_path.name}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load project: {e}")
    
    def save_current_project(self) -> None:
        """Save the current project to file."""
        if not self.scan_project:
            QMessageBox.warning(self, "No Project", "No project loaded.")
            return
        
        self._sync_project_from_table()
        
        if not self.scan_project_path:
            # Ask user where to save
            tag = self.scan_project.scanProjectTag or "untitled"
            suggested_name = f"negaccel-scan-project-{tag}.json"
            selected, _ = QFileDialog.getSaveFileName(
                self,
                "Save scan project",
                str(Path.home() / suggested_name),
                "JSON files (*.json)",
            )
            if not selected:
                return
            self.scan_project_path = Path(selected)
        
        try:
            self.scan_project.save(self.scan_project_path)
            self._refresh_metric_selectors()
            self.status_label.setText(f"Saved: {self.scan_project_path.name}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save project: {e}")

    def _base_authoring_spec(self) -> dict[str, Any]:
        authoring_spec = self.parent_window.build_authoring_spec()
        geometry_workspace = self.parent_window.widgets.get("geometry.solidsEditor")
        current_geometry_document = getattr(geometry_workspace, "current_geometry_document", None)
        if callable(current_geometry_document):
            authoring_spec["geometry"] = current_geometry_document()
        if self.parent_window.authoring_path is not None:
            authoring_spec["__source_path"] = self.parent_window.authoring_path
        return authoring_spec

    def _case_file_stem(self, case_index: int) -> str:
        return f"case_{case_index + 1:03d}"

    def _materialize_scan_cases(self) -> None:
        if self.scan_project is None:
            raise WorkflowError("No scan project is loaded")

        base_authoring_spec = self._base_authoring_spec()
        base_case_tag = self.parent_window.widgets["caseTag"].text().strip() or "negaccel_scan"
        output_root = self._scan_output_root()

        for case in self.scan_project.cases:
            authoring_spec = copy.deepcopy(base_authoring_spec)
            for parameter, value in zip(self.scan_project.parameters, case.parameterValues):
                set_path(authoring_spec, parameter.path, value)

            case_stem = self._case_file_stem(case.caseIndex)
            case_tag = f"{base_case_tag}_{case_stem}"
            metadata = authoring_spec.setdefault("metadata", {})
            metadata["caseTag"] = case_tag

            runtime_case = authored_to_runtime_case(authoring_spec)
            case_directory = output_root / case_stem
            case_directory.mkdir(parents=True, exist_ok=True)
            runtime_path = case_directory / f"{case_stem}.json"
            ensure_case_metadata(runtime_case, case_tag, case_directory)
            write_json(runtime_path, runtime_case)

            case.runtimeJsonPath = str(runtime_path)
            case.outputDirectoryPath = str(case_directory)
            case.status = "pending"
            case.errorMessage = None
            case.completedAt = None

        self._sync_table_from_project()

    def _update_case_status_cell(self, case_index: int, status: str) -> None:
        if self.scan_project is None or case_index >= len(self.scan_project.cases):
            return
        status_column = len(self.scan_project.parameters)
        item = self.parameters_table.item(case_index, status_column)
        if item is None:
            item = QTableWidgetItem()
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.parameters_table.setItem(case_index, status_column, item)
        item.setText(status)

    def _append_case_output(self, case_index: int, process: QProcess) -> None:
        data = bytes(process.readAllStandardOutput()).decode("utf-8", errors="replace")
        if not data:
            return
        self._process_buffer[case_index] = self._process_buffer.get(case_index, "") + data
        self.results_display.append(f"[case {case_index + 1:03d}]\n{data.rstrip()}\n")

    def _start_case_process(self, case_index: int) -> None:
        if self.scan_project is None:
            return
        case = self.scan_project.cases[case_index]
        runtime_path = Path(case.runtimeJsonPath or "")
        simulator_path = DEFAULT_SIMULATOR_PATH.resolve()
        if not runtime_path.exists():
            raise WorkflowError(f"Runtime JSON not found for case {case_index + 1:03d}: {runtime_path}")
        if not simulator_path.exists():
            raise WorkflowError(f"Simulator executable not found: {simulator_path}")

        program, arguments, display_command = self.parent_window._build_launch_invocation(runtime_path, simulator_path)

        process = QProcess(self)
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        process.setWorkingDirectory(str(REPO_ROOT))
        process.readyReadStandardOutput.connect(lambda idx=case_index, proc=process: self._append_case_output(idx, proc))
        process.finished.connect(lambda exit_code, exit_status, idx=case_index: self._on_case_finished(idx, exit_code, exit_status))

        case.status = "running"
        self._update_case_status_cell(case_index, case.status)
        self._running_processes[case_index] = process
        self.results_display.append(f"$ {display_command}\n")
        process.start(program, arguments)

    def _launch_available_cases(self) -> None:
        while self._pending_case_indices and len(self._running_processes) < self._worker_limit:
            case_index = self._pending_case_indices.pop(0)
            self._start_case_process(case_index)

    def _on_case_finished(self, case_index: int, exit_code: int, exit_status: QProcess.ExitStatus) -> None:
        if self.scan_project is None:
            return
        case = self.scan_project.cases[case_index]
        self._running_processes.pop(case_index, None)
        if exit_status == QProcess.ExitStatus.NormalExit and exit_code == 0:
            case.status = "completed"
            case.errorMessage = None
        else:
            case.status = "failed"
            case.errorMessage = f"exit code {exit_code}"
        case.completedAt = datetime.utcnow().isoformat() + "Z"
        self._completed_case_count += 1
        self._update_case_status_cell(case_index, case.status)

        total = len(self.scan_project.cases)
        progress = int((100 * self._completed_case_count) / total) if total else 0
        self.progress_bar.setValue(progress)
        self.status_label.setText(f"Completed {self._completed_case_count}/{total} cases")

        if self.scan_project_path is not None:
            self.scan_project.save(self.scan_project_path)

        if self._pending_case_indices:
            self._launch_available_cases()
            return

        if not self._running_processes:
            self.run_scan_button.setEnabled(True)
            self.stop_scan_button.setEnabled(False)
            self._refresh_metric_selectors()
            self.status_label.setText(f"Scan finished: {self._completed_case_count}/{total} cases processed")
    
    def run_scan(self) -> None:
        """Execute the scan."""
        if not self.scan_project or not self.scan_project.cases:
            QMessageBox.warning(self, "No Cases", "Add cases to the scan table first.")
            return
        
        # Save current table state before running
        self._sync_project_from_table()
        self.save_current_project()
        
        mode = self.execution_mode_combo.currentText().lower()
        workers = self.worker_count_spin.value() if mode == "parallel" else 1

        try:
            self._materialize_scan_cases()
        except Exception as exc:
            QMessageBox.critical(self, "Case generation failed", str(exc))
            return
        
        self.run_scan_button.setEnabled(False)
        self.stop_scan_button.setEnabled(True)
        self.status_label.setText(f"Starting scan ({mode}, {workers} workers)...")
        self.progress_bar.setValue(0)
        self._worker_limit = max(1, workers)
        self._completed_case_count = 0
        self._pending_case_indices = [case.caseIndex for case in self.scan_project.cases]

        try:
            self._launch_available_cases()
        except Exception as exc:
            self.run_scan_button.setEnabled(True)
            self.stop_scan_button.setEnabled(False)
            QMessageBox.critical(self, "Scan launch failed", str(exc))
    
    def stop_scan(self) -> None:
        """Stop the current scan execution."""
        self._pending_case_indices = []
        for case_index, process in list(self._running_processes.items()):
            if process.state() != QProcess.ProcessState.NotRunning:
                process.terminate()
                if not process.waitForFinished(1000):
                    process.kill()
            if self.scan_project is not None and case_index < len(self.scan_project.cases):
                self.scan_project.cases[case_index].status = "skipped"
                self._update_case_status_cell(case_index, "skipped")
        self._running_processes.clear()
        self.status_label.setText("Scan stopped")
        self.run_scan_button.setEnabled(True)
        self.stop_scan_button.setEnabled(False)
    
    def generate_plot(self) -> None:
        """Generate a plot from results."""
        x_metric = self.x_metric_combo.currentText().strip()
        y_metric = self.y_metric_combo.currentText().strip()
        if not x_metric or not y_metric:
            self._plot_placeholder("No numeric metrics are available yet.")
            return

        metric_rows = self._collect_case_metric_rows()
        plotted_rows = [
            row
            for row in metric_rows
            if isinstance(row.get(x_metric), (int, float)) and isinstance(row.get(y_metric), (int, float))
        ]
        if not plotted_rows:
            self._plot_placeholder("Selected metrics are unavailable for the current completed cases.")
            return

        x_values = [float(row[x_metric]) for row in plotted_rows]
        y_values = [float(row[y_metric]) for row in plotted_rows]
        labels = [str(row.get("caseLabel") or row.get("caseIndex")) for row in plotted_rows]

        self.results_figure.clear()
        axis = self.results_figure.add_subplot(111)
        axis.plot(x_values, y_values, marker="o", linewidth=1.6, color="#c24d2c")
        for x_value, y_value, label in zip(x_values, y_values, labels):
            axis.annotate(label, (x_value, y_value), textcoords="offset points", xytext=(5, 5), fontsize=8)
        axis.set_xlabel(x_metric)
        axis.set_ylabel(y_metric)
        axis.set_title(f"{y_metric} vs {x_metric}")
        axis.grid(True, alpha=0.3)
        self.results_canvas.draw_idle()
    
    def export_results_csv(self) -> None:
        """Export scan results to CSV."""
        metric_rows = self._collect_case_metric_rows()
        if not metric_rows:
            QMessageBox.information(self, "No Data", "No scan data is available to export.")
            return

        selected, _ = QFileDialog.getSaveFileName(
            self,
            "Export scan results CSV",
            str((self._scan_output_root() / "scan_results.csv").resolve()),
            "CSV files (*.csv)",
        )
        if not selected:
            return

        fieldnames: list[str] = []
        for row in metric_rows:
            for key in row.keys():
                if key not in fieldnames:
                    fieldnames.append(key)

        with Path(selected).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(metric_rows)

        self.status_label.setText(f"Exported CSV: {Path(selected).name}")
