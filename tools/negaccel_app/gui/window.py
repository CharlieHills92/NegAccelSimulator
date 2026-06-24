"""Main window composition for the NegAccel GUI."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .authoring import AuthoringMixin
from .common import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFontDatabase,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QProgressBar,
    QPushButton,
    QProcess,
    QSpinBox,
    QSplitter,
    QTabWidget,
    QTextEdit,
    QTimer,
    QVBoxLayout,
    QWidget,
    Qt,
    EXAMPLE_AUTHORING_PATH,
    build_runtime_path,
    create_scrollable_form,
    install_qt_message_filter,
    runtime_path_parent_seed,
)
from .execution import ExecutionMixin
from .parameter_flagging import ParameterFlagRegistry, ScanCaseDefinition, ScanProjectData
from .results import ResultsMixin
from .sections import FORM_SECTIONS
from .scan_manager import ScanManagerWindow


class NegAccelMainWindow(AuthoringMixin, ExecutionMixin, ResultsMixin, QMainWindow):
    def __init__(
        self,
        authoring_path: Path | None = None,
        runtime_path: Path | None = None,
    ) -> None:
        super().__init__()
        install_qt_message_filter()
        self.setWindowTitle("NegAccel Studio")

        self.authoring_path: Path | None = None
        self.loaded_authoring_template = self._load_template(EXAMPLE_AUTHORING_PATH)
        self.current_runtime_case: dict[str, Any] | None = None
        self.current_runtime_path: Path | None = None
        self.preview_updates_blocked = False
        self.default_output_parent_directory = runtime_path_parent_seed(runtime_path)
        self.initial_runtime_case_tag = runtime_path.stem if runtime_path is not None else None
        self.runtime_path_auto = True

        self.preview_timer = QTimer(self)
        self.preview_timer.setSingleShot(True)
        self.preview_timer.timeout.connect(self.refresh_runtime_preview)

        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.append_process_output)
        self.process.started.connect(self.on_process_started)
        self.process.finished.connect(self.on_process_finished)
        self._process_output_buffer = ""
        self._simulation_current_iteration = 0
        self._simulation_total_iterations = 0
        self._stop_requested = False

        self.widgets: dict[str, QWidget] = {}
        self.output_tabs: QTabWidget | None = None
        self.scan_manager_window: ScanManagerWindow | None = None
        self.parameter_flag_registry = ParameterFlagRegistry()
        self.parameter_widgets_by_path: dict[str, QWidget] = {}
        self.parameter_checkboxes_by_path: dict[str, QCheckBox] = {}
        self.parameter_labels_by_path: dict[str, str] = {}

        self._build_ui(runtime_path)
        self._resize_to_available_screen()
        self._apply_style()

        if authoring_path is not None:
            self.load_authoring_from_path(authoring_path)
        else:
            self.populate_from_authoring_spec(self.loaded_authoring_template, EXAMPLE_AUTHORING_PATH)
            if self.initial_runtime_case_tag is not None:
                self.widgets["caseTag"].setText(self.initial_runtime_case_tag)
                self.on_case_tag_edited(self.initial_runtime_case_tag)
                self.refresh_runtime_preview()

    def _resize_to_available_screen(self) -> None:
        screen = self.screen() or QApplication.primaryScreen()
        if screen is None:
            self.resize(1560, 980)
            return

        available = screen.availableGeometry()
        self.setMaximumSize(available.width(), available.height())
        width = min(1560, max(840, available.width() - 40))
        height = min(980, max(640, available.height() - 60))
        width = min(width, available.width())
        height = min(height, available.height())

        self.resize(width, height)

        frame = self.frameGeometry()
        frame.moveCenter(available.center())
        self.move(frame.topLeft())

    def showEvent(self, event) -> None:
        super().showEvent(event)
        QTimer.singleShot(0, self._resize_to_available_screen)

    def _apply_style(self) -> None:
        fixed_font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        self.authoring_preview.setFont(fixed_font)
        self.runtime_preview.setFont(fixed_font)
        self.process_output.setFont(fixed_font)
        QApplication.instance().setStyle("Fusion")
        self.setStyleSheet(
            """
            QWidget {
                background: #f4efe7;
                color: #1f2933;
                font-size: 12px;
            }
            QWidget[matplotlibCanvas="true"] {
                background: transparent;
            }
            QToolBar[matplotlibToolbar="true"] {
                background: transparent;
                border: none;
            }
            QToolBar[matplotlibToolbar="true"] QToolButton {
                background: transparent;
                border: 1px solid transparent;
                border-radius: 4px;
                padding: 4px;
            }
            QToolBar[matplotlibToolbar="true"] QToolButton:hover {
                background: rgba(194, 77, 44, 0.10);
                border-color: rgba(194, 77, 44, 0.30);
            }
            QMainWindow, QTabWidget::pane, QGroupBox, QPlainTextEdit, QTextEdit, QLineEdit, QComboBox,
            QSpinBox, QDoubleSpinBox {
                background: #fffaf2;
            }
            QGroupBox {
                border: 1px solid #d9c7b7;
                border-radius: 8px;
                margin-top: 12px;
                padding-top: 10px;
                font-weight: 600;
            }
            QGroupBox::title {
                left: 10px;
                padding: 0 4px;
            }
            QPushButton {
                background: #c24d2c;
                color: white;
                border: none;
                border-radius: 6px;
                padding: 8px 14px;
                font-weight: 600;
            }
            QPushButton:disabled {
                background: #b8aca1;
                color: #f7f1ea;
            }
            QPushButton:hover:!disabled {
                background: #a63f23;
            }
            QTabBar::tab {
                background: #e8ddd0;
                color: #1f2933;
                border: 1px solid #d9c7b7;
                border-bottom: none;
                padding: 7px 12px;
                margin-right: 1px;
            }
            QTabBar::tab:selected {
                background: #fffaf2;
            }
            QTabBar::tab:!selected {
                margin-top: 2px;
            }
            QTabBar::scroller {
                width: 18px;
            }
            QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QPlainTextEdit, QTextEdit {
                border: 1px solid #ccb7a5;
                border-radius: 6px;
                padding: 5px;
                selection-background-color: #c24d2c;
            }
            """
        )

    def _build_ui(self, runtime_path: Path | None) -> None:
        central = QWidget()
        root_layout = QVBoxLayout(central)

        self.authoring_path_display = QLineEdit()
        self.authoring_path_display.setReadOnly(True)
        self.authoring_path_display.setToolTip(
            "Currently opened authoring JSON file. This is the source document being previewed, materialized, and saved."
        )
        self.widgets["execution.outputParentDirectory"] = QLineEdit(str(self.default_output_parent_directory))
        self.widgets["execution.outputParentDirectory"].setMinimumWidth(280)
        self.widgets["execution.outputParentDirectory"].setToolTip(
            "Parent directory under which the GUI creates the case-tag folder and writes the runtime JSON plus all solver outputs."
        )
        self.widgets["execution.outputParentDirectory"].textChanged.connect(self.on_output_parent_directory_changed)
        self.output_parent_directory_browse_button = QPushButton("Browse")
        self.output_parent_directory_browse_button.clicked.connect(self.choose_output_parent_directory)
        self.runtime_path_edit = QLineEdit(str(runtime_path or build_runtime_path("negaccel_case", self.default_output_parent_directory)))
        self.runtime_path_edit.setReadOnly(True)
        self.runtime_path_edit.setToolTip(
            "Generated output path for the solver-ready runtime JSON. This path is derived from the selected output parent folder and the current case tag."
        )
        self.load_button = QPushButton("Load config")
        self.load_button.clicked.connect(self.choose_authoring_file)
        self.save_button = QPushButton("Save config")
        self.save_button.clicked.connect(self.save_authoring_file)
        self.widgets["caseTag"] = QLineEdit()
        self.widgets["caseTag"].textEdited.connect(self.on_case_tag_edited)
        self.widgets["caseTag"].setMinimumWidth(220)
        self.case_tag_from_json_checkbox = QCheckBox("Get case tag from JSON file")
        self.case_tag_from_json_checkbox.setEnabled(False)
        self.case_tag_from_json_checkbox.toggled.connect(self.on_case_tag_from_json_toggled)
        self.case_tag_inline_label = QLabel("Case tag")
        self.materialize_button = QPushButton("Write Runtime JSON")
        self.materialize_button.clicked.connect(self.materialize_runtime_case)
        self.materialize_button.setToolTip(
            "Generate the solver-ready runtime JSON from the current authoring form without starting the simulation."
        )
        self.run_button = QPushButton("Run Simulation")
        self.run_button.clicked.connect(self.run_simulation)
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop_simulation)
        self.stop_button.setEnabled(False)
        self.simulation_status_bar = QProgressBar()
        self.simulation_status_bar.setRange(0, 1)
        self.simulation_status_bar.setValue(0)
        self.simulation_status_bar.setTextVisible(True)
        self.simulation_status_bar.setMinimumWidth(420)
        self.simulation_status_bar.setFormat("Ready")
        self.simulation_status_bar.setToolTip("Current simulation phase and coarse run progress.")
        self.scan_manager_button = QPushButton("Scan Manager")
        self.scan_manager_button.clicked.connect(self.open_scan_manager)
        self.scan_manager_button.setToolTip("Open the parameter scan manager for batch simulations.")

        general_box = QGroupBox("General")
        general_form = QFormLayout(general_box)
        authoring_output_row = QHBoxLayout()
        authoring_output_row.setContentsMargins(0, 0, 0, 0)
        authoring_output_row.addWidget(QLabel("Authoring JSON"))
        authoring_output_row.addWidget(self.authoring_path_display, 1)
        authoring_output_row.addWidget(self.load_button)
        authoring_output_row.addWidget(self.save_button)
        authoring_output_row.addSpacing(18)
        authoring_output_row.addWidget(QLabel("Output parent folder"))
        authoring_output_row.addWidget(self.widgets["execution.outputParentDirectory"], 1)
        authoring_output_row.addWidget(self.output_parent_directory_browse_button)
        authoring_output_wrapper = QWidget()
        authoring_output_wrapper.setLayout(authoring_output_row)
        general_form.addRow("", authoring_output_wrapper)

        run_button_row = QHBoxLayout()
        run_button_row.addWidget(self.case_tag_inline_label)
        run_button_row.addWidget(self.widgets["caseTag"])
        run_button_row.addWidget(self.case_tag_from_json_checkbox)
        run_button_row.addWidget(self.run_button)
        run_button_row.addWidget(self.stop_button)
        run_button_row.addWidget(self.simulation_status_bar, 1)
        run_button_row.addWidget(self.scan_manager_button)
        run_button_row.addStretch(1)
        run_button_wrapper = QWidget()
        run_button_wrapper.setLayout(run_button_row)
        general_form.addRow("", run_button_wrapper)

        root_layout.addWidget(general_box)

        self.output_tabs = self._build_main_tabs()
        root_layout.addWidget(self.output_tabs, 1)

        self.setCentralWidget(central)
        self.statusBar().showMessage("Ready")

    def _build_main_tabs(self) -> QTabWidget:
        tabs = QTabWidget()
        for section in FORM_SECTIONS:
            existing_widget_keys = set(self.widgets)
            if section.title == "Metadata":
                tabs.addTab(self._build_metadata_tab(section, existing_widget_keys), section.title)
            elif section.build_workspace is not None:
                tabs.addTab(section.build_workspace(self), section.title)
            else:
                form_layout = section.build_form(self)
                self._decorate_new_numeric_rows(form_layout, existing_widget_keys)
                tabs.addTab(create_scrollable_form(form_layout), section.title)
        tabs.addTab(self._build_run_log_tab(), "Run Log")
        tabs.addTab(self._build_visualization_tab(), "Visualization")
        self._connect_preview_updates()
        return tabs

    def _build_metadata_tab(self, section, existing_widget_keys: set[str]) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)

        form_layout = section.build_form(self)
        self._decorate_new_numeric_rows(form_layout, existing_widget_keys)

        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(create_scrollable_form(form_layout))
        splitter.addWidget(self._build_json_preview_panel())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([250, 640])

        layout.addWidget(splitter)
        return page

    def open_scan_manager(self) -> None:
        """Open or bring focus to the Scan Manager window."""
        if self.scan_manager_window is None or not self.scan_manager_window.isVisible():
            # Derive scan project path from authoring path
            scan_project_path = None
            if self.authoring_path is not None:
                # Convention: scan-project-{authoring_stem}.json in same directory as authoring file
                scan_project_path = self.authoring_path.parent / f"negaccel-scan-project-{self.authoring_path.stem}.json"
            
            self.scan_manager_window = ScanManagerWindow(
                parent_window=self,
                authoring_path=self.authoring_path,
                scan_project_path=scan_project_path if scan_project_path and scan_project_path.exists() else None,
            )
            if self.scan_manager_window.scan_project is None:
                self.scan_manager_window.seed_from_main_window()
            self.scan_manager_window.show()
        else:
            self.scan_manager_window.raise_()
            self.scan_manager_window.activateWindow()

    def _decorate_new_numeric_rows(self, form_layout: QFormLayout, existing_widget_keys: set[str]) -> None:
        path_by_widget_id: dict[int, str] = {}
        for path in sorted(set(self.widgets) - existing_widget_keys):
            widget = self.widgets.get(path)
            if isinstance(widget, (QSpinBox, QDoubleSpinBox)):
                path_by_widget_id[id(widget)] = path

        for row in range(form_layout.rowCount()):
            label_item = form_layout.itemAt(row, QFormLayout.ItemRole.LabelRole)
            field_item = form_layout.itemAt(row, QFormLayout.ItemRole.FieldRole)
            if label_item is None or field_item is None:
                continue

            field_widget = field_item.widget()
            label_widget = label_item.widget()
            if field_widget is None or not isinstance(field_widget, (QSpinBox, QDoubleSpinBox)):
                continue

            path = path_by_widget_id.get(id(field_widget))
            if not path:
                continue

            label_text = label_widget.text() if isinstance(label_widget, QLabel) else path
            self.parameter_widgets_by_path[path] = field_widget
            self.parameter_labels_by_path[path] = label_text
            wrapped = self._build_parameterized_field(path, label_text, field_widget)
            form_layout.setWidget(row, QFormLayout.ItemRole.FieldRole, wrapped)

    def _build_parameterized_field(self, path: str, label_text: str, field_widget: QWidget) -> QWidget:
        wrapper = QWidget()
        row = QHBoxLayout(wrapper)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        row.addWidget(field_widget, 1)

        checkbox = self._build_parameter_checkbox(path, label_text, field_widget)
        row.addWidget(checkbox, 0)
        return wrapper

    def _build_parameter_checkbox(self, path: str, label_text: str, widget: QWidget) -> QCheckBox:
        self.parameter_widgets_by_path[path] = widget
        self.parameter_labels_by_path[path] = label_text
        checkbox = self.parameter_checkboxes_by_path.get(path)
        if checkbox is None:
            checkbox = QCheckBox("P")
            checkbox.setMaximumWidth(28)
            checkbox.toggled.connect(
                lambda checked, p=path, l=label_text, w=widget: self._on_parameter_flag_toggled(p, l, w, checked)
            )
            self.parameter_checkboxes_by_path[path] = checkbox
        checkbox.setToolTip(f"Include '{label_text}' in the Scan Manager parameter list")
        checkbox.blockSignals(True)
        checkbox.setChecked(self.parameter_flag_registry.is_marked(path))
        checkbox.blockSignals(False)
        return checkbox

    def _on_parameter_flag_toggled(self, path: str, label_text: str, widget: QWidget, checked: bool) -> None:
        if checked:
            if isinstance(widget, QSpinBox):
                self.parameter_flag_registry.mark_parameter(
                    path,
                    label_text,
                    "integer",
                    min_val=float(widget.minimum()),
                    max_val=float(widget.maximum()),
                    step=float(widget.singleStep()),
                )
            elif isinstance(widget, QDoubleSpinBox):
                self.parameter_flag_registry.mark_parameter(
                    path,
                    label_text,
                    "number",
                    min_val=float(widget.minimum()),
                    max_val=float(widget.maximum()),
                    step=float(widget.singleStep()),
                )
        else:
            self.parameter_flag_registry.unmark_parameter(path)

    def build_scan_project_seed(self) -> ScanProjectData:
        project = ScanProjectData()
        case_tag = self.widgets["caseTag"].text().strip() or "negaccel_scan"
        project.scanProjectTag = f"{case_tag}_scan"
        project.authoringCasePath = str(self.authoring_path) if self.authoring_path is not None else None
        project.parameters = self.parameter_flag_registry.get_flagged_parameters()

        if project.parameters:
            values = [self._current_parameter_value(parameter.path) for parameter in project.parameters]
            project.cases = [
                ScanCaseDefinition(
                    caseIndex=0,
                    caseLabel="Current GUI values",
                    parameterValues=values,
                )
            ]
        return project

    def _current_parameter_value(self, path: str) -> Any:
        widget = self.parameter_widgets_by_path.get(path)
        if isinstance(widget, QSpinBox):
            return int(widget.value())
        if isinstance(widget, QDoubleSpinBox):
            return float(widget.value())
        return None

    def _double_spin(self, minimum: float, maximum: float, decimals: int, step: float) -> QDoubleSpinBox:
        spin = QDoubleSpinBox()
        spin.setRange(minimum, maximum)
        spin.setDecimals(decimals)
        spin.setSingleStep(step)
        return spin

    def _connect_preview_updates(self) -> None:
        for widget in self.widgets.values():
            if isinstance(widget, QLineEdit):
                widget.textChanged.connect(self.schedule_preview_refresh)
            elif isinstance(widget, QTextEdit):
                widget.textChanged.connect(self.schedule_preview_refresh)
            elif isinstance(widget, QComboBox):
                widget.currentTextChanged.connect(self.schedule_preview_refresh)
            elif isinstance(widget, (QSpinBox, QDoubleSpinBox)):
                widget.valueChanged.connect(self.schedule_preview_refresh)
            elif isinstance(widget, QCheckBox):
                widget.toggled.connect(self.schedule_preview_refresh)