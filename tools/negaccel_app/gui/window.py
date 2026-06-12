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
    QLineEdit,
    QMainWindow,
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
    create_scrollable_form,
    default_runtime_path,
)
from .execution import ExecutionMixin
from .results import ResultsMixin
from .sections import FORM_SECTIONS, geometry as geometry_section


class NegAccelMainWindow(AuthoringMixin, ExecutionMixin, ResultsMixin, QMainWindow):
    def __init__(
        self,
        authoring_path: Path | None = None,
        runtime_path: Path | None = None,
    ) -> None:
        super().__init__()
        self.setWindowTitle("NegAccel Studio")
        self.resize(1560, 980)

        self.authoring_path: Path | None = None
        self.loaded_authoring_template = self._load_template(EXAMPLE_AUTHORING_PATH)
        self.current_runtime_case: dict[str, Any] | None = None
        self.current_runtime_path: Path | None = None
        self.preview_updates_blocked = False
        self.runtime_path_auto = runtime_path is None

        self.preview_timer = QTimer(self)
        self.preview_timer.setSingleShot(True)
        self.preview_timer.timeout.connect(self.refresh_runtime_preview)

        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.append_process_output)
        self.process.started.connect(self.on_process_started)
        self.process.finished.connect(self.on_process_finished)

        self.widgets: dict[str, QWidget] = {}
        self.output_tabs: QTabWidget | None = None

        self._build_ui(runtime_path)
        self._apply_style()

        if authoring_path is not None:
            self.load_authoring_from_path(authoring_path)
        else:
            self.populate_from_authoring_spec(self.loaded_authoring_template, EXAMPLE_AUTHORING_PATH)

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
                padding: 8px 12px;
                border-top-left-radius: 6px;
                border-top-right-radius: 6px;
                margin-right: 2px;
            }
            QTabBar::tab:selected {
                background: #fffaf2;
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

        execution_box = QGroupBox("Execution")
        execution_form = QFormLayout(execution_box)

        runtime_row = QHBoxLayout()
        self.runtime_path_edit = QLineEdit(str(runtime_path or default_runtime_path("negaccel_case")))
        self.runtime_path_edit.textEdited.connect(self.on_runtime_path_edited)
        runtime_browse = QPushButton("Browse")
        runtime_browse.clicked.connect(self.choose_runtime_path)
        runtime_row.addWidget(self.runtime_path_edit, 1)
        runtime_row.addWidget(runtime_browse)
        runtime_wrapper = QWidget()
        runtime_wrapper.setLayout(runtime_row)
        execution_form.addRow("Runtime case JSON", runtime_wrapper)

        self.load_existing_checkbox = QCheckBox("Load existing .dat outputs instead of solving a fresh case")
        execution_form.addRow("Run mode", self.load_existing_checkbox)

        root_layout.addWidget(execution_box)

        button_row = QHBoxLayout()
        self.load_button = QPushButton("Load Authoring JSON")
        self.load_button.clicked.connect(self.choose_authoring_file)
        self.save_button = QPushButton("Save Authoring JSON")
        self.save_button.clicked.connect(self.save_authoring_file)
        self.materialize_button = QPushButton("Materialize Runtime Case")
        self.materialize_button.clicked.connect(self.materialize_runtime_case)
        self.run_button = QPushButton("Run Simulation")
        self.run_button.clicked.connect(self.run_simulation)
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop_simulation)
        self.stop_button.setEnabled(False)
        for button in [self.load_button, self.save_button, self.materialize_button, self.run_button, self.stop_button]:
            button_row.addWidget(button)
        button_row.addStretch(1)
        root_layout.addLayout(button_row)

        self.output_tabs = self._build_main_tabs()
        root_layout.addWidget(self.output_tabs, 1)

        self.setCentralWidget(central)
        self.statusBar().showMessage("Ready")

    def _build_main_tabs(self) -> QTabWidget:
        tabs = QTabWidget()
        for section in FORM_SECTIONS:
            if section.title == "Metadata":
                tabs.addTab(self._build_metadata_tab(section), section.title)
            elif section.title == "Geometry":
                tabs.addTab(geometry_section.build_workspace(self), section.title)
            else:
                tabs.addTab(create_scrollable_form(section.build_form(self)), section.title)
        tabs.addTab(self._build_run_log_tab(), "Run Log")
        tabs.addTab(self._build_visualization_tab(), "Visualization")
        self._connect_preview_updates()
        return tabs

    def _build_metadata_tab(self, section) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)

        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(create_scrollable_form(section.build_form(self)))
        splitter.addWidget(self._build_json_preview_panel())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([250, 640])

        layout.addWidget(splitter)
        return page

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