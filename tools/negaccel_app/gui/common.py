"""Shared GUI helpers and dependency imports for NegAccel tooling."""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Any

try:
    from PySide6.QtCore import QProcess, Qt, QTimer, qInstallMessageHandler
    from PySide6.QtGui import QFontDatabase, QTextCursor
    from PySide6.QtWidgets import (
        QApplication,
        QCheckBox,
        QComboBox,
        QDoubleSpinBox,
        QFileDialog,
        QFormLayout,
        QGroupBox,
        QHBoxLayout,
        QLabel,
        QLineEdit,
        QListWidget,
        QListWidgetItem,
        QMainWindow,
        QMessageBox,
        QPlainTextEdit,
        QProgressBar,
        QPushButton,
        QSizePolicy,
        QScrollArea,
        QSpinBox,
        QStackedWidget,
        QSplitter,
        QTabWidget,
        QTextEdit,
        QVBoxLayout,
        QWidget,
    )
except ImportError as exc:  # pragma: no cover - depends on the local environment
    raise SystemExit(
        "PySide6 is required for the GUI. Install the packages listed in tools/gui_requirements.txt."
    ) from exc

from negaccel_app.workflow import (
    REPO_ROOT,
    WorkflowError,
    authored_to_runtime_case,
    ensure_case_metadata,
    load_json,
    serialize_geometry_file_path,
    sanitize_case_tag,
    write_json,
)


EXAMPLE_AUTHORING_PATH = REPO_ROOT / "negaccel-authoring.example.json"
DEFAULT_SIMULATOR_PATH = REPO_ROOT / "NegAccelExec"
DEFAULT_SETUP_SCRIPT = REPO_ROOT / "setup_environment.sh"
LOG_LEVELS = ["critical", "error", "warning", "info", "debug", "trace"]
MAGNETIC_SOURCE_MODES = ["directory", "file", "none"]
_qt_message_handler_installed = False
_previous_qt_message_handler = None


def _filtered_qt_message_handler(message_type, context, message) -> None:
    if message.startswith("QPainter::end: Painter ended with ") and message.endswith(" saved states"):
        return

    if _previous_qt_message_handler is not None:
        _previous_qt_message_handler(message_type, context, message)
        return

    sys.stderr.write(f"{message}\n")


def install_qt_message_filter() -> None:
    global _qt_message_handler_installed, _previous_qt_message_handler
    if _qt_message_handler_installed:
        return

    _previous_qt_message_handler = qInstallMessageHandler(_filtered_qt_message_handler)
    _qt_message_handler_installed = True


def nested_get(document: dict[str, Any], *path: str, default: Any = None) -> Any:
    current: Any = document
    for token in path:
        if not isinstance(current, dict) or token not in current:
            return default
        current = current[token]
    return current


def parse_number_list(raw_text: str, field_name: str) -> list[float]:
    values: list[float] = []
    for token in raw_text.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        try:
            values.append(float(stripped))
        except ValueError as exc:
            raise WorkflowError(f"{field_name} must contain only comma-separated numeric values") from exc
    return values


def format_number_list(values: list[Any]) -> str:
    return ", ".join(f"{float(value):g}" for value in values)


def resolve_runtime_parent_directory(parent_directory: str | Path | None = None) -> Path:
    if isinstance(parent_directory, Path):
        candidate = parent_directory
    elif isinstance(parent_directory, str) and parent_directory.strip():
        candidate = Path(parent_directory.strip()).expanduser()
    else:
        candidate = REPO_ROOT

    if not candidate.is_absolute():
        candidate = (REPO_ROOT / candidate).resolve()
    else:
        candidate = candidate.resolve()
    return candidate


def build_runtime_path(case_tag: str, parent_directory: str | Path | None = None) -> Path:
    sanitized = sanitize_case_tag(case_tag or "negaccel_case")
    root_directory = resolve_runtime_parent_directory(parent_directory)
    return root_directory / sanitized / f"{sanitized}.json"


def default_runtime_path(case_tag: str) -> Path:
    return build_runtime_path(case_tag)


def runtime_path_parent_seed(runtime_path: Path | None) -> Path:
    if runtime_path is None:
        return REPO_ROOT

    resolved_path = runtime_path.resolve()
    if resolved_path.parent.name == resolved_path.stem:
        return resolved_path.parent.parent
    return resolved_path.parent


def create_scrollable_form(form_layout: QFormLayout) -> QScrollArea:
    container = QWidget()
    container.setLayout(form_layout)

    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setWidget(container)
    return scroll


def set_combo_value(combo: QComboBox, value: str) -> None:
    index = combo.findText(value)
    if index >= 0:
        combo.setCurrentIndex(index)
    elif combo.isEditable():
        combo.setEditText(value)
