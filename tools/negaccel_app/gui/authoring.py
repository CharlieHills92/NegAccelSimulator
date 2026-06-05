"""Authoring, preview, and materialization helpers for the NegAccel GUI."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

from .common import (
    QFileDialog,
    QMessageBox,
    REPO_ROOT,
    WorkflowError,
    authored_to_runtime_case,
    default_runtime_path,
    ensure_case_metadata,
    load_json,
    sanitize_case_tag,
    write_json,
)
from .sections import FORM_SECTIONS


class AuthoringMixin:
    def _load_template(self, path: Path) -> dict[str, Any]:
        payload = load_json(path)
        if not isinstance(payload, dict):
            raise WorkflowError(f"Authoring template must be a JSON object: {path}")
        return payload

    def choose_authoring_file(self) -> None:
        selected, _ = QFileDialog.getOpenFileName(
            self, "Open authoring JSON", str(REPO_ROOT), "JSON files (*.json)"
        )
        if selected:
            self.load_authoring_from_path(Path(selected))

    def load_authoring_from_path(self, path: Path) -> None:
        try:
            payload = self._load_template(path)
        except Exception as exc:  # pragma: no cover - UI error path
            QMessageBox.critical(self, "Load failed", str(exc))
            return
        self.populate_from_authoring_spec(payload, path)
        self.append_log(f"Loaded authoring specification: {path}")

    def populate_from_authoring_spec(self, spec: dict[str, Any], source_path: Path | None) -> None:
        self.preview_updates_blocked = True
        self.loaded_authoring_template = copy.deepcopy(spec)
        self.authoring_path = source_path
        self.authoring_path_label.setText(str(source_path) if source_path else "Unsaved authoring specification")

        for section in FORM_SECTIONS:
            section.populate(self, spec)

        if self.runtime_path_auto:
            self.runtime_path_edit.setText(str(default_runtime_path(self.widgets["caseTag"].text())))

        self.preview_updates_blocked = False
        self.refresh_runtime_preview()

    def build_authoring_spec(self) -> dict[str, Any]:
        spec = copy.deepcopy(self.loaded_authoring_template)
        for section in FORM_SECTIONS:
            section.collect(self, spec)
        return spec

    def refresh_runtime_preview(self) -> None:
        if self.preview_updates_blocked:
            return

        try:
            authoring_spec = self.build_authoring_spec()
            runtime_case, runtime_path = self.build_runtime_case_preview(authoring_spec)
        except Exception as exc:
            self.authoring_preview.setPlainText(f"Preview unavailable:\n{exc}")
            self.runtime_preview.setPlainText(f"Preview unavailable:\n{exc}")
            self.statusBar().showMessage("Runtime preview contains validation errors")
            return

        self.current_runtime_case = runtime_case
        self.current_runtime_path = runtime_path
        self.authoring_preview.setPlainText(json.dumps(authoring_spec, indent=2))
        self.runtime_preview.setPlainText(json.dumps(runtime_case, indent=2))
        self.statusBar().showMessage(f"Runtime preview ready: {runtime_path}")

    def build_runtime_case_preview(self, authoring_spec: dict[str, Any] | None = None) -> tuple[dict[str, Any], Path]:
        working_spec = authoring_spec or self.build_authoring_spec()
        runtime_case = authored_to_runtime_case(working_spec)
        runtime_path = Path(
            self.runtime_path_edit.text().strip() or default_runtime_path(working_spec["metadata"]["caseTag"])
        )
        if not runtime_path.is_absolute():
            runtime_path = (REPO_ROOT / runtime_path).resolve()
        ensure_case_metadata(runtime_case, runtime_case["metadata"]["caseTag"], runtime_path.parent)
        return runtime_case, runtime_path

    def choose_runtime_path(self) -> None:
        suggested = self.runtime_path_edit.text().strip() or str(default_runtime_path(self.widgets["caseTag"].text()))
        selected, _ = QFileDialog.getSaveFileName(
            self, "Select runtime case JSON", suggested, "JSON files (*.json)"
        )
        if selected:
            self.runtime_path_auto = False
            self.runtime_path_edit.setText(selected)
            self.refresh_runtime_preview()

    def save_authoring_file(self) -> None:
        default_target = self.authoring_path or default_runtime_path(self.widgets["caseTag"].text()).with_name(
            f"{sanitize_case_tag(self.widgets['caseTag'].text())}_authoring.json"
        )
        selected, _ = QFileDialog.getSaveFileName(
            self, "Save authoring JSON", str(default_target), "JSON files (*.json)"
        )
        if not selected:
            return

        try:
            payload = self.build_authoring_spec()
            write_json(Path(selected), payload)
        except Exception as exc:  # pragma: no cover - UI error path
            QMessageBox.critical(self, "Save failed", str(exc))
            return

        self.authoring_path = Path(selected)
        self.authoring_path_label.setText(selected)
        self.loaded_authoring_template = copy.deepcopy(payload)
        self.authoring_preview.setPlainText(json.dumps(payload, indent=2))
        self.append_log(f"Saved authoring specification: {selected}")
        self.statusBar().showMessage(f"Saved authoring JSON to {selected}")

    def schedule_preview_refresh(self) -> None:
        if self.preview_updates_blocked:
            return
        self.preview_timer.start(160)

    def on_runtime_path_edited(self, _text: str) -> None:
        self.runtime_path_auto = False

    def on_case_tag_edited(self, _text: str) -> None:
        if self.runtime_path_auto:
            self.runtime_path_edit.setText(str(default_runtime_path(self.widgets["caseTag"].text())))

    def _set_optional_text(self, document: dict[str, Any], key: str, value: str) -> None:
        stripped = value.strip()
        if stripped:
            document[key] = stripped
        elif key in document:
            document.pop(key)