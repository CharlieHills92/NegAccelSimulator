"""Authoring, preview, and materialization helpers for the NegAccel GUI."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

from .common import (
    build_runtime_path,
    QFileDialog,
    QMessageBox,
    nested_get,
    REPO_ROOT,
    resolve_runtime_parent_directory,
    WorkflowError,
    authored_to_runtime_case,
    ensure_case_metadata,
    load_json,
    sanitize_case_tag,
    write_json,
)
from .sections import FORM_SECTIONS


class AuthoringMixin:
    def _resolved_output_parent_directory(self) -> Path:
        widget = self.widgets.get("execution.outputParentDirectory") if isinstance(self.widgets, dict) else None
        raw_value = widget.text().strip() if widget is not None and hasattr(widget, "text") else ""
        if raw_value:
            return resolve_runtime_parent_directory(raw_value)
        return resolve_runtime_parent_directory(getattr(self, "default_output_parent_directory", REPO_ROOT))

    def _derived_runtime_path(self, case_tag: str | None = None) -> Path:
        resolved_case_tag = case_tag if case_tag is not None else self.widgets["caseTag"].text()
        return build_runtime_path(resolved_case_tag, self._resolved_output_parent_directory())

    def _sync_runtime_path_preview(self, case_tag: str | None = None) -> Path:
        runtime_path = self._derived_runtime_path(case_tag)
        self.runtime_path_edit.setText(str(runtime_path))
        return runtime_path

    def choose_output_parent_directory(self) -> None:
        start_dir = str(self._resolved_output_parent_directory())
        selected = QFileDialog.getExistingDirectory(self, "Select simulation parent folder", start_dir)
        if selected:
            self.widgets["execution.outputParentDirectory"].setText(str(Path(selected).resolve()))

    def on_output_parent_directory_changed(self, _text: str) -> None:
        self._sync_runtime_path_preview()

    def _authoring_filename_case_tag(self, path: Path | None = None) -> str | None:
        target_path = self.authoring_path if path is None else path
        if target_path is None:
            return None
        return target_path.stem

    def _refresh_case_tag_from_json_state(self) -> None:
        checkbox = getattr(self, "case_tag_from_json_checkbox", None)
        case_tag_widget = self.widgets.get("caseTag") if isinstance(self.widgets, dict) else None
        if checkbox is None or case_tag_widget is None:
            return

        has_authoring_path = self.authoring_path is not None
        checkbox.blockSignals(True)
        checkbox.setEnabled(has_authoring_path)
        if not has_authoring_path:
            checkbox.setChecked(False)
        checkbox.blockSignals(False)
        case_tag_widget.setReadOnly(bool(checkbox.isChecked() and has_authoring_path))

    def _apply_case_tag_from_authoring_filename(self) -> None:
        case_tag = self._authoring_filename_case_tag()
        case_tag_widget = self.widgets.get("caseTag") if isinstance(self.widgets, dict) else None
        if case_tag is None or case_tag_widget is None:
            return
        case_tag_widget.setText(case_tag)
        self.on_case_tag_edited(case_tag)

    def on_case_tag_from_json_toggled(self, enabled: bool) -> None:
        case_tag_widget = self.widgets.get("caseTag") if isinstance(self.widgets, dict) else None
        if case_tag_widget is None:
            return

        if enabled:
            self._apply_case_tag_from_authoring_filename()
            case_tag_widget.setReadOnly(True)
        else:
            case_tag_widget.setReadOnly(False)

        self.schedule_preview_refresh()

    def _update_authoring_path_display(self) -> None:
        display = getattr(self, "authoring_path_display", None)
        if display is None:
            return
        if self.authoring_path is not None:
            display.setText(str(self.authoring_path))
        else:
            display.setText("Unsaved in-memory authoring spec")

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
        self._update_authoring_path_display()

        for section in FORM_SECTIONS:
            section.populate(self, spec)

        self._refresh_case_tag_from_json_state()
        if getattr(self, "case_tag_from_json_checkbox", None) is not None and self.case_tag_from_json_checkbox.isChecked():
            self._apply_case_tag_from_authoring_filename()

        configured_parent_directory = nested_get(spec, "execution", "outputParentDirectory")
        if isinstance(configured_parent_directory, str) and configured_parent_directory.strip():
            self.widgets["execution.outputParentDirectory"].setText(
                str(resolve_runtime_parent_directory(configured_parent_directory))
            )
        else:
            self.widgets["execution.outputParentDirectory"].setText(str(self.default_output_parent_directory))

        if self.runtime_path_auto:
            self._sync_runtime_path_preview()

        self.preview_updates_blocked = False
        self.refresh_runtime_preview()

    def build_authoring_spec(self) -> dict[str, Any]:
        spec = copy.deepcopy(self.loaded_authoring_template)
        for section in FORM_SECTIONS:
            section.collect(self, spec)
        execution = spec.setdefault("execution", {})
        execution["outputParentDirectory"] = str(self._resolved_output_parent_directory())
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
        working_spec = copy.deepcopy(authoring_spec) if authoring_spec is not None else self.build_authoring_spec()
        geometry_workspace = self.widgets.get("geometry.solidsEditor")
        current_geometry_document = getattr(geometry_workspace, "current_geometry_document", None)
        if callable(current_geometry_document):
            working_spec["geometry"] = current_geometry_document()
        if self.authoring_path is not None:
            working_spec["__source_path"] = self.authoring_path
        runtime_case = authored_to_runtime_case(working_spec)
        runtime_path = self._derived_runtime_path(working_spec["metadata"]["caseTag"])
        ensure_case_metadata(runtime_case, runtime_case["metadata"]["caseTag"], runtime_path.parent)
        return runtime_case, runtime_path

    def _default_authoring_save_target(self) -> Path:
        case_tag = self.widgets["caseTag"].text()
        return self.authoring_path or self._derived_runtime_path(case_tag).with_name(
            f"{sanitize_case_tag(case_tag)}_authoring.json"
        )

    def _build_authoring_payload_for_path(self, target_path: Path) -> dict[str, Any]:
        payload = self.build_authoring_spec()
        if getattr(self, "case_tag_from_json_checkbox", None) is not None and self.case_tag_from_json_checkbox.isChecked():
            metadata = payload.setdefault("metadata", {})
            metadata["caseTag"] = self._authoring_filename_case_tag(target_path) or metadata.get("caseTag", "")
        return payload

    def _write_authoring_file(self, target_path: Path) -> dict[str, Any]:
        payload = self._build_authoring_payload_for_path(target_path)
        write_json(target_path, payload)

        self.authoring_path = target_path
        self.loaded_authoring_template = copy.deepcopy(payload)
        self._update_authoring_path_display()
        self._refresh_case_tag_from_json_state()
        if getattr(self, "case_tag_from_json_checkbox", None) is not None and self.case_tag_from_json_checkbox.isChecked():
            self._apply_case_tag_from_authoring_filename()
        self.authoring_preview.setPlainText(json.dumps(payload, indent=2))
        return payload

    def auto_save_authoring_file(self) -> Path:
        target_path = self._default_authoring_save_target()
        self._write_authoring_file(target_path)
        self.append_log(f"Saved authoring specification: {target_path}")
        self.statusBar().showMessage(f"Saved authoring JSON to {target_path}")
        return target_path

    def save_authoring_file(self) -> None:
        default_target = self._default_authoring_save_target()
        selected, _ = QFileDialog.getSaveFileName(
            self, "Save authoring JSON", str(default_target), "JSON files (*.json)"
        )
        if not selected:
            return

        try:
            selected_path = Path(selected)
            self._write_authoring_file(selected_path)
        except Exception as exc:  # pragma: no cover - UI error path
            QMessageBox.critical(self, "Save failed", str(exc))
            return

        self.append_log(f"Saved authoring specification: {selected}")
        self.statusBar().showMessage(f"Saved authoring JSON to {selected}")

    def schedule_preview_refresh(self) -> None:
        if self.preview_updates_blocked:
            return
        self.preview_timer.start(160)

    def on_case_tag_edited(self, _text: str) -> None:
        if self.runtime_path_auto:
            self._sync_runtime_path_preview()

    def _set_optional_text(self, document: dict[str, Any], key: str, value: str) -> None:
        stripped = value.strip()
        if stripped:
            document[key] = stripped
        elif key in document:
            document.pop(key)