"""Magnetic-field builders for workflow materialization."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..common import REPO_ROOT, WorkflowError, require_string


RESERVED_FIELD_NAMES = {"all", "none", "default", "sum", "total"}


def _validate_name(name: str, context: str) -> str:
    normalized = name.strip()
    if not normalized:
        raise WorkflowError(f"{context}.name must be a non-empty string")
    lowered = normalized.lower()
    if lowered in RESERVED_FIELD_NAMES:
        raise WorkflowError(f"{context}.name uses reserved value: {normalized}")
    for char in normalized:
        if not (char.isalnum() or char in {"_", "-"}):
            raise WorkflowError(f"{context}.name may contain only letters, digits, '_' or '-'")
    return normalized


def _serialize_repo_relative_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return resolved.as_posix()


def _resolve_field_file_path(file_path: str, directory: str, base_path: Path | None, context: str) -> str:
    raw_file_path = file_path.strip()
    if not raw_file_path:
        raise WorkflowError(f"{context}.filePath is required")

    candidate_path = Path(raw_file_path).expanduser()
    if candidate_path.is_absolute():
        resolved = candidate_path.resolve()
        if not resolved.exists():
            raise WorkflowError(f"{context}.filePath does not exist: {resolved}")
        return _serialize_repo_relative_path(resolved)

    directory_path = Path(directory.strip()) if directory.strip() else None
    if directory_path is not None:
        try:
            candidate_path.relative_to(directory_path)
            combined = candidate_path
        except ValueError:
            combined = directory_path / candidate_path
    else:
        combined = candidate_path

    search_roots: list[Path] = []
    if base_path is not None:
        search_roots.append(base_path.parent)
    search_roots.append(REPO_ROOT)

    for root in search_roots:
        resolved = (root / combined).resolve()
        if resolved.exists():
            return _serialize_repo_relative_path(resolved)

    raise WorkflowError(f"{context}.filePath does not exist: {combined.as_posix()}")


def _build_new_model(authoring: dict[str, Any], base_path: Path | None) -> dict[str, Any]:
    directory = str(authoring.get("directory", "")).strip()
    enabled = bool(authoring.get("enabled", False))
    authored_fields = authoring.get("fields")
    if authored_fields is None:
        authored_fields = []
    if not isinstance(authored_fields, list):
        raise WorkflowError("magneticField.fields must be an array")

    runtime_fields: list[dict[str, Any]] = []
    seen_names: set[str] = set()
    for index, field in enumerate(authored_fields):
        context = f"magneticField.fields[{index}]"
        if not isinstance(field, dict):
            raise WorkflowError(f"{context} must be an object")
        name = _validate_name(require_string(field, "name", context), context)
        lowered_name = name.lower()
        if lowered_name in seen_names:
            raise WorkflowError(f"Duplicate magnetic field name: {name}")
        seen_names.add(lowered_name)

        source_type = require_string(field, "sourceType", context)
        if source_type not in {"constant", "file"}:
            raise WorkflowError(f"{context}.sourceType must be one of constant or file")

        scale = float(field.get("scale", 1.0))
        if scale < 0.0:
            raise WorkflowError(f"{context}.scale must be >= 0")

        runtime_entry: dict[str, Any] = {
            "name": name,
            "sourceType": source_type,
            "scale": scale,
        }
        if source_type == "constant":
            constant = field.get("constantValue")
            if not isinstance(constant, dict):
                raise WorkflowError(f"{context}.constantValue must be an object")
            runtime_entry["constantValue"] = {
                "bx": float(constant.get("bx", 0.0)),
                "by": float(constant.get("by", 0.0)),
                "bz": float(constant.get("bz", 0.0)),
            }
        else:
            runtime_entry["filePath"] = _resolve_field_file_path(
                require_string(field, "filePath", context),
                directory,
                base_path,
                context,
            )
        runtime_fields.append(runtime_entry)

    runtime_enabled = enabled and bool(runtime_fields)
    runtime: dict[str, Any] = {
        "enabled": runtime_enabled,
        "fields": runtime_fields,
    }
    return runtime


def _build_legacy_model(authoring: dict[str, Any], base_path: Path | None) -> dict[str, Any]:
    enabled = bool(authoring.get("enabled", False))
    if not enabled:
        return {
            "enabled": False,
            "fields": [],
        }

    source_mode = require_string(authoring, "sourceMode", "magneticField")
    if source_mode not in {"directory", "file"}:
        raise WorkflowError("magneticField.sourceMode must be one of directory or file")

    scale = float(authoring.get("scale", 1.0))
    if scale < 0.0:
        raise WorkflowError("magneticField.scale must be >= 0")

    runtime: dict[str, Any] = {
        "enabled": True,
        "fields": [],
    }
    if source_mode == "directory":
        directory = require_string(authoring, "directory", "magneticField")
        runtime["fields"].append(
            {
                "name": "external",
                "sourceType": "file",
                "scale": scale,
                "filePath": _resolve_field_file_path("EXTfield.fld", directory, base_path, "magneticField"),
            }
        )
    else:
        runtime["fields"].append(
            {
                "name": "external",
                "sourceType": "file",
                "scale": scale,
                "filePath": _resolve_field_file_path(
                    require_string(authoring, "file", "magneticField"),
                    "",
                    base_path,
                    "magneticField",
                ),
            }
        )
    return runtime


def build_magnetic_field(authoring: dict[str, Any], base_path: Path | None = None) -> dict[str, Any]:
    if not authoring:
        return {
            "enabled": False,
            "fields": [],
        }

    if isinstance(authoring.get("fields"), list):
        return _build_new_model(authoring, base_path)
    return _build_legacy_model(authoring, base_path)
