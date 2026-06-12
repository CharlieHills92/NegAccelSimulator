"""Magnetic-field builders for workflow materialization."""

from __future__ import annotations

from typing import Any

from ..common import WorkflowError, require_string


def build_magnetic_field(authoring: dict[str, Any]) -> dict[str, Any]:
    if not authoring:
        return {
            "enabled": False,
            "sourceMode": "none",
            "scale": 1.0,
        }

    enabled = bool(authoring.get("enabled", False))
    if not enabled:
        return {
            "enabled": False,
            "sourceMode": "none",
            "scale": float(authoring.get("scale", 1.0)),
        }

    source_mode = require_string(authoring, "sourceMode", "magneticField")
    if source_mode not in {"directory", "file"}:
        raise WorkflowError(
            "magneticField.sourceMode must be one of directory or file"
        )

    runtime_field: dict[str, Any] = {
        "enabled": True,
        "sourceMode": source_mode,
        "scale": float(authoring.get("scale", 1.0)),
    }
    if source_mode == "directory" and "case" in authoring:
        runtime_field["case"] = int(authoring["case"])

    if source_mode == "directory":
        runtime_field["directory"] = require_string(authoring, "directory", "magneticField")
    if source_mode == "file":
        runtime_field["file"] = require_string(authoring, "file", "magneticField")

    return runtime_field
