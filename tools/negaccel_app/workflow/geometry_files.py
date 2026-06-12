"""Geometry file helpers for NegAccel authoring."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .common import REPO_ROOT, WorkflowError, load_json, write_json
from .domains.geometry import _normalize_generated_solids


TEMPLATES_DIR = REPO_ROOT / "tools" / "Templates"


def resolve_geometry_file_path(reference: str, base_path: Path | None = None) -> Path:
    raw_reference = str(reference or "").strip()
    if not raw_reference:
        raise WorkflowError("geometry.path must not be empty")

    path = Path(raw_reference)
    if path.is_absolute():
        return path.resolve()
    if base_path is not None:
        return (base_path.parent / path).resolve()
    return (REPO_ROOT / path).resolve()


def serialize_geometry_file_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return resolved.as_posix()


def _positive_number(payload: dict[str, Any], key: str, context: str) -> float:
    value = payload.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise WorkflowError(f"{context}.{key} must be numeric")
    number = float(value)
    if number <= 0.0:
        raise WorkflowError(f"{context}.{key} must be > 0")
    return number


def _number(payload: dict[str, Any], key: str, context: str) -> float:
    value = payload.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise WorkflowError(f"{context}.{key} must be numeric")
    return float(value)


def normalize_geometry_file_document(payload: Any, context: str) -> dict[str, Any]:
    if not isinstance(payload, dict):
        raise WorkflowError(f"{context} must be a JSON object")

    normalized = {
        "name": str(payload.get("name") or "Geometry"),
        "meshSizeMeters": _positive_number(payload, "meshSizeMeters", context),
        "exportGeometryVtk": bool(payload.get("exportGeometryVtk", True)),
    }

    domain = payload.get("domain")
    if not isinstance(domain, dict):
        raise WorkflowError(f"{context}.domain must be an object")
    normalized["domain"] = {
        "xSizeMeters": _positive_number(domain, "xSizeMeters", f"{context}.domain"),
        "ySizeMeters": _positive_number(domain, "ySizeMeters", f"{context}.domain"),
        "zSizeMeters": _positive_number(domain, "zSizeMeters", f"{context}.domain"),
    }
    if "zStartMeters" in domain:
        normalized["domain"]["zStartMeters"] = _number(domain, "zStartMeters", f"{context}.domain")

    gaps = payload.get("gaps")
    if gaps is not None:
        if not isinstance(gaps, dict):
            raise WorkflowError(f"{context}.gaps must be an object when provided")
        if "accelerationGapMeters" in gaps:
            raise WorkflowError(
                f"{context}.gaps.accelerationGapMeters is no longer supported; remove it from the geometry file"
            )
        normalized_gaps: dict[str, float] = {}
        for key in ("extractionGapMeters",):
            if key not in gaps:
                continue
            normalized_gaps[key] = _positive_number(gaps, key, f"{context}.gaps")
        if normalized_gaps:
            normalized["gaps"] = normalized_gaps

    solids = payload.get("solids", [])
    normalized["solids"] = _normalize_generated_solids({"solids": solids}) if solids else []
    return normalized


def load_geometry_file(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    return normalize_geometry_file_document(load_json(resolved), resolved.as_posix())


def load_geometry_reference(reference: str, base_path: Path | None = None) -> tuple[Path, dict[str, Any]]:
    path = resolve_geometry_file_path(reference, base_path)
    return path, load_geometry_file(path)


def write_geometry_file(path: Path, payload: dict[str, Any]) -> Path:
    normalized = normalize_geometry_file_document(payload, path.resolve().as_posix())
    write_json(path.resolve(), normalized)
    return path.resolve()