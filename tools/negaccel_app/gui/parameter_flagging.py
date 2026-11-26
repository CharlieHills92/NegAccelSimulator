"""Parameter flagging infrastructure for GUI-driven scan projects."""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

from .common import WorkflowError, load_json, write_json


@dataclass
class ParameterDefinition:
    """Represents a single flagged parameter."""
    path: str  # JSON pointer or dot notation
    label: str
    parameterType: str  # "number", "integer", "boolean", "string"
    min: Optional[float] = None
    max: Optional[float] = None
    step: Optional[float] = None
    unit: Optional[str] = None
    
    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ScanCaseDefinition:
    """Represents one scan case (one row in the scan table)."""
    caseIndex: int
    parameterValues: list[Any]
    caseLabel: str = ""
    status: str = "pending"  # pending, running, completed, failed, skipped
    runtimeJsonPath: Optional[str] = None
    outputDirectoryPath: Optional[str] = None
    errorMessage: Optional[str] = None
    completedAt: Optional[str] = None
    
    def to_dict(self) -> dict[str, Any]:
        return {k: v for k, v in asdict(self).items() if v is not None or k == "parameterValues"}


@dataclass
class ExecutionHistoryEntry:
    """Historical record of a scan execution."""
    runTag: str
    startedAt: str
    completedAt: Optional[str] = None
    mode: str = "sequential"  # sequential, parallel
    workerCount: int = 1
    caseIndices: list[int] = field(default_factory=list)
    
    def to_dict(self) -> dict[str, Any]:
        return {k: v for k, v in asdict(self).items() if v is not None or k == "caseIndices"}


@dataclass
class ScanProjectData:
    """In-memory representation of a scan project file."""
    schemaVersion: str = "1.0.0"
    scanProjectTag: str = "untitled_scan"
    authoringCasePath: Optional[str] = None
    parameters: list[ParameterDefinition] = field(default_factory=list)
    cases: list[ScanCaseDefinition] = field(default_factory=list)
    executionHistory: list[ExecutionHistoryEntry] = field(default_factory=list)
    defaultXMetric: Optional[str] = None
    defaultYMetric: Optional[str] = None
    favoriteMetrics: list[str] = field(default_factory=list)
    
    @staticmethod
    def from_file(path: Path) -> ScanProjectData:
        """Load scan project from JSON file."""
        try:
            data = load_json(path)
        except Exception as e:
            raise WorkflowError(f"Failed to load scan project {path}: {e}")
        
        project = ScanProjectData()
        if "metadata" in data:
            meta = data["metadata"]
            project.scanProjectTag = meta.get("scanProjectTag", "untitled_scan")
            project.authoringCasePath = meta.get("authoringCasePath")
        
        if "parameters" in data:
            project.parameters = [
                ParameterDefinition(**p) for p in data["parameters"]
            ]
        
        if "scanTable" in data and "cases" in data["scanTable"]:
            project.cases = [
                ScanCaseDefinition(**c) for c in data["scanTable"]["cases"]
            ]
        
        if "scanTable" in data and "executionHistory" in data["scanTable"]:
            project.executionHistory = [
                ExecutionHistoryEntry(**eh) for eh in data["scanTable"]["executionHistory"]
            ]
        
        if "visualization" in data:
            viz = data["visualization"]
            project.defaultXMetric = viz.get("defaultXMetric")
            project.defaultYMetric = viz.get("defaultYMetric")
            project.favoriteMetrics = viz.get("favoriteMetrics", [])
        
        return project
    
    def to_dict(self) -> dict[str, Any]:
        """Convert to JSON-serializable dictionary."""
        return {
            "schemaVersion": self.schemaVersion,
            "metadata": {
                "scanProjectTag": self.scanProjectTag,
                "authoringCasePath": self.authoringCasePath,
                "lastModified": datetime.now().isoformat() + "Z",
            },
            "parameters": [p.to_dict() for p in self.parameters],
            "scanTable": {
                "cases": [c.to_dict() for c in self.cases],
                "executionHistory": [eh.to_dict() for eh in self.executionHistory],
            },
            "visualization": {
                "defaultXMetric": self.defaultXMetric,
                "defaultYMetric": self.defaultYMetric,
                "favoriteMetrics": self.favoriteMetrics,
            }
        }
    
    def save(self, path: Path) -> None:
        """Write scan project to JSON file."""
        data = self.to_dict()
        write_json(path, data)


class ParameterFlagRegistry:
    """In-memory registry of flagged parameters for the current authoring case."""
    
    def __init__(self):
        self._flags: dict[str, ParameterDefinition] = {}
    
    def mark_parameter(self, path: str, label: str, parameterType: str, 
                       min_val: Optional[float] = None, max_val: Optional[float] = None,
                       step: Optional[float] = None, unit: Optional[str] = None) -> None:
        """Mark a field as a scannable parameter."""
        self._flags[path] = ParameterDefinition(
            path=path,
            label=label,
            parameterType=parameterType,
            min=min_val,
            max=max_val,
            step=step,
            unit=unit,
        )
    
    def unmark_parameter(self, path: str) -> None:
        """Unmark a previously flagged parameter."""
        self._flags.pop(path, None)
    
    def is_marked(self, path: str) -> bool:
        """Check if a field is marked as a parameter."""
        return path in self._flags
    
    def get_flagged_parameters(self) -> list[ParameterDefinition]:
        """Get all flagged parameters in insertion order."""
        return list(self._flags.values())
    
    def clear_all(self) -> None:
        """Clear all flagged parameters."""
        self._flags.clear()
    
    def load_from_project(self, project: ScanProjectData) -> None:
        """Load parameter definitions from a scan project."""
        self._flags = {p.path: p for p in project.parameters}
    
    def to_project_parameters(self) -> list[dict[str, Any]]:
        """Export flagged parameters as project data."""
        return [p.to_dict() for p in self._flags.values()]
