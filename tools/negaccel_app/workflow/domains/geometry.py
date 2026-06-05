"""Geometry builders for workflow materialization."""

from __future__ import annotations

import copy
from typing import Any

from ..common import WorkflowError, default_domain_for_template, require_number


def build_runtime_geometry(authoring_geometry: dict[str, Any], template: str, outputs: dict[str, Any]) -> dict[str, Any]:
    mesh_size = require_number(authoring_geometry, "meshSizeMeters", "geometry")

    domain = authoring_geometry.get("domain")
    if domain is None:
        runtime_domain = default_domain_for_template(template)
    else:
        if not isinstance(domain, dict):
            raise WorkflowError("geometry.domain must be an object when provided")
        runtime_domain = {
            "xSizeMeters": require_number(domain, "xSizeMeters", "geometry.domain"),
            "ySizeMeters": require_number(domain, "ySizeMeters", "geometry.domain"),
            "zSizeMeters": require_number(domain, "zSizeMeters", "geometry.domain"),
        }

    geometry_export_vtk = bool(
        authoring_geometry.get("exportGeometryVtk", outputs.get("vtk", {}).get("exportGeometry", True))
    )

    runtime_geometry = {
        "source": {
            "mode": "builtin-generator",
            "template": template,
        },
        "mesh": {
            "sizeMeters": mesh_size,
            "exportGeometryVtk": geometry_export_vtk,
        },
        "domain": runtime_domain,
    }

    if "gaps" in authoring_geometry:
        if not isinstance(authoring_geometry["gaps"], dict):
            raise WorkflowError("geometry.gaps must be an object when provided")
        runtime_geometry["gaps"] = copy.deepcopy(authoring_geometry["gaps"])

    return runtime_geometry
