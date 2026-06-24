"""Run-setup builders for workflow materialization."""

from __future__ import annotations

from math import floor
from typing import Any

from ..common import WorkflowError, require_number


def _copy_solver_block(solver: dict[str, Any], block_name: str) -> dict[str, Any] | None:
    block = solver.get(block_name)
    if block is None:
        return None
    if not isinstance(block, dict):
        raise WorkflowError(f"run.solver.{block_name} must be an object when provided")
    return dict(block)


def _mesh_node_count(size_meters: float, mesh_size_meters: float) -> int:
    return int(floor(size_meters / mesh_size_meters)) + 1


def multigrid_level_support(
    mesh_size_meters: float,
    x_size_meters: float,
    y_size_meters: float,
    z_size_meters: float,
) -> dict[str, Any]:
    node_counts = (
        _mesh_node_count(x_size_meters, mesh_size_meters),
        _mesh_node_count(y_size_meters, mesh_size_meters),
        _mesh_node_count(z_size_meters, mesh_size_meters),
    )

    max_levels = 1
    sizes = node_counts
    while True:
        if any(size % 2 == 0 for size in sizes):
            break
        next_sizes = tuple((size + 1) // 2 for size in sizes)
        if next_sizes == sizes:
            break
        max_levels += 1
        sizes = next_sizes

    return {
        "nodeCounts": node_counts,
        "maxLevels": max_levels,
    }


def multigrid_level_support_for_runtime_geometry(runtime_geometry: dict[str, Any]) -> dict[str, Any]:
    mesh = runtime_geometry.get("mesh")
    if not isinstance(mesh, dict):
        raise WorkflowError("geometry.mesh must be an object")
    domain = runtime_geometry.get("domain")
    if not isinstance(domain, dict):
        raise WorkflowError("geometry.domain must be an object")

    mesh_size_meters = require_number(mesh, "sizeMeters", "geometry.mesh")
    x_size_meters = require_number(domain, "xSizeMeters", "geometry.domain")
    y_size_meters = require_number(domain, "ySizeMeters", "geometry.domain")
    z_size_meters = require_number(domain, "zSizeMeters", "geometry.domain")
    return multigrid_level_support(mesh_size_meters, x_size_meters, y_size_meters, z_size_meters)


def build_simulation(
    run: dict[str, Any],
    particle_count: int | None = None,
    plasma: dict[str, Any] | None = None,
    runtime_geometry: dict[str, Any] | None = None,
) -> dict[str, Any]:
    domain_decomposition = run.get("domainDecomposition", {})
    if not isinstance(domain_decomposition, dict):
        raise WorkflowError("run.domainDecomposition must be an object when provided")

    solver = run.get("solver", {})
    if not isinstance(solver, dict):
        raise WorkflowError("run.solver must be an object when provided")

    simulation = {
        "particleCount": int(particle_count) if particle_count is not None else int(require_number(run, "particleCount", "run")),
        "iterations": int(require_number(run, "iterations", "run")),
        "domainDecomposition": {
            "splitDomain": bool(domain_decomposition.get("splitDomain", False)),
        },
        "solver": {
            "type": solver.get("type", "bicgstab"),
        },
        "spaceCharge": {
            "alphaCoeff": float(run.get("spaceCharge", {}).get("alphaCoeff", 0.5)),
        },
        "convergence": {
            "currentDensityTolerance": float(run.get("convergence", {}).get("currentDensityTolerance", 1.0)),
        },
    }
    if plasma is None:
        simulation["solver"]["shieldModel"] = solver.get("shieldModel", "nsimp")
    else:
        simulation["solver"].update(plasma)

    selected_solver_type = str(simulation["solver"]["type"])
    if selected_solver_type == "bicgstab":
        bicgstab = _copy_solver_block(solver, "bicgstab")
        if bicgstab:
            simulation["solver"]["bicgstab"] = bicgstab
    elif selected_solver_type == "multigrid":
        multigrid = _copy_solver_block(solver, "multigrid")
        if multigrid:
            if runtime_geometry is not None and "levels" in multigrid:
                support = multigrid_level_support_for_runtime_geometry(runtime_geometry)
                requested_levels = int(multigrid["levels"])
                max_levels = int(support["maxLevels"])
                if requested_levels > max_levels:
                    node_counts = support["nodeCounts"]
                    raise WorkflowError(
                        "run.solver.multigrid.levels="
                        f"{requested_levels} exceeds the maximum allowed {max_levels} "
                        f"for mesh {node_counts[0]} x {node_counts[1]} x {node_counts[2]}"
                    )
            simulation["solver"]["multigrid"] = multigrid

    if "domainIndex" in domain_decomposition:
        simulation["domainDecomposition"]["domainIndex"] = int(domain_decomposition["domainIndex"])
    return simulation
