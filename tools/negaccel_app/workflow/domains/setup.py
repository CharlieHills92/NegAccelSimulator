"""Run-setup builders for workflow materialization."""

from __future__ import annotations

from typing import Any

from ..common import WorkflowError, require_number


def build_simulation(run: dict[str, Any]) -> dict[str, Any]:
    domain_decomposition = run.get("domainDecomposition", {})
    if not isinstance(domain_decomposition, dict):
        raise WorkflowError("run.domainDecomposition must be an object when provided")

    simulation = {
        "particleCount": int(require_number(run, "particleCount", "run")),
        "iterations": int(require_number(run, "iterations", "run")),
        "domainDecomposition": {
            "splitDomain": bool(domain_decomposition.get("splitDomain", False)),
        },
        "solver": {
            "type": run.get("solver", {}).get("type", "bicgstab"),
            "shieldModel": run.get("solver", {}).get("shieldModel", "nsimp"),
        },
        "spaceCharge": {
            "alphaCoeff": float(run.get("spaceCharge", {}).get("alphaCoeff", 0.5)),
            "pgFilterScale": float(run.get("spaceCharge", {}).get("pgFilterScale", 0.0)),
            "cesmadcmScale": float(run.get("spaceCharge", {}).get("cesmadcmScale", 0.0)),
        },
        "convergence": {
            "currentDensityTolerance": float(run.get("convergence", {}).get("currentDensityTolerance", 1.0)),
        },
    }
    if "domainIndex" in domain_decomposition:
        simulation["domainDecomposition"]["domainIndex"] = int(domain_decomposition["domainIndex"])
    return simulation
