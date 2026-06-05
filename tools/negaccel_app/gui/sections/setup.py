"""Run-setup form section for the NegAccel GUI."""

from __future__ import annotations

from ..common import QCheckBox, QComboBox, QFormLayout, QSpinBox, nested_get, set_combo_value


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    window.widgets["run.particleCount"] = QSpinBox()
    window.widgets["run.particleCount"].setRange(1, 100000000)
    window.widgets["run.iterations"] = QSpinBox()
    window.widgets["run.iterations"].setRange(1, 10000)
    solver_type = QComboBox()
    solver_type.setEditable(True)
    solver_type.addItems(["bicgstab"])
    window.widgets["run.solverType"] = solver_type
    shield_model = QComboBox()
    shield_model.setEditable(True)
    shield_model.addItems(["nsimp"])
    window.widgets["run.shieldModel"] = shield_model
    window.widgets["run.alphaCoeff"] = window._double_spin(0.0, 10.0, 4, 0.05)
    window.widgets["run.pgFilterScale"] = window._double_spin(0.0, 100.0, 4, 0.1)
    window.widgets["run.cesmadcmScale"] = window._double_spin(0.0, 100.0, 4, 0.1)
    window.widgets["run.currentDensityTolerance"] = window._double_spin(0.0, 1000.0, 4, 0.1)
    window.widgets["run.splitDomain"] = QCheckBox("Enable domain split")
    layout.addRow("Particle count", window.widgets["run.particleCount"])
    layout.addRow("Iterations", window.widgets["run.iterations"])
    layout.addRow("Solver type", window.widgets["run.solverType"])
    layout.addRow("Shield model", window.widgets["run.shieldModel"])
    layout.addRow("Alpha coefficient", window.widgets["run.alphaCoeff"])
    layout.addRow("PG filter scale", window.widgets["run.pgFilterScale"])
    layout.addRow("CESMADCm scale", window.widgets["run.cesmadcmScale"])
    layout.addRow("Current density tolerance", window.widgets["run.currentDensityTolerance"])
    layout.addRow(window.widgets["run.splitDomain"])
    return layout


def populate(window, spec: dict[str, object]) -> None:
    window.widgets["run.particleCount"].setValue(int(nested_get(spec, "run", "particleCount", default=50000)))
    window.widgets["run.iterations"].setValue(int(nested_get(spec, "run", "iterations", default=4)))
    set_combo_value(window.widgets["run.solverType"], str(nested_get(spec, "run", "solver", "type", default="bicgstab")))
    set_combo_value(window.widgets["run.shieldModel"], str(nested_get(spec, "run", "solver", "shieldModel", default="nsimp")))
    window.widgets["run.alphaCoeff"].setValue(float(nested_get(spec, "run", "spaceCharge", "alphaCoeff", default=0.5)))
    window.widgets["run.pgFilterScale"].setValue(float(nested_get(spec, "run", "spaceCharge", "pgFilterScale", default=0.0)))
    window.widgets["run.cesmadcmScale"].setValue(float(nested_get(spec, "run", "spaceCharge", "cesmadcmScale", default=0.0)))
    window.widgets["run.currentDensityTolerance"].setValue(
        float(nested_get(spec, "run", "convergence", "currentDensityTolerance", default=1.0))
    )
    window.widgets["run.splitDomain"].setChecked(
        bool(nested_get(spec, "run", "domainDecomposition", "splitDomain", default=False))
    )


def collect(window, spec: dict[str, object]) -> None:
    run = spec.setdefault("run", {})
    run["particleCount"] = int(window.widgets["run.particleCount"].value())
    run["iterations"] = int(window.widgets["run.iterations"].value())
    solver = run.setdefault("solver", {})
    solver["type"] = window.widgets["run.solverType"].currentText().strip() or "bicgstab"
    solver["shieldModel"] = window.widgets["run.shieldModel"].currentText().strip() or "nsimp"
    space_charge = run.setdefault("spaceCharge", {})
    space_charge["alphaCoeff"] = float(window.widgets["run.alphaCoeff"].value())
    space_charge["pgFilterScale"] = float(window.widgets["run.pgFilterScale"].value())
    space_charge["cesmadcmScale"] = float(window.widgets["run.cesmadcmScale"].value())
    convergence = run.setdefault("convergence", {})
    convergence["currentDensityTolerance"] = float(window.widgets["run.currentDensityTolerance"].value())
    domain_decomposition = run.setdefault("domainDecomposition", {})
    domain_decomposition["splitDomain"] = window.widgets["run.splitDomain"].isChecked()
