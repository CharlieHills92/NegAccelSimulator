"""Particle-source form section for the NegAccel GUI."""

from __future__ import annotations

from ..common import QComboBox, QFormLayout, QLineEdit, nested_get, set_combo_value


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    window.widgets["source.species"] = QLineEdit()
    window.widgets["source.chargeState"] = window._double_spin(-5.0, 5.0, 3, 0.1)
    window.widgets["source.massU"] = window._double_spin(0.0, 100.0, 4, 0.1)
    window.widgets["source.currentDensityAm2"] = window._double_spin(0.0, 1.0e6, 3, 1.0)
    window.widgets["source.perpendicularTemperatureEV"] = window._double_spin(0.0, 1.0e4, 4, 0.1)
    window.widgets["source.parallelTemperatureEV"] = window._double_spin(0.0, 1.0e4, 4, 0.1)
    window.widgets["source.axialEnergyEV"] = window._double_spin(0.0, 1.0e6, 4, 0.1)
    window.widgets["source.plasmaPotentialVolts"] = window._double_spin(-1.0e6, 1.0e6, 4, 1.0)
    window.widgets["source.electronsModelWeight"] = window._double_spin(0.0, 1.0e6, 4, 0.1)
    source_model = QComboBox()
    source_model.setEditable(True)
    source_model.addItems(["uniform"])
    window.widgets["source.sourceModel"] = source_model
    layout.addRow("Species", window.widgets["source.species"])
    layout.addRow("Charge state", window.widgets["source.chargeState"])
    layout.addRow("Mass [u]", window.widgets["source.massU"])
    layout.addRow("Current density [A/m^2]", window.widgets["source.currentDensityAm2"])
    layout.addRow("Perpendicular temperature [eV]", window.widgets["source.perpendicularTemperatureEV"])
    layout.addRow("Parallel temperature [eV]", window.widgets["source.parallelTemperatureEV"])
    layout.addRow("Axial energy [eV]", window.widgets["source.axialEnergyEV"])
    layout.addRow("Plasma potential [V]", window.widgets["source.plasmaPotentialVolts"])
    layout.addRow("Electrons model weight", window.widgets["source.electronsModelWeight"])
    layout.addRow("Source model", window.widgets["source.sourceModel"])
    return layout


def populate(window, spec: dict[str, object]) -> None:
    window.widgets["source.species"].setText(str(nested_get(spec, "particleSource", "species", default="H-")))
    window.widgets["source.chargeState"].setValue(float(nested_get(spec, "particleSource", "chargeState", default=-1.0)))
    window.widgets["source.massU"].setValue(float(nested_get(spec, "particleSource", "massU", default=1.0)))
    window.widgets["source.currentDensityAm2"].setValue(
        float(nested_get(spec, "particleSource", "currentDensityAm2", default=270.0))
    )
    window.widgets["source.perpendicularTemperatureEV"].setValue(
        float(nested_get(spec, "particleSource", "perpendicularTemperatureEV", default=0.0))
    )
    window.widgets["source.parallelTemperatureEV"].setValue(
        float(nested_get(spec, "particleSource", "parallelTemperatureEV", default=0.0))
    )
    window.widgets["source.axialEnergyEV"].setValue(float(nested_get(spec, "particleSource", "axialEnergyEV", default=3.0)))
    window.widgets["source.plasmaPotentialVolts"].setValue(
        float(nested_get(spec, "particleSource", "plasmaPotentialVolts", default=0.0))
    )
    window.widgets["source.electronsModelWeight"].setValue(
        float(nested_get(spec, "particleSource", "electronsModelWeight", default=0.0))
    )
    set_combo_value(window.widgets["source.sourceModel"], str(nested_get(spec, "particleSource", "sourceModel", default="uniform")))


def collect(window, spec: dict[str, object]) -> None:
    source = spec.setdefault("particleSource", {})
    source["species"] = window.widgets["source.species"].text().strip() or "H-"
    source["chargeState"] = float(window.widgets["source.chargeState"].value())
    source["massU"] = float(window.widgets["source.massU"].value())
    source["currentDensityAm2"] = float(window.widgets["source.currentDensityAm2"].value())
    source["perpendicularTemperatureEV"] = float(window.widgets["source.perpendicularTemperatureEV"].value())
    source["parallelTemperatureEV"] = float(window.widgets["source.parallelTemperatureEV"].value())
    source["axialEnergyEV"] = float(window.widgets["source.axialEnergyEV"].value())
    source["plasmaPotentialVolts"] = float(window.widgets["source.plasmaPotentialVolts"].value())
    source["electronsModelWeight"] = float(window.widgets["source.electronsModelWeight"].value())
    source["sourceModel"] = window.widgets["source.sourceModel"].currentText().strip() or "uniform"
