"""Surface-interaction form section for the NegAccel GUI."""

from __future__ import annotations

from ..common import QCheckBox, QComboBox, QFormLayout, QSpinBox, nested_get, set_combo_value


def build_form(window) -> QFormLayout:
    layout = QFormLayout()

    window.widgets["surfaceInteractions.enabled"] = QCheckBox("Enable surface interactions")

    model = QComboBox()
    model.setEditable(True)
    model.addItems(["EAMCC"])
    window.widgets["surfaceInteractions.model"] = model

    minimum_impact_z = window._double_spin(0.0, 10.0, 6, 0.001)
    window.widgets["surfaceInteractions.minimumImpactZMeters"] = minimum_impact_z

    max_secondary_electrons = QSpinBox()
    max_secondary_electrons.setRange(0, 100000)
    window.widgets["surfaceInteractions.maximumSecondaryElectrons"] = max_secondary_electrons

    secondary_energy = window._double_spin(0.0, 1.0e5, 3, 0.5)
    window.widgets["surfaceInteractions.secondaryElectronEnergyEV"] = secondary_energy

    window.widgets["surfaceInteractions.debug"] = QCheckBox("Surface debug")

    layout.addRow(window.widgets["surfaceInteractions.enabled"])
    layout.addRow("Model", window.widgets["surfaceInteractions.model"])
    layout.addRow("Minimum impact z [m]", window.widgets["surfaceInteractions.minimumImpactZMeters"])
    layout.addRow("Max secondary electrons", window.widgets["surfaceInteractions.maximumSecondaryElectrons"])
    layout.addRow("Secondary electron energy [eV]", window.widgets["surfaceInteractions.secondaryElectronEnergyEV"])
    layout.addRow(window.widgets["surfaceInteractions.debug"])

    window.widgets["surfaceInteractions.enabled"].toggled.connect(window.schedule_preview_refresh)
    window.widgets["surfaceInteractions.model"].currentTextChanged.connect(window.schedule_preview_refresh)
    window.widgets["surfaceInteractions.minimumImpactZMeters"].valueChanged.connect(window.schedule_preview_refresh)
    window.widgets["surfaceInteractions.maximumSecondaryElectrons"].valueChanged.connect(window.schedule_preview_refresh)
    window.widgets["surfaceInteractions.secondaryElectronEnergyEV"].valueChanged.connect(window.schedule_preview_refresh)
    window.widgets["surfaceInteractions.debug"].toggled.connect(window.schedule_preview_refresh)
    return layout


def populate(window, spec: dict[str, object]) -> None:
    surface = nested_get(spec, "surfaceInteractions", default={})
    if not isinstance(surface, dict):
        surface = {}

    window.widgets["surfaceInteractions.enabled"].setChecked(bool(surface.get("enabled", False)))
    set_combo_value(window.widgets["surfaceInteractions.model"], str(surface.get("model") or "EAMCC"))
    window.widgets["surfaceInteractions.minimumImpactZMeters"].setValue(
        float(surface.get("minimumImpactZMeters") or 0.0)
    )
    window.widgets["surfaceInteractions.maximumSecondaryElectrons"].setValue(
        int(surface.get("maximumSecondaryElectrons") or 0)
    )
    window.widgets["surfaceInteractions.secondaryElectronEnergyEV"].setValue(
        float(surface.get("secondaryElectronEnergyEV") or 0.0)
    )
    window.widgets["surfaceInteractions.debug"].setChecked(bool(surface.get("debug", False)))


def collect(window, spec: dict[str, object]) -> None:
    spec["surfaceInteractions"] = {
        "enabled": window.widgets["surfaceInteractions.enabled"].isChecked(),
        "model": window.widgets["surfaceInteractions.model"].currentText().strip() or "EAMCC",
        "minimumImpactZMeters": float(window.widgets["surfaceInteractions.minimumImpactZMeters"].value()),
        "maximumSecondaryElectrons": int(window.widgets["surfaceInteractions.maximumSecondaryElectrons"].value()),
        "secondaryElectronEnergyEV": float(window.widgets["surfaceInteractions.secondaryElectronEnergyEV"].value()),
        "debug": window.widgets["surfaceInteractions.debug"].isChecked(),
    }