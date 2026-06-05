"""Gas and surface-interaction form section for the NegAccel GUI."""

from __future__ import annotations

from ..common import QCheckBox, QComboBox, QFormLayout, QLabel, QLineEdit, QSpinBox, WorkflowError, nested_get, set_combo_value


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    window.widgets["gas.profileName"] = QLineEdit()
    window.widgets["gas.species"] = QLineEdit()
    window.widgets["gas.path"] = QLineEdit()
    window.widgets["gas.useFor"] = QLineEdit()
    window.widgets["gas.useFor"].setPlaceholderText("stripping, secondaryEmission, diagnostics")
    window.widgets["stripping.enabled"] = QCheckBox("Enable stripping")
    window.widgets["stripping.generateSecondaries"] = QCheckBox("Generate secondaries")
    window.widgets["stripping.minimumZ"] = window._double_spin(0.0, 10.0, 6, 0.001)
    window.widgets["surface.enabled"] = QCheckBox("Enable surface interactions")
    surface_model = QComboBox()
    surface_model.setEditable(True)
    surface_model.addItems(["EAMCC"])
    window.widgets["surface.model"] = surface_model
    window.widgets["surface.minimumImpactZ"] = window._double_spin(0.0, 10.0, 6, 0.001)
    window.widgets["surface.maximumSecondaryElectrons"] = QSpinBox()
    window.widgets["surface.maximumSecondaryElectrons"].setRange(0, 100000)
    window.widgets["surface.secondaryElectronEnergyEV"] = window._double_spin(0.0, 1.0e5, 3, 0.5)
    window.widgets["surface.debug"] = QCheckBox("Surface debug")
    layout.addRow(QLabel("Density profile"))
    layout.addRow("Profile name", window.widgets["gas.profileName"])
    layout.addRow("Gas species", window.widgets["gas.species"])
    layout.addRow("Profile path", window.widgets["gas.path"])
    layout.addRow("Use for", window.widgets["gas.useFor"])
    layout.addRow(window.widgets["stripping.enabled"])
    layout.addRow(window.widgets["stripping.generateSecondaries"])
    layout.addRow("Stripping minimum z [m]", window.widgets["stripping.minimumZ"])
    layout.addRow(QLabel("Surface collisions"))
    layout.addRow(window.widgets["surface.enabled"])
    layout.addRow("Model", window.widgets["surface.model"])
    layout.addRow("Minimum impact z [m]", window.widgets["surface.minimumImpactZ"])
    layout.addRow("Max secondary electrons", window.widgets["surface.maximumSecondaryElectrons"])
    layout.addRow("Secondary electron energy [eV]", window.widgets["surface.secondaryElectronEnergyEV"])
    layout.addRow(window.widgets["surface.debug"])
    return layout


def populate(window, spec: dict[str, object]) -> None:
    window.widgets["gas.profileName"].setText(
        str(nested_get(spec, "gasInteractions", "densityProfile", "name", default="background"))
    )
    window.widgets["gas.species"].setText(str(nested_get(spec, "gasInteractions", "densityProfile", "species", default="H2")))
    window.widgets["gas.path"].setText(
        str(nested_get(spec, "gasInteractions", "densityProfile", "path", default="densprofiles/MTF_dens.txt"))
    )
    window.widgets["gas.useFor"].setText(", ".join(list(nested_get(spec, "gasInteractions", "densityProfile", "useFor", default=[]))))
    window.widgets["stripping.enabled"].setChecked(bool(nested_get(spec, "gasInteractions", "stripping", "enabled", default=False)))
    window.widgets["stripping.generateSecondaries"].setChecked(
        bool(nested_get(spec, "gasInteractions", "stripping", "generateSecondaries", default=False))
    )
    window.widgets["stripping.minimumZ"].setValue(
        float(nested_get(spec, "gasInteractions", "stripping", "minimumZMeters", default=0.0))
    )
    window.widgets["surface.enabled"].setChecked(bool(nested_get(spec, "surfaceInteractions", "enabled", default=False)))
    set_combo_value(window.widgets["surface.model"], str(nested_get(spec, "surfaceInteractions", "model", default="EAMCC")))
    window.widgets["surface.minimumImpactZ"].setValue(
        float(nested_get(spec, "surfaceInteractions", "minimumImpactZMeters", default=0.0))
    )
    window.widgets["surface.maximumSecondaryElectrons"].setValue(
        int(nested_get(spec, "surfaceInteractions", "maximumSecondaryElectrons", default=0))
    )
    window.widgets["surface.secondaryElectronEnergyEV"].setValue(
        float(nested_get(spec, "surfaceInteractions", "secondaryElectronEnergyEV", default=0.0))
    )
    window.widgets["surface.debug"].setChecked(bool(nested_get(spec, "surfaceInteractions", "debug", default=False)))


def collect(window, spec: dict[str, object]) -> None:
    gas = spec.setdefault("gasInteractions", {})
    density_profile = gas.setdefault("densityProfile", {})
    density_profile["name"] = window.widgets["gas.profileName"].text().strip() or "background"
    density_profile["species"] = window.widgets["gas.species"].text().strip() or "H2"
    density_profile["path"] = window.widgets["gas.path"].text().strip()
    use_for = [token.strip() for token in window.widgets["gas.useFor"].text().split(",") if token.strip()]
    if not use_for:
        raise WorkflowError("gasInteractions.densityProfile.useFor must not be empty")
    density_profile["useFor"] = use_for
    stripping = gas.setdefault("stripping", {})
    stripping["enabled"] = window.widgets["stripping.enabled"].isChecked()
    stripping["generateSecondaries"] = window.widgets["stripping.generateSecondaries"].isChecked()
    stripping["minimumZMeters"] = float(window.widgets["stripping.minimumZ"].value())

    surface = spec.setdefault("surfaceInteractions", {})
    surface["enabled"] = window.widgets["surface.enabled"].isChecked()
    surface["model"] = window.widgets["surface.model"].currentText().strip() or "EAMCC"
    surface["minimumImpactZMeters"] = float(window.widgets["surface.minimumImpactZ"].value())
    surface["maximumSecondaryElectrons"] = int(window.widgets["surface.maximumSecondaryElectrons"].value())
    surface["secondaryElectronEnergyEV"] = float(window.widgets["surface.secondaryElectronEnergyEV"].value())
    surface["debug"] = window.widgets["surface.debug"].isChecked()
