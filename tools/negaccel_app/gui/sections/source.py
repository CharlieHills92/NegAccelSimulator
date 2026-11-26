"""Particles form section for the NegAccel GUI."""

from __future__ import annotations

from copy import deepcopy

from negaccel_app.particles import (
    PARTICLE_FAMILIES,
    build_family_particle_types,
    detect_particle_family_from_kinds,
    get_default_source_particle_type_id,
    get_particle_kind_definition,
    get_particle_kind_family,
    get_particle_kind_for_type_id,
    map_particle_kind_to_family,
    particle_type_id_for_kind,
)

from ..common import (
    ParameterBindingToggle,
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QSizePolicy,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QVBoxLayout,
    QWidget,
    Qt,
    nested_get,
    set_combo_value,
)
from ...workflow.domains.beamlets import (
    beamlet_offsets,
    beamlet_pattern_warnings,
    normalize_beamlet_pattern,
)


def _build_default_source(family: str = "H") -> dict[str, object]:
    return {
        "id": "source-1",
        "name": "Main extraction source",
        "particleTypeId": get_default_source_particle_type_id(family),
        "sourceModel": "uniform",
        "particleCount": 50000,
        "currentDensityAm2": 270.0,
        "perpendicularTemperatureEV": 0.0,
        "parallelTemperatureEV": 0.0,
        "axialEnergyEV": 3.0,
        "uniform": {
            "centerMeters": [0.0, 0.0, 0.0],
            "mainDirection": [0.0, 0.0, 1.0],
            "inPlaneReferenceDirection": [1.0, 0.0, 0.0],
            "widthMeters": 0.02,
            "heightMeters": 0.02,
        },
    }


def _ispin(min_: int, max_: int, value: int = 0) -> QSpinBox:
    widget = QSpinBox()
    widget.setRange(min_, max_)
    widget.setValue(value)
    return widget


def _default_particles_spec(family: str = "H") -> dict[str, object]:
    return {
        "types": build_family_particle_types(family),
        "sources": [_build_default_source(family)],
        "plasma": {
            "model": "nsimp",
            "initialPlasmaMaxZMeters": 7.0e-3,
            "positiveIonTemperatureEV": 0.8,
            "plasmaPotentialVolts": 0.0,
        },
    }


def _detect_particles_family(particles: dict[str, object]) -> str:
    type_lookup: dict[str, str] = {}
    raw_types = particles.get("types") if isinstance(particles.get("types"), list) else []
    raw_sources = particles.get("sources") if isinstance(particles.get("sources"), list) else []

    for raw_type in raw_types:
        if not isinstance(raw_type, dict):
            continue
        type_id = str(raw_type.get("id") or "").strip()
        kind = str(raw_type.get("kind") or "").strip()
        if type_id and kind:
            type_lookup[type_id] = kind

    kinds: list[str] = []
    for raw_source in raw_sources:
        if not isinstance(raw_source, dict):
            continue
        particle_type_id = str(raw_source.get("particleTypeId") or "").strip()
        source_kind = type_lookup.get(particle_type_id)
        if not source_kind:
            candidate_kind = str(raw_source.get("kind") or "").strip()
            if candidate_kind:
                source_kind = candidate_kind
        if source_kind:
            kinds.append(source_kind)

    for raw_type in raw_types:
        if isinstance(raw_type, dict):
            candidate_kind = str(raw_type.get("kind") or "").strip()
            if candidate_kind:
                kinds.append(candidate_kind)

    return detect_particle_family_from_kinds(kinds)


def _normalize_source_particle_type_id(
    source: dict[str, object],
    type_lookup: dict[str, str],
    family: str,
) -> str:
    particle_type_id = str(source.get("particleTypeId") or "").strip()
    kind = type_lookup.get(particle_type_id)
    if not kind:
        raw_kind = str(source.get("kind") or "").strip()
        if raw_kind:
            kind = raw_kind
    if not kind and particle_type_id:
        kind = get_particle_kind_for_type_id(particle_type_id)

    if not kind:
        return get_default_source_particle_type_id(family)

    return particle_type_id_for_kind(map_particle_kind_to_family(kind, family))


def _legacy_particles_spec(spec: dict[str, object]) -> dict[str, object]:
    particle_source = nested_get(spec, "particleSource", default={})
    if not isinstance(particle_source, dict):
        particle_source = {}

    kind = str(particle_source.get("species") or "H-")
    family = get_particle_kind_family(kind) or "H"
    default = _default_particles_spec(family)
    default["sources"][0]["particleTypeId"] = particle_type_id_for_kind(map_particle_kind_to_family(kind, family))
    default["sources"][0]["particleCount"] = int(nested_get(spec, "run", "particleCount", default=50000))
    default["sources"][0]["currentDensityAm2"] = float(particle_source.get("currentDensityAm2", 270.0))
    default["sources"][0]["perpendicularTemperatureEV"] = float(
        particle_source.get("perpendicularTemperatureEV", 0.0)
    )
    default["sources"][0]["parallelTemperatureEV"] = float(
        particle_source.get("parallelTemperatureEV", 0.0)
    )
    default["sources"][0]["axialEnergyEV"] = float(particle_source.get("axialEnergyEV", 3.0))
    plasma = default["plasma"]
    model = str(
        nested_get(
            spec,
            "run",
            "solver",
            "plasmaModel",
            default=nested_get(spec, "run", "solver", "shieldModel", default="nsimp"),
        )
    )
    plasma["model"] = model if model in {"nsimp", "shield"} else "nsimp"
    plasma["initialPlasmaMaxZMeters"] = float(
        nested_get(spec, "particles", "plasma", "initialPlasmaMaxZMeters", default=7.0e-3)
    )
    if plasma["model"] == "nsimp":
        plasma["positiveIonTemperatureEV"] = float(
            nested_get(spec, "particles", "plasma", "positiveIonTemperatureEV", default=0.8)
        )
        plasma["plasmaPotentialVolts"] = float(
            nested_get(
                spec,
                "particles",
                "plasma",
                "plasmaPotentialVolts",
                default=nested_get(spec, "particleSource", "plasmaPotentialVolts", default=0.0),
            )
        )
    else:
        plasma.pop("positiveIonTemperatureEV", None)
        plasma.pop("plasmaPotentialVolts", None)
        plasma["tanhWidthEV"] = float(nested_get(spec, "particles", "plasma", "tanhWidthEV", default=1.0))
        plasma["meniscusVoltageVolts"] = float(
            nested_get(spec, "particles", "plasma", "meniscusVoltageVolts", default=0.0)
        )
    return default


def _normalized_particles_spec(spec: dict[str, object]) -> dict[str, object]:
    particles = nested_get(spec, "particles", default=None)
    if not isinstance(particles, dict):
        return _legacy_particles_spec(spec)

    family = _detect_particles_family(particles)
    normalized = _default_particles_spec(family)

    raw_types = particles.get("types") if isinstance(particles.get("types"), list) else []
    type_lookup: dict[str, str] = {}
    for raw_type in raw_types:
        if not isinstance(raw_type, dict):
            continue
        type_id = str(raw_type.get("id") or "").strip()
        kind = str(raw_type.get("kind") or "").strip()
        if type_id and kind:
            type_lookup[type_id] = kind

    raw_sources = particles.get("sources") if isinstance(particles.get("sources"), list) else []
    normalized_sources: list[dict[str, object]] = []
    for raw_source in raw_sources:
        if not isinstance(raw_source, dict):
            continue
        source = deepcopy(raw_source)
        source["particleTypeId"] = _normalize_source_particle_type_id(source, type_lookup, family)
        normalized_sources.append(source)

    if normalized_sources:
        normalized["sources"] = normalized_sources

    plasma = particles.get("plasma") if isinstance(particles.get("plasma"), dict) else None
    if plasma:
        normalized["plasma"] = deepcopy(plasma)
        normalized["plasma"].setdefault("initialPlasmaMaxZMeters", 7.0e-3)

    # Carried through verbatim. Both blocks are optional and absent means off, so copying
    # them only when present keeps a case that never asked for either exactly as it was --
    # and keeps a case that did ask from quietly losing the setting on reload.
    for optional_block in ("halo", "beamlets"):
        block = particles.get(optional_block)
        if isinstance(block, dict):
            normalized[optional_block] = deepcopy(block)

    return normalized


class _ParticleTypesFamilyWidget(QWidget):
    def __init__(self, change_callback=None, parent=None):
        super().__init__(parent)
        self._cb = change_callback
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self._family_combo = QComboBox()
        self._family_combo.addItems(list(PARTICLE_FAMILIES))
        self._family_combo.setMaximumWidth(64)
        self._species_label = QLabel()
        self._species_label.setWordWrap(False)
        layout.addWidget(
            _build_inline_row(
                ("Family", self._family_combo, 0),
                ("Available species", self._species_label, 1),
            )
        )
        self._family_combo.currentTextChanged.connect(self._on_family_changed)
        self._refresh_species_summary()

    def _fire(self) -> None:
        if self._cb:
            self._cb()

    def _refresh_species_summary(self) -> None:
        family = self.get_particle_family()
        names = [
            str(particle_type.get("name") or particle_type.get("kind") or "")
            for particle_type in build_family_particle_types(family)
        ]
        self._species_label.setText(", ".join(name for name in names if name))

    def _on_family_changed(self) -> None:
        self._refresh_species_summary()
        self._fire()

    def set_particles(self, particles: dict[str, object]) -> None:
        family = _detect_particles_family(particles)
        self._family_combo.blockSignals(True)
        set_combo_value(self._family_combo, family)
        self._family_combo.blockSignals(False)
        self._refresh_species_summary()

    def get_particle_family(self) -> str:
        family = self._family_combo.currentText().strip()
        return family if family in PARTICLE_FAMILIES else PARTICLE_FAMILIES[0]

    def get_particle_types(self) -> list[dict[str, object]]:
        return build_family_particle_types(self.get_particle_family())


class _ParticleSourceDetailWidget(QWidget):
    def __init__(self, window, change_callback=None, parent=None):
        super().__init__(parent)
        self._window = window
        self._cb = change_callback
        self._particle_type_lookup: dict[str, dict[str, object]] = {}
        self._editable_widgets = []
        self._source_row = 0
        self._parameter_bindings: list[ParameterBindingToggle] = []
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(8)

        meta_layout = QFormLayout()
        self._id_edit = QLineEdit()
        self._name_edit = QLineEdit()
        self._source_model_combo = QComboBox()
        self._source_model_combo.addItems(["uniform"])
        meta_layout.addRow("Id", self._id_edit)
        meta_layout.addRow("Name", self._name_edit)
        meta_layout.addRow("Source model", self._source_model_combo)
        root.addLayout(meta_layout)

        self._particle_type_combo = QComboBox()
        self._charge_label = QLabel("-1")
        self._mass_label = QLabel("1")
        self._particle_count_spin = QSpinBox()
        self._particle_count_spin.setRange(1, 1000000000)
        self._current_density_spin = self._window._double_spin(0.0, 1.0e9, 3, 1.0)
        self._axial_energy_spin = self._window._double_spin(0.0, 1.0e6, 4, 0.1)
        self._perp_spin = self._window._double_spin(0.0, 1.0e5, 4, 0.1)
        self._par_spin = self._window._double_spin(0.0, 1.0e5, 4, 0.1)
        self._center_x = self._window._double_spin(-1.0e3, 1.0e3, 6, 0.001)
        self._center_y = self._window._double_spin(-1.0e3, 1.0e3, 6, 0.001)
        self._center_z = self._window._double_spin(-1.0e3, 1.0e3, 6, 0.001)
        self._main_x = self._window._double_spin(-1.0, 1.0, 6, 0.1)
        self._main_y = self._window._double_spin(-1.0, 1.0, 6, 0.1)
        self._main_z = self._window._double_spin(-1.0, 1.0, 6, 0.1)
        self._ref_x = self._window._double_spin(-1.0, 1.0, 6, 0.1)
        self._ref_y = self._window._double_spin(-1.0, 1.0, 6, 0.1)
        self._ref_z = self._window._double_spin(-1.0, 1.0, 6, 0.1)
        self._width_spin = self._window._double_spin(0.0, 1.0e3, 6, 0.001)
        self._height_spin = self._window._double_spin(0.0, 1.0e3, 6, 0.001)

        root.addWidget(
            _build_inline_row(
                ("Particle type", self._particle_type_combo, 3),
                ("Charge state", self._charge_label, 1),
                ("Mass [u]", self._mass_label, 1),
            )
        )
        root.addWidget(
            _build_inline_row(
                (
                    "Particle count",
                    self._particle_count_spin,
                    1,
                    self._make_parameter_checkbox(
                        "Particle count",
                        self._particle_count_spin,
                        lambda: self._source_parameter_path("particleCount"),
                    ),
                ),
                (
                    "Current density [A/m^2]",
                    self._current_density_spin,
                    1,
                    self._make_parameter_checkbox(
                        "Current density [A/m^2]",
                        self._current_density_spin,
                        lambda: self._source_parameter_path("currentDensityAm2"),
                    ),
                ),
            )
        )
        root.addWidget(
            _build_inline_row(
                (
                    "Axial energy [eV]",
                    self._axial_energy_spin,
                    1,
                    self._make_parameter_checkbox(
                        "Axial energy [eV]",
                        self._axial_energy_spin,
                        lambda: self._source_parameter_path("axialEnergyEV"),
                    ),
                ),
                (
                    "Perpendicular temperature [eV]",
                    self._perp_spin,
                    1,
                    self._make_parameter_checkbox(
                        "Perpendicular temperature [eV]",
                        self._perp_spin,
                        lambda: self._source_parameter_path("perpendicularTemperatureEV"),
                    ),
                ),
                (
                    "Parallel temperature [eV]",
                    self._par_spin,
                    1,
                    self._make_parameter_checkbox(
                        "Parallel temperature [eV]",
                        self._par_spin,
                        lambda: self._source_parameter_path("parallelTemperatureEV"),
                    ),
                ),
            )
        )
        root.addWidget(
            _build_inline_row(
                (
                    "Center X [m]",
                    self._center_x,
                    1,
                    self._make_parameter_checkbox(
                        "Center X [m]",
                        self._center_x,
                        lambda: self._source_parameter_path("uniform.centerMeters[0]"),
                    ),
                ),
                (
                    "Center Y [m]",
                    self._center_y,
                    1,
                    self._make_parameter_checkbox(
                        "Center Y [m]",
                        self._center_y,
                        lambda: self._source_parameter_path("uniform.centerMeters[1]"),
                    ),
                ),
                (
                    "Center Z [m]",
                    self._center_z,
                    1,
                    self._make_parameter_checkbox(
                        "Center Z [m]",
                        self._center_z,
                        lambda: self._source_parameter_path("uniform.centerMeters[2]"),
                    ),
                ),
            )
        )
        root.addWidget(
            _build_inline_row(
                ("Main direction X", self._main_x, 1),
                ("Main direction Y", self._main_y, 1),
                ("Main direction Z", self._main_z, 1),
            )
        )
        root.addWidget(
            _build_inline_row(
                ("In-plane ref X", self._ref_x, 1),
                ("In-plane ref Y", self._ref_y, 1),
                ("In-plane ref Z", self._ref_z, 1),
            )
        )
        root.addWidget(
            _build_inline_row(
                (
                    "Width [m]",
                    self._width_spin,
                    1,
                    self._make_parameter_checkbox(
                        "Width [m]",
                        self._width_spin,
                        lambda: self._source_parameter_path("uniform.widthMeters"),
                    ),
                ),
                (
                    "Height [m]",
                    self._height_spin,
                    1,
                    self._make_parameter_checkbox(
                        "Height [m]",
                        self._height_spin,
                        lambda: self._source_parameter_path("uniform.heightMeters"),
                    ),
                ),
            )
        )
        root.addStretch(1)

        self._editable_widgets = [
            self._id_edit,
            self._name_edit,
            self._source_model_combo,
            self._particle_count_spin,
            self._current_density_spin,
            self._axial_energy_spin,
            self._perp_spin,
            self._par_spin,
            self._center_x,
            self._center_y,
            self._center_z,
            self._main_x,
            self._main_y,
            self._main_z,
            self._ref_x,
            self._ref_y,
            self._ref_z,
            self._width_spin,
            self._height_spin,
        ]

        for widget in self._editable_widgets:
            _connect_widget_change(widget, self._fire)
        self._particle_type_combo.currentIndexChanged.connect(self._on_particle_type_changed)

    def _source_parameter_path(self, suffix: str) -> str:
        row = self._source_row if self._source_row >= 0 else 0
        return f"particles.sources[{row}].{suffix}"

    def _source_value_reader(self, suffix: str):
        def _reader() -> Any:
            particles_editor = self._window.widgets.get("particles.editor")
            sources = particles_editor.get_sources() if hasattr(particles_editor, "get_sources") else None
            if not isinstance(sources, list) or self._source_row >= len(sources):
                return None
            source = sources[self._source_row]
            current: Any = source
            normalized = suffix.replace("[", ".").replace("]", "")
            for token in normalized.split("."):
                if not token:
                    continue
                if isinstance(current, list):
                    try:
                        current = current[int(token)]
                    except (ValueError, IndexError):
                        return None
                elif isinstance(current, dict):
                    current = current.get(token)
                else:
                    return None
            return current

        return _reader

    def _make_parameter_checkbox(self, label_text: str, widget, path_provider) -> ParameterBindingToggle:
        raw_path = path_provider()
        suffix = raw_path.split("].", 1)[1] if "]." in raw_path else raw_path
        checkbox = self._window._build_parameter_toggle(
            path_provider,
            label_text,
            widget,
            value_reader=self._source_value_reader(suffix),
        )
        self._parameter_bindings.append(checkbox)
        return checkbox

    def _refresh_parameter_checkboxes(self) -> None:
        for checkbox in self._parameter_bindings:
            checkbox.refresh_binding_state()

    def set_source_row(self, row: int) -> None:
        self._source_row = row if row >= 0 else 0
        self._refresh_parameter_checkboxes()

    def _current_particle_type_id(self) -> str:
        current = self._particle_type_combo.currentData()
        if isinstance(current, str) and current:
            return current
        return str(current or "")

    def _refresh_particle_summary(self, particle_type_id: str | None = None) -> None:
        current_id = particle_type_id or self._current_particle_type_id()
        particle_type = self._particle_type_lookup.get(current_id)
        kind = str((particle_type or {}).get("kind") or "H-")
        definition = get_particle_kind_definition(kind)
        self._charge_label.setText(f"{definition.charge_state:g}")
        self._mass_label.setText(f"{definition.mass_u:g}")

    def _on_particle_type_changed(self) -> None:
        self._refresh_particle_summary()
        self._fire()

    def set_available_particle_types(self, particle_types: list[dict[str, object]]) -> None:
        available_particle_types = [deepcopy(item) for item in particle_types if isinstance(item, dict)]
        if not available_particle_types:
            available_particle_types = build_family_particle_types("H")

        current_id = self._current_particle_type_id()
        self._particle_type_lookup = {
            str(item.get("id") or ""): item
            for item in available_particle_types
            if str(item.get("id") or "")
        }

        self._particle_type_combo.blockSignals(True)
        self._particle_type_combo.clear()
        for particle_type in available_particle_types:
            particle_type_id = str(particle_type.get("id") or "")
            if not particle_type_id:
                continue
            particle_type_name = str(particle_type.get("name") or particle_type.get("kind") or particle_type_id)
            self._particle_type_combo.addItem(particle_type_name, particle_type_id)
        self._particle_type_combo.blockSignals(False)

        if current_id not in self._particle_type_lookup:
            current_id = next(iter(self._particle_type_lookup), get_default_source_particle_type_id("H"))
        _set_combo_current_data(self._particle_type_combo, current_id)
        self._refresh_particle_summary(current_id)

    def _fire(self) -> None:
        if self._cb:
            self._cb()

    def set_source(self, source: dict[str, object]) -> None:
        _set_signal_blocked(self._editable_widgets + [self._particle_type_combo], True)
        self._id_edit.setText(str(source.get("id") or "source-1"))
        self._name_edit.setText(str(source.get("name") or ""))
        set_combo_value(self._source_model_combo, str(source.get("sourceModel") or "uniform"))
        self._particle_count_spin.setValue(int(source.get("particleCount") or 1))
        self._current_density_spin.setValue(float(source.get("currentDensityAm2") or 0.0))
        self._perp_spin.setValue(float(source.get("perpendicularTemperatureEV") or 0.0))
        self._par_spin.setValue(float(source.get("parallelTemperatureEV") or 0.0))
        self._axial_energy_spin.setValue(float(source.get("axialEnergyEV") or 0.0))
        uniform = source.get("uniform") if isinstance(source.get("uniform"), dict) else {}
        center = uniform.get("centerMeters") if isinstance(uniform, dict) else None
        main = uniform.get("mainDirection") if isinstance(uniform, dict) else None
        reference = uniform.get("inPlaneReferenceDirection") if isinstance(uniform, dict) else None
        _set_spin_triplet([self._center_x, self._center_y, self._center_z], center, [0.0, 0.0, 0.0])
        _set_spin_triplet([self._main_x, self._main_y, self._main_z], main, [0.0, 0.0, 1.0])
        _set_spin_triplet([self._ref_x, self._ref_y, self._ref_z], reference, [1.0, 0.0, 0.0])
        self._width_spin.setValue(float(uniform.get("widthMeters") or 0.02))
        self._height_spin.setValue(float(uniform.get("heightMeters") or 0.02))

        particle_type_id = str(source.get("particleTypeId") or "").strip()
        if particle_type_id not in self._particle_type_lookup:
            particle_type_id = next(iter(self._particle_type_lookup), get_default_source_particle_type_id("H"))
        _set_combo_current_data(self._particle_type_combo, particle_type_id)
        _set_signal_blocked(self._editable_widgets + [self._particle_type_combo], False)
        self._refresh_particle_summary(particle_type_id)
        self._refresh_parameter_checkboxes()

    def get_source(self) -> dict[str, object]:
        particle_type_id = self._current_particle_type_id() or next(
            iter(self._particle_type_lookup),
            get_default_source_particle_type_id("H"),
        )
        source_id = self._id_edit.text().strip() or "source-1"
        return {
            "id": source_id,
            "name": self._name_edit.text().strip() or source_id,
            "particleTypeId": particle_type_id,
            "sourceModel": self._source_model_combo.currentText().strip() or "uniform",
            "particleCount": int(self._particle_count_spin.value()),
            "currentDensityAm2": float(self._current_density_spin.value()),
            "perpendicularTemperatureEV": float(self._perp_spin.value()),
            "parallelTemperatureEV": float(self._par_spin.value()),
            "axialEnergyEV": float(self._axial_energy_spin.value()),
            "uniform": {
                "centerMeters": [
                    float(self._center_x.value()),
                    float(self._center_y.value()),
                    float(self._center_z.value()),
                ],
                "mainDirection": [
                    float(self._main_x.value()),
                    float(self._main_y.value()),
                    float(self._main_z.value()),
                ],
                "inPlaneReferenceDirection": [
                    float(self._ref_x.value()),
                    float(self._ref_y.value()),
                    float(self._ref_z.value()),
                ],
                "widthMeters": float(self._width_spin.value()),
                "heightMeters": float(self._height_spin.value()),
            },
        }


class _ParticleSourcesEditorWidget(QWidget):
    def __init__(self, window, change_callback=None, parent=None):
        super().__init__(parent)
        self._window = window
        self._cb = change_callback
        self._prev_row = -1
        self._available_particle_types: list[dict[str, object]] = []
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        root.addWidget(splitter, 1)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 4, 0)
        left_layout.setSpacing(4)
        buttons = QHBoxLayout()
        add_btn = QPushButton("+ Add")
        rem_btn = QPushButton("− Remove")
        add_btn.clicked.connect(self._add_source)
        rem_btn.clicked.connect(self._remove_source)
        buttons.addWidget(add_btn)
        buttons.addWidget(rem_btn)
        left_layout.addLayout(buttons)
        self._list = QListWidget()
        self._list.setMinimumWidth(84)
        self._list.setMaximumWidth(144)
        self._list.currentRowChanged.connect(self._on_row_changed)
        left_layout.addWidget(self._list)

        self._detail = _ParticleSourceDetailWidget(self._window, change_callback=self._cb)
        detail_scroll = QScrollArea()
        detail_scroll.setWidget(self._detail)
        detail_scroll.setWidgetResizable(True)
        detail_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        splitter.addWidget(left)
        splitter.addWidget(detail_scroll)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([120, 390])
        self.setMinimumHeight(260)

    def _default_family(self) -> str:
        kinds = [
            str(particle_type.get("kind") or "")
            for particle_type in self._available_particle_types
            if isinstance(particle_type, dict)
        ]
        return detect_particle_family_from_kinds(kinds)

    def _load_current(self) -> None:
        row = self._list.currentRow()
        self._detail.set_source_row(row)
        if row < 0:
            self._detail.set_source(_build_default_source(self._default_family()))
            return
        item = self._list.item(row)
        self._detail.set_source(item.data(Qt.ItemDataRole.UserRole) or {})

    def refresh_particle_types(self, particle_types: list[dict[str, object]]) -> None:
        self._flush_current()
        self._available_particle_types = [deepcopy(item) for item in particle_types if isinstance(item, dict)]
        self._detail.set_available_particle_types(self._available_particle_types)
        self._load_current()

    def remap_particle_family(self, old_family: str, new_family: str) -> None:
        if old_family == new_family:
            return

        self._flush_current()
        old_types = {
            str(particle_type.get("id") or ""): str(particle_type.get("kind") or "")
            for particle_type in build_family_particle_types(old_family)
        }
        for index in range(self._list.count()):
            item = self._list.item(index)
            source = deepcopy(item.data(Qt.ItemDataRole.UserRole) or {})
            old_id = str(source.get("particleTypeId") or "").strip()
            kind = old_types.get(old_id) or get_particle_kind_for_type_id(old_id)
            if kind:
                source["particleTypeId"] = particle_type_id_for_kind(map_particle_kind_to_family(kind, new_family))
                item.setData(Qt.ItemDataRole.UserRole, source)

    def _fire(self) -> None:
        if self._cb:
            self._cb()

    def _flush_current(self) -> None:
        if 0 <= self._prev_row < self._list.count():
            item = self._list.item(self._prev_row)
            updated = self._detail.get_source()
            item.setData(Qt.ItemDataRole.UserRole, updated)
            item.setText(str(updated.get("name") or updated.get("id") or "unnamed"))

    def _on_row_changed(self, row: int) -> None:
        self._flush_current()
        self._prev_row = row
        self._load_current()

    def _add_source(self) -> None:
        family = self._default_family()
        new_source = deepcopy(_build_default_source(family))
        new_source["id"] = f"source-{self._list.count() + 1}"
        new_source["name"] = f"Source {self._list.count() + 1}"
        item = QListWidgetItem(str(new_source["name"]))
        item.setData(Qt.ItemDataRole.UserRole, new_source)
        self._list.addItem(item)
        self._list.setCurrentItem(item)
        self._fire()

    def _remove_source(self) -> None:
        row = self._list.currentRow()
        if row < 0:
            return
        self._list.takeItem(row)
        self._prev_row = -1
        self._load_current()
        self._fire()

    def set_sources(self, sources: list[dict[str, object]]) -> None:
        self._list.blockSignals(True)
        self._list.clear()
        for source in sources:
            item = QListWidgetItem(str(source.get("name") or source.get("id") or "unnamed"))
            item.setData(Qt.ItemDataRole.UserRole, deepcopy(source))
            self._list.addItem(item)
        self._list.blockSignals(False)
        self._prev_row = -1
        if self._list.count() > 0:
            self._list.setCurrentRow(0)
        else:
            self._load_current()

    def get_sources(self) -> list[dict[str, object]]:
        self._flush_current()
        sources: list[dict[str, object]] = []
        for index in range(self._list.count()):
            item = self._list.item(index)
            data = item.data(Qt.ItemDataRole.UserRole)
            if data:
                sources.append(deepcopy(data))
        return sources


class _ParticlesWorkspaceWidget(QWidget):
    def __init__(self, window, parent=None):
        super().__init__(parent)
        self._window = window
        self._particle_family = PARTICLE_FAMILIES[0]
        self._build_ui()

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(10)

        types_box = QGroupBox("Particle types")
        types_layout = QVBoxLayout(types_box)
        self._types_editor = _ParticleTypesFamilyWidget(change_callback=self._on_particle_family_changed)
        types_layout.addWidget(self._types_editor)
        root.addWidget(types_box)

        sources_box = QGroupBox("Sources")
        sources_box.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        sources_layout = QVBoxLayout(sources_box)
        self._sources_editor = _ParticleSourcesEditorWidget(self._window, change_callback=self._notify_change)
        sources_layout.addWidget(self._sources_editor, 1)
        sources_layout.addWidget(self._build_beamlet_box())

        plasma_box = QGroupBox("Plasma model")
        plasma_box.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        plasma_layout = QVBoxLayout(plasma_box)
        plasma_layout.setContentsMargins(8, 8, 8, 8)
        plasma_layout.setSpacing(8)
        self._plasma_model = QComboBox()
        self._plasma_model.addItems(["nsimp", "shield"])
        self._plasma_model.setMaximumWidth(120)
        self._initial_plasma_max_z = self._window._double_spin(-10.0, 10.0, 6, 0.001)
        plasma_layout.addWidget(
            _build_inline_row(
                ("Model", self._plasma_model, 0),
                (
                    "Initial plasma max z [m]",
                    self._initial_plasma_max_z,
                    1,
                    self._make_plasma_parameter_checkbox(
                        "Initial plasma max z [m]",
                        self._initial_plasma_max_z,
                        "particles.plasma.initialPlasmaMaxZMeters",
                    ),
                ),
            )
        )

        self._positive_ion_temp = self._window._double_spin(0.0, 1.0e4, 4, 0.1)
        self._plasma_potential = self._window._double_spin(-1.0e6, 1.0e6, 4, 1.0)
        self._tanh_width = self._window._double_spin(0.0, 1.0e4, 4, 0.1)
        self._meniscus_voltage = self._window._double_spin(-1.0e6, 1.0e6, 4, 1.0)
        self._nsimp_fields = _build_form_panel(
            (
                "Positive ion temperature [eV]",
                self._positive_ion_temp,
                self._make_plasma_parameter_checkbox(
                    "Positive ion temperature [eV]",
                    self._positive_ion_temp,
                    "particles.plasma.positiveIonTemperatureEV",
                ),
            ),
            (
                "Plasma potential [V]",
                self._plasma_potential,
                self._make_plasma_parameter_checkbox(
                    "Plasma potential [V]",
                    self._plasma_potential,
                    "particles.plasma.plasmaPotentialVolts",
                ),
            ),
        )
        self._shield_fields = _build_form_panel(
            (
                "Tanh width [eV]",
                self._tanh_width,
                self._make_plasma_parameter_checkbox(
                    "Tanh width [eV]",
                    self._tanh_width,
                    "particles.plasma.tanhWidthEV",
                ),
            ),
            (
                "Meniscus voltage [V]",
                self._meniscus_voltage,
                self._make_plasma_parameter_checkbox(
                    "Meniscus voltage [V]",
                    self._meniscus_voltage,
                    "particles.plasma.meniscusVoltageVolts",
                ),
            ),
        )
        plasma_layout.addWidget(self._nsimp_fields)
        plasma_layout.addWidget(self._shield_fields)
        plasma_layout.addStretch(1)

        # --- Beam halo -------------------------------------------------------------
        # Pimazzoni et al., Fus. Eng. Des. 192 (2023) 113621, Sec. 2: halo particles are born
        # on the downstream chamfer of the plasma grid and aimed along its normal. Off unless
        # asked for, so an existing case is unaffected.
        halo_box = QGroupBox("Beam halo")
        halo_layout = QVBoxLayout(halo_box)
        halo_layout.setContentsMargins(8, 8, 8, 8)
        halo_layout.setSpacing(8)

        self._halo_enabled = QCheckBox("Emit beamlet halo from the plasma-grid chamfer")
        self._halo_enabled.setToolTip(
            "Adds halo particles born on the downstream chamfer of the plasma grid, aimed "
            "along the local surface normal, following Pimazzoni et al., Fus. Eng. Des. 192 "
            "(2023) 113621. Leave unchecked and no halo is emitted."
        )
        halo_layout.addWidget(self._halo_enabled)

        self._halo_fraction = self._window._double_spin(0.1, 50.0, 2, 0.5)
        self._halo_fraction.setToolTip(
            "Halo current as a percentage of the total core source current. The reference "
            "work fixes this at 8 %, and states that holding it independent of beamlet "
            "perveance is the largest single uncertainty in its results. The halo "
            "macroparticle count is the same fraction of the core count, so each halo "
            "particle carries the same charge as a core one and the extra run time is just "
            "this percentage."
        )
        self._halo_energy = self._window._double_spin(0.01, 1000.0, 3, 0.5)
        self._halo_energy.setToolTip(
            "Emission energy normal to the chamfer. The reference uses 4 eV and notes the "
            "initial velocity is of minor relevance, because the field at the emission point "
            "focuses these particles strongly."
        )
        self._halo_fields = _build_form_panel(
            ("Halo current [% of core]", self._halo_fraction, None),
            ("Halo axial energy [eV]", self._halo_energy, None),
        )
        halo_layout.addWidget(self._halo_fields)

        self._halo_note = QLabel()
        self._halo_note.setWordWrap(True)
        self._halo_note.setStyleSheet("color: #666; font-size: 11px;")
        halo_layout.addWidget(self._halo_note)
        halo_layout.addStretch(1)

        _connect_widget_change(self._halo_enabled, self._on_halo_toggled)
        _connect_widget_change(self._halo_fraction, self._on_halo_toggled)
        _connect_widget_change(self._halo_energy, self._notify_change)

        lower_panel = QWidget()
        lower_row = QHBoxLayout(lower_panel)
        lower_row.setContentsMargins(0, 0, 0, 0)
        lower_row.setSpacing(10)
        lower_row.addWidget(sources_box, 3)
        lower_row.addWidget(plasma_box, 2)
        lower_row.addWidget(halo_box, 2)
        root.addWidget(lower_panel, 1)

        _connect_widget_change(self._plasma_model, self._on_plasma_model_changed)
        _connect_widget_change(self._initial_plasma_max_z, self._notify_change)
        _connect_widget_change(self._positive_ion_temp, self._notify_change)
        _connect_widget_change(self._plasma_potential, self._notify_change)
        _connect_widget_change(self._tanh_width, self._notify_change)
        _connect_widget_change(self._meniscus_voltage, self._notify_change)

        self._sources_editor.refresh_particle_types(self._types_editor.get_particle_types())
        self._on_plasma_model_changed()
        self._on_beamlets_changed()

    # --- Multi-beamlet extraction ---------------------------------------------------
    # The grid solids already tile their bore over an aperture pattern. This tiles the beam
    # over the matching lattice, so a patterned geometry is fed one beamlet per aperture
    # instead of one hand-placed source per aperture. Only the source centres move: every
    # beamlet keeps the current density, particle count, energy and temperatures of the
    # source it was replicated from, so the physics of a single beamlet is unchanged.

    def _build_beamlet_box(self) -> QGroupBox:
        box = QGroupBox("Multi-beamlet extraction")
        layout = QVBoxLayout(box)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        self._beamlets_enabled = QCheckBox("Extract one beamlet per grid aperture")
        self._beamlets_enabled.setToolTip(
            "Replicates every source over a lattice of beamlets. Leave unchecked and a "
            "single beamlet is extracted, exactly as before.\n\n"
            "Each beamlet carries the full per-source current and particle count, so a "
            "beamlet in an N-beamlet run is statistically identical to the same beamlet run "
            "on its own -- but the total current, particle count, memory and run time all "
            "scale with the number of beamlets."
        )
        layout.addWidget(self._beamlets_enabled)

        self._beamlets_layout_cb = QComboBox()
        self._beamlets_layout_cb.addItems(["rectangular-grid", "staggered-grid"])
        self._beamlets_count_x = _ispin(1, 199, 5)
        self._beamlets_count_y = _ispin(1, 199, 5)
        self._beamlets_pitch_x = self._window._double_spin(0.0, 1.0, 6, 0.001)
        self._beamlets_pitch_y = self._window._double_spin(0.0, 1.0, 6, 0.001)
        self._beamlets_x_offset = self._window._double_spin(-1.0, 1.0, 6, 0.001)
        self._beamlets_y_offset = self._window._double_spin(-1.0, 1.0, 6, 0.001)
        self._beamlets_row_shift = self._window._double_spin(-1.0, 1.0, 6, 0.001)
        self._beamlets_pitch_x.setValue(0.019)
        self._beamlets_pitch_y.setValue(0.021)

        self._beamlets_row_shift.setToolTip(
            "Alternating row offset, following the grid geometry: even rows are shifted by "
            "-this and odd rows by +this, so the row-to-row stagger is TWICE this value. "
            "For a half-pitch stagger enter a quarter of the pitch."
        )
        self._beamlets_row_shift_row = _build_form_panel(
            ("Row shift X [m]", self._beamlets_row_shift, None),
        )
        self._beamlets_fields = QWidget()
        fields_layout = QVBoxLayout(self._beamlets_fields)
        fields_layout.setContentsMargins(0, 0, 0, 0)
        fields_layout.setSpacing(4)
        fields_layout.addWidget(
            _build_inline_row(
                ("Layout", self._beamlets_layout_cb, 1),
                ("Count X", self._beamlets_count_x, 0),
                ("Count Y", self._beamlets_count_y, 0),
            )
        )
        fields_layout.addWidget(
            _build_inline_row(
                ("Pitch X [m]", self._beamlets_pitch_x, 1),
                ("Pitch Y [m]", self._beamlets_pitch_y, 1),
            )
        )
        fields_layout.addWidget(
            _build_inline_row(
                ("X offset [m]", self._beamlets_x_offset, 1),
                ("Y offset [m]", self._beamlets_y_offset, 1),
            )
        )
        fields_layout.addWidget(self._beamlets_row_shift_row)

        sync_button = QPushButton("Sync from grid aperture pattern")
        sync_button.setToolTip(
            "Copies the aperture pattern of the plasma grid into these fields, so every "
            "beamlet is emitted into a real aperture. Nothing else keeps the two in step: "
            "a mismatch fires beamlets at solid metal and the run will still complete."
        )
        sync_button.clicked.connect(self._sync_beamlets_from_geometry)
        fields_layout.addWidget(sync_button)
        layout.addWidget(self._beamlets_fields)

        self._beamlets_note = QLabel()
        self._beamlets_note.setWordWrap(True)
        self._beamlets_note.setStyleSheet("color: #666; font-size: 11px;")
        layout.addWidget(self._beamlets_note)

        _connect_widget_change(self._beamlets_enabled, self._on_beamlets_changed)
        for widget in (
            self._beamlets_layout_cb,
            self._beamlets_count_x,
            self._beamlets_count_y,
            self._beamlets_pitch_x,
            self._beamlets_pitch_y,
            self._beamlets_x_offset,
            self._beamlets_y_offset,
            self._beamlets_row_shift,
        ):
            _connect_widget_change(widget, self._on_beamlets_changed)

        return box

    def _beamlet_pattern(self) -> dict[str, object]:
        """The block as the authoring spec carries it, whether enabled or not."""
        layout_name = self._beamlets_layout_cb.currentText().strip() or "rectangular-grid"
        pattern: dict[str, object] = {
            "enabled": bool(self._beamlets_enabled.isChecked()),
            "layout": layout_name,
            "countX": int(self._beamlets_count_x.value()),
            "countY": int(self._beamlets_count_y.value()),
            "pitchXMeters": float(self._beamlets_pitch_x.value()),
            "pitchYMeters": float(self._beamlets_pitch_y.value()),
            "xOffsetMeters": float(self._beamlets_x_offset.value()),
            "yOffsetMeters": float(self._beamlets_y_offset.value()),
        }
        if layout_name == "staggered-grid":
            pattern["rowShiftXMeters"] = float(self._beamlets_row_shift.value())
        return pattern

    def _plasma_grid_aperture_pattern(self) -> dict[str, object] | None:
        """The plasma grid's aperture pattern, found the way the runtime finds it.

        Name first, then the solid that starts furthest upstream -- the same rule as
        findPlasmaGrid in ParticleManager.cpp, so the GUI and the run agree on which solid
        the beamlets belong to.
        """
        editor = self._window.widgets.get("geometry.solidsEditor")
        if editor is None or not hasattr(editor, "get_solids"):
            return None
        try:
            solids = editor.get_solids()
        except Exception:                                       # noqa: BLE001
            return None

        upstream: tuple[float, dict] | None = None
        for solid in solids:
            if not isinstance(solid, dict):
                continue
            z_profile = solid.get("zProfileMeters")
            if not isinstance(z_profile, list) or len(z_profile) < 2:
                continue
            name = str(solid.get("name") or "").strip().upper()
            if name in {"PG", "PLASMA GRID", "PLASMAGRID"} or solid.get("role") == "plasmaGrid":
                return solid.get("aperturePattern") if isinstance(solid.get("aperturePattern"), dict) else {}
            try:
                z_front = float(z_profile[0])
            except (TypeError, ValueError):
                continue
            if upstream is None or z_front < upstream[0]:
                upstream = (z_front, solid)

        if upstream is None:
            return None
        pattern = upstream[1].get("aperturePattern")
        return pattern if isinstance(pattern, dict) else {}

    def _sync_beamlets_from_geometry(self) -> None:
        pattern = self._plasma_grid_aperture_pattern()
        if pattern is None:
            self._beamlets_note.setText(
                "No grid solid found to sync from. Define the geometry first."
            )
            return
        if str(pattern.get("layout", "single")) == "single":
            self._beamlets_note.setText(
                "The plasma grid has a single aperture, so there is nothing to tile. Give it "
                "a rectangular or staggered aperture pattern in the Geometry section first."
            )
            return

        widgets = (
            self._beamlets_layout_cb,
            self._beamlets_count_x,
            self._beamlets_count_y,
            self._beamlets_pitch_x,
            self._beamlets_pitch_y,
            self._beamlets_x_offset,
            self._beamlets_y_offset,
            self._beamlets_row_shift,
        )
        _set_signal_blocked(widgets, True)
        try:
            set_combo_value(self._beamlets_layout_cb, str(pattern.get("layout")))
            self._beamlets_count_x.setValue(int(pattern.get("countX", 1) or 1))
            self._beamlets_count_y.setValue(int(pattern.get("countY", 1) or 1))
            self._beamlets_pitch_x.setValue(float(pattern.get("pitchXMeters", 0.0) or 0.0))
            self._beamlets_pitch_y.setValue(float(pattern.get("pitchYMeters", 0.0) or 0.0))
            self._beamlets_x_offset.setValue(float(pattern.get("xOffsetMeters", 0.0) or 0.0))
            self._beamlets_y_offset.setValue(float(pattern.get("yOffsetMeters", 0.0) or 0.0))
            self._beamlets_row_shift.setValue(float(pattern.get("rowShiftXMeters", 0.0) or 0.0))
        finally:
            _set_signal_blocked(widgets, False)
        self._on_beamlets_changed()

    def _on_beamlets_changed(self) -> None:
        enabled = self._beamlets_enabled.isChecked()
        self._beamlets_fields.setEnabled(enabled)
        staggered = self._beamlets_layout_cb.currentText().strip() == "staggered-grid"
        self._beamlets_row_shift_row.setVisible(staggered)

        if not enabled:
            self._beamlets_note.setText("Single beamlet on the machine axis.")
            self._on_halo_toggled()
            return

        try:
            pattern = normalize_beamlet_pattern(self._beamlet_pattern())
            offsets = beamlet_offsets(pattern) if pattern else []
            warnings = beamlet_pattern_warnings(pattern, self._plasma_grid_aperture_pattern())
        except Exception as error:                              # noqa: BLE001
            self._beamlets_note.setText(str(error))
            self._on_halo_toggled()
            return

        try:
            per_beamlet = sum(
                int(source.get("particleCount") or 0)
                for source in self._sources_editor.get_sources()
            )
        except Exception:                                       # noqa: BLE001
            per_beamlet = 0

        xs = [offset[0] for offset in offsets] or [0.0]
        ys = [offset[1] for offset in offsets] or [0.0]
        lines = [
            "{count} beamlets, {span_x:.1f} x {span_y:.1f} mm across. "
            "{total:,} core macroparticles in total ({per:,} per beamlet).".format(
                count=len(offsets),
                span_x=1.0e3 * (max(xs) - min(xs)),
                span_y=1.0e3 * (max(ys) - min(ys)),
                total=per_beamlet * len(offsets),
                per=per_beamlet,
            ),
            "Field mesh and memory scale with the beamlet count, not the particle count: "
            "run this on Slurm, not the login node.",
        ]
        lines.extend(warnings)
        mismatch = self._beamlet_geometry_mismatch()
        if mismatch:
            lines.append(mismatch)
        domain = self._beamlet_domain_warning(offsets)
        if domain:
            lines.append(domain)
        self._beamlets_note.setText(" ".join(lines))
        self._on_halo_toggled()

    def _beamlet_domain_warning(self, offsets: list[tuple[float, float]]) -> str:
        """Warn when the lattice will not fit the field domain.

        The runtime rejects a source whose footprint leaves geometry.domain, so without this
        the first sign of a too-small domain is a failed run. Reported here instead, while
        the domain is still one field away.
        """
        try:
            domain_x = float(self._window.widgets["geometry.domain.x"].value())
            domain_y = float(self._window.widgets["geometry.domain.y"].value())
        except Exception:                                       # noqa: BLE001
            return ""

        half_width = 0.0
        half_height = 0.0
        try:
            for source in self._sources_editor.get_sources():
                uniform = source.get("uniform")
                if isinstance(uniform, dict):
                    half_width = max(half_width, 0.5 * float(uniform.get("widthMeters") or 0.0))
                    half_height = max(half_height, 0.5 * float(uniform.get("heightMeters") or 0.0))
        except Exception:                                       # noqa: BLE001
            return ""

        needed_x = 2.0 * (max(abs(offset[0]) for offset in offsets) + half_width)
        needed_y = 2.0 * (max(abs(offset[1]) for offset in offsets) + half_height)
        if needed_x <= domain_x + 1.0e-12 and needed_y <= domain_y + 1.0e-12:
            return ""
        return (
            "The lattice needs a domain of at least {nx:.1f} x {ny:.1f} mm and the geometry "
            "domain is {dx:.1f} x {dy:.1f} mm, so the run will be rejected. Enlarge "
            "geometry.domain.".format(
                nx=1.0e3 * needed_x, ny=1.0e3 * needed_y,
                dx=1.0e3 * domain_x, dy=1.0e3 * domain_y,
            )
        )

    def _beamlet_geometry_mismatch(self) -> str:
        """Warn when the beamlet lattice and the grid aperture pattern have drifted apart."""
        pattern = self._plasma_grid_aperture_pattern()
        if not pattern or str(pattern.get("layout", "single")) == "single":
            return (
                "The plasma grid has no aperture pattern, so these beamlets are being fired "
                "into a single-aperture grid."
            )

        differences: list[str] = []
        for label, key, value in (
            ("layout", "layout", self._beamlets_layout_cb.currentText().strip()),
            ("count X", "countX", int(self._beamlets_count_x.value())),
            ("count Y", "countY", int(self._beamlets_count_y.value())),
            ("pitch X", "pitchXMeters", float(self._beamlets_pitch_x.value())),
            ("pitch Y", "pitchYMeters", float(self._beamlets_pitch_y.value())),
        ):
            grid_value = pattern.get(key)
            if isinstance(value, float):
                if grid_value is None or abs(float(grid_value) - value) > 1.0e-9:
                    differences.append(label)
            elif str(grid_value) != str(value):
                differences.append(label)

        if not differences:
            return ""
        return (
            "Does not match the grid aperture pattern ({}): some beamlets will be fired at "
            "solid metal. Use 'Sync from grid aperture pattern'.".format(", ".join(differences))
        )

    def _notify_change(self) -> None:
        self._window.schedule_preview_refresh()

    def _make_plasma_parameter_checkbox(self, label_text: str, widget, path: str) -> ParameterBindingToggle:
        return self._window._build_parameter_toggle(lambda p=path: p, label_text, widget)

    def _on_particle_family_changed(self) -> None:
        new_family = self._types_editor.get_particle_family()
        if new_family != self._particle_family:
            self._sources_editor.remap_particle_family(self._particle_family, new_family)
            self._particle_family = new_family
        self._sources_editor.refresh_particle_types(self._types_editor.get_particle_types())
        self._notify_change()

    def _on_halo_toggled(self) -> None:
        """Grey out the halo inputs when disabled and show what it will cost."""
        enabled = self._halo_enabled.isChecked()
        self._halo_fields.setEnabled(enabled)
        if not enabled:
            self._halo_note.setText("No halo particles will be emitted.")
        else:
            try:
                core = sum(int(src.get("particleCount") or 0)
                           for src in self._sources_editor.get_sources())
            except Exception:                                   # noqa: BLE001
                core = 0
            pct = float(self._halo_fraction.value())
            # Per beamlet: each beamlet gets its own halo ring on its own aperture chamfer,
            # carrying its own share of the current (ParticleManager.cpp, addBeamHalo).
            beamlets = self._beamlet_multiplier()
            extra = int(round(0.01 * pct * core)) * beamlets
            core_total = core * beamlets
            per_beamlet = (
                " One ring per beamlet, on its own aperture chamfer." if beamlets > 1 else ""
            )
            self._halo_note.setText(
                "Adds about {:,} macroparticles to the {:,} core particles "
                "(+{:.1f}% run time and memory).{}".format(extra, core_total, pct, per_beamlet)
            )
        self._notify_change()

    def _beamlet_multiplier(self) -> int:
        """How many beamlets the current settings produce; 1 when disabled or invalid."""
        if not getattr(self, "_beamlets_enabled", None) or not self._beamlets_enabled.isChecked():
            return 1
        try:
            pattern = normalize_beamlet_pattern(self._beamlet_pattern())
            return max(1, len(beamlet_offsets(pattern))) if pattern else 1
        except Exception:                                       # noqa: BLE001
            return 1

    def _on_plasma_model_changed(self) -> None:
        model = self._plasma_model.currentText().strip() or "nsimp"
        is_nsimp = model == "nsimp"
        self._nsimp_fields.setVisible(is_nsimp)
        self._shield_fields.setVisible(not is_nsimp)
        self._notify_change()

    def set_particles(self, particles: dict[str, object]) -> None:
        self._types_editor.set_particles(particles)
        self._particle_family = self._types_editor.get_particle_family()
        self._sources_editor.refresh_particle_types(self._types_editor.get_particle_types())
        sources = particles.get("sources") if isinstance(particles.get("sources"), list) else []
        plasma = particles.get("plasma") if isinstance(particles.get("plasma"), dict) else {}
        self._sources_editor.set_sources([deepcopy(item) for item in sources])
        set_combo_value(self._plasma_model, str(plasma.get("model") or "nsimp"))
        initial_plasma_max_z = plasma.get("initialPlasmaMaxZMeters")
        self._initial_plasma_max_z.setValue(float(initial_plasma_max_z if initial_plasma_max_z is not None else 7.0e-3))
        self._positive_ion_temp.setValue(float(plasma.get("positiveIonTemperatureEV") or 0.8))
        self._plasma_potential.setValue(float(plasma.get("plasmaPotentialVolts") or 0.0))
        self._tanh_width.setValue(float(plasma.get("tanhWidthEV") or 1.0))
        self._meniscus_voltage.setValue(float(plasma.get("meniscusVoltageVolts") or 0.0))

        halo = particles.get("halo") if isinstance(particles.get("halo"), dict) else {}
        self._halo_enabled.setChecked(bool(halo.get("enabled", False)))
        self._halo_fraction.setValue(float(halo.get("currentFractionPercent") or 8.0))
        self._halo_energy.setValue(float(halo.get("axialEnergyEV") or 4.0))

        beamlets = particles.get("beamlets") if isinstance(particles.get("beamlets"), dict) else {}
        self._beamlets_enabled.setChecked(bool(beamlets.get("enabled", False)))
        set_combo_value(self._beamlets_layout_cb, str(beamlets.get("layout") or "rectangular-grid"))
        self._beamlets_count_x.setValue(int(beamlets.get("countX") or 5))
        self._beamlets_count_y.setValue(int(beamlets.get("countY") or 5))
        self._beamlets_pitch_x.setValue(float(beamlets.get("pitchXMeters") or 0.019))
        self._beamlets_pitch_y.setValue(float(beamlets.get("pitchYMeters") or 0.021))
        self._beamlets_x_offset.setValue(float(beamlets.get("xOffsetMeters") or 0.0))
        self._beamlets_y_offset.setValue(float(beamlets.get("yOffsetMeters") or 0.0))
        self._beamlets_row_shift.setValue(float(beamlets.get("rowShiftXMeters") or 0.0))

        self._on_halo_toggled()
        self._on_plasma_model_changed()
        self._on_beamlets_changed()

    def get_particle_family(self) -> str:
        return self._types_editor.get_particle_family()

    def get_particle_types(self) -> list[dict[str, object]]:
        return self._types_editor.get_particle_types()

    def get_sources(self) -> list[dict[str, object]]:
        self._sources_editor.refresh_particle_types(self._types_editor.get_particle_types())
        return self._sources_editor.get_sources()

    def get_particles(self) -> dict[str, object]:
        particles = {
            "types": self.get_particle_types(),
            "sources": self.get_sources(),
            "plasma": {
                "model": self._plasma_model.currentText().strip() or "nsimp",
                "initialPlasmaMaxZMeters": float(self._initial_plasma_max_z.value()),
            },
        }
        if particles["plasma"]["model"] == "nsimp":
            particles["plasma"]["positiveIonTemperatureEV"] = float(self._positive_ion_temp.value())
            particles["plasma"]["plasmaPotentialVolts"] = float(self._plasma_potential.value())
        else:
            particles["plasma"]["tanhWidthEV"] = float(self._tanh_width.value())
            particles["plasma"]["meniscusVoltageVolts"] = float(self._meniscus_voltage.value())
        # Written only when enabled: an absent block means no halo, all the way down to the
        # runtime, so a case cannot acquire one by defaulting.
        if self._halo_enabled.isChecked():
            particles["halo"] = {
                "enabled": True,
                "currentFractionPercent": float(self._halo_fraction.value()),
                "axialEnergyEV": float(self._halo_energy.value()),
            }
        # Same rule as the halo: written only when enabled, so a single-beamlet case stays
        # exactly the case it was and cannot acquire a lattice by defaulting.
        if self._beamlets_enabled.isChecked():
            particles["beamlets"] = self._beamlet_pattern()
        return particles


def _connect_widget_change(widget, callback) -> None:
    if isinstance(widget, QLineEdit):
        widget.textChanged.connect(callback)
    elif isinstance(widget, QComboBox):
        widget.currentTextChanged.connect(callback)
    elif isinstance(widget, QCheckBox):
        # Checkboxes have no valueChanged. Use toggled(bool) rather than stateChanged: it is
        # the signal on QAbstractButton and is stable across Qt versions, whereas
        # stateChanged is deprecated from Qt 6.7 in favour of checkStateChanged.
        widget.toggled.connect(callback)
    elif hasattr(widget, "valueChanged"):
        widget.valueChanged.connect(callback)
    else:
        # Fail loudly and usefully rather than with a bare AttributeError from a signal that
        # happens not to exist on whatever was passed in.
        raise TypeError(
            "_connect_widget_change does not know how to observe %s; add a branch for it"
            % type(widget).__name__
        )


def _set_signal_blocked(widgets, blocked: bool) -> None:
    for widget in widgets:
        widget.blockSignals(blocked)


def _set_spin_triplet(spins, values, fallback) -> None:
    vector = values if isinstance(values, list) and len(values) == 3 else fallback
    for spin, value in zip(spins, vector):
        spin.setValue(float(value))


def _set_combo_current_data(combo: QComboBox, value: str) -> None:
    for index in range(combo.count()):
        if str(combo.itemData(index) or "") == value:
            combo.setCurrentIndex(index)
            return
    if combo.count() > 0:
        combo.setCurrentIndex(0)


def _build_inline_control(label_text: str, control, stretch: int = 0, extra_widget=None) -> tuple[QWidget, int]:
    container = QWidget()
    layout = QVBoxLayout(container)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(2)
    if extra_widget is None:
        caption = QLabel(label_text)
        caption.setWordWrap(True)
        layout.addWidget(caption)
    else:
        caption_row = QWidget()
        caption_layout = QHBoxLayout(caption_row)
        caption_layout.setContentsMargins(0, 0, 0, 0)
        caption_layout.setSpacing(4)
        caption = QLabel(label_text)
        caption.setWordWrap(True)
        caption_layout.addWidget(caption, 1)
        caption_layout.addWidget(extra_widget, 0)
        layout.addWidget(caption_row)
    layout.addWidget(control)
    return container, stretch


def _build_inline_row(*items) -> QWidget:
    row = QWidget()
    layout = QHBoxLayout(row)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(8)
    for item in items:
        if len(item) == 3:
            label_text, control, stretch = item
            extra_widget = None
        else:
            label_text, control, stretch, extra_widget = item
        container, container_stretch = _build_inline_control(label_text, control, stretch, extra_widget)
        layout.addWidget(container, container_stretch)
    return row


def _build_form_panel(*rows) -> QWidget:
    panel = QWidget()
    layout = QFormLayout(panel)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(6)
    for row in rows:
        if len(row) == 2:
            label_text, control = row
            extra_widget = None
        else:
            label_text, control, extra_widget = row
        if extra_widget is None:
            layout.addRow(label_text, control)
            continue
        wrapper = QWidget()
        wrapper_layout = QHBoxLayout(wrapper)
        wrapper_layout.setContentsMargins(0, 0, 0, 0)
        wrapper_layout.setSpacing(6)
        wrapper_layout.addWidget(control, 1)
        wrapper_layout.addWidget(extra_widget, 0)
        layout.addRow(label_text, wrapper)
    return panel


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    editor = _ParticlesWorkspaceWidget(window)
    window.widgets["particles.editor"] = editor
    layout.addRow(editor)
    return layout


def populate(window, spec: dict[str, object]) -> None:
    editor = window.widgets["particles.editor"]
    editor.set_particles(_normalized_particles_spec(spec))


def collect(window, spec: dict[str, object]) -> None:
    editor = window.widgets["particles.editor"]
    spec["particles"] = editor.get_particles()
    spec.pop("particleSource", None)