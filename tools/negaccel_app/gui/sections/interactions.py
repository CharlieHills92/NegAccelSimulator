"""Gas interactions workspace section for the NegAccel GUI."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from negaccel_app.particles import build_family_particle_types, map_particle_kind_to_family

from ..common import ParameterBindingToggle, REPO_ROOT, WorkflowError, nested_get
from ..visualization import configure_matplotlib_canvas, configure_matplotlib_toolbar


@dataclass(frozen=True)
class _ReactionDefinition:
    reaction_id: str
    label: str
    projectile_templates: tuple[str, ...]
    default_source_path: str
    default_fit_degree: int
    projectile_fate: str = "consume"
    product_templates: tuple[tuple[str, int, str], ...] = ()
    minimum_energy_ev: float | None = None
    maximum_energy_ev: float | None = None
    scale_energy_by_ion_mass: bool = True


_DEFAULT_DENSITY_PROFILE_PATH = "densprofiles/MTF_dens.dens"
_GAS_SPECIES_BY_FAMILY = {
    "H": "H2",
    "D": "D2",
}
_SUPPORTED_PRODUCT_SPEED_CLASSES = {"fast", "slow"}
_REACTION_DEFINITIONS = (
    _ReactionDefinition(
        "negative_ion_single_stripping",
        "Negative ion single stripping",
        ("H-",),
        "Cross-sections/negative_ion_single_stripping.cross",
        6,
        projectile_fate="consume",
        product_templates=(("H0", 1, "fast"), ("e-", 1, "fast")),
        minimum_energy_ev=2.4,
    ),
    _ReactionDefinition(
        "negative_ion_double_stripping",
        "Negative ion double stripping",
        ("H-",),
        "Cross-sections/negative_ion_double_stripping.cross",
        6,
        projectile_fate="consume",
        product_templates=(("H+", 1, "fast"), ("e-", 2, "fast")),
        minimum_energy_ev=1.0e3,
    ),
    _ReactionDefinition(
        "background_gas_ionization",
        "Background gas ionization",
        ("H-", "H0"),
        "Cross-sections/background_gas_ionization.cross",
        4,
        projectile_fate="survive",
        product_templates=(("H2+", 1, "slow"), ("e-", 1, "slow")),
        minimum_energy_ev=10.0,
        maximum_energy_ev=1.0e7,
    ),
    _ReactionDefinition(
        "neutral_projectile_ionization",
        "Neutral projectile ionization",
        ("H0",),
        "Cross-sections/neutral_projectile_ionization.cross",
        6,
        projectile_fate="consume",
        product_templates=(("H+", 1, "fast"), ("e-", 1, "fast")),
    ),
    _ReactionDefinition(
        "positive_ion_charge_exchange",
        "Positive ion charge exchange",
        ("H+",),
        "Cross-sections/positive_ion_charge_exchange.cross",
        6,
        projectile_fate="consume",
        product_templates=(("H2+", 1, "slow"), ("H0", 1, "fast")),
    ),
)
_REACTION_DEFINITIONS_BY_ID = {
    definition.reaction_id: definition
    for definition in _REACTION_DEFINITIONS
}


def _process_id_suffix(particle_kind: str) -> str:
    return particle_kind.lower().replace("+", "plus").replace("-", "minus")


def _default_process_id(definition: _ReactionDefinition, projectile_kind: str) -> str:
    if len(definition.projectile_templates) == 1:
        return definition.reaction_id
    return f"{definition.reaction_id}_{_process_id_suffix(projectile_kind)}"


def _double_spin(minimum: float, maximum: float, decimals: int, step: float) -> QDoubleSpinBox:
    spin = QDoubleSpinBox()
    spin.setRange(minimum, maximum)
    spin.setDecimals(decimals)
    spin.setSingleStep(step)
    return spin


def _optional_energy_widget(window=None, path_provider=None, label_text: str = "") -> tuple[QWidget, QCheckBox, QDoubleSpinBox, ParameterBindingToggle | None]:
    container = QWidget()
    layout = QHBoxLayout(container)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(6)

    enabled = QCheckBox("Enabled")
    spin = _double_spin(0.0, 1.0e12, 6, 0.1)
    spin.setEnabled(False)
    enabled.toggled.connect(spin.setEnabled)

    layout.addWidget(enabled)
    if window is not None and path_provider is not None and label_text:
        wrapped = window._build_parameterized_editor(path_provider, label_text, spin)
        layout.addWidget(wrapped, 1)
        return container, enabled, spin, wrapped.parameter_toggle()
    layout.addWidget(spin, 1)
    return container, enabled, spin, None


class _PathSelector(QWidget):
    def __init__(
        self,
        window,
        dialog_title: str,
        file_filter: str,
        placeholder: str,
        changed_callback=None,
    ) -> None:
        super().__init__()
        self._window = window
        self._dialog_title = dialog_title
        self._file_filter = file_filter
        self._changed_callback = changed_callback

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        self._edit = QLineEdit()
        self._edit.setPlaceholderText(placeholder)
        self._browse = QPushButton("Browse")
        layout.addWidget(self._edit, 1)
        layout.addWidget(self._browse)

        self._edit.textChanged.connect(self._notify_change)
        self._browse.clicked.connect(self._browse_file)

    def _notify_change(self) -> None:
        if self._changed_callback is not None:
            self._changed_callback()

    def _normalize_path_for_display(self, path: Path) -> str:
        roots: list[Path] = []
        if self._window.authoring_path is not None:
            roots.append(self._window.authoring_path.parent.resolve())
        roots.append(REPO_ROOT.resolve())
        for root in roots:
            try:
                return path.relative_to(root).as_posix()
            except ValueError:
                continue
        return str(path)

    def _resolve_candidate(self, raw_path: str) -> Path:
        candidate = Path(raw_path).expanduser()
        if candidate.is_absolute():
            return candidate
        roots: list[Path] = []
        if self._window.authoring_path is not None:
            roots.append(self._window.authoring_path.parent.resolve())
        roots.append(REPO_ROOT.resolve())
        for root in roots:
            resolved = (root / candidate).resolve()
            if resolved.exists():
                return resolved
        return (roots[0] / candidate).resolve()

    def _browse_file(self) -> None:
        current_text = self._edit.text().strip()
        start_dir = str(REPO_ROOT)
        if current_text:
            start_dir = str(self._resolve_candidate(current_text).parent)
        elif self._window.authoring_path is not None:
            start_dir = str(self._window.authoring_path.parent)

        selected, _ = QFileDialog.getOpenFileName(
            self,
            self._dialog_title,
            start_dir,
            self._file_filter,
        )
        if selected:
            self._edit.setText(self._normalize_path_for_display(Path(selected).resolve()))

    def text(self) -> str:
        return self._edit.text()

    def setText(self, value: str) -> None:
        self._edit.setText(value)


class CrossSectionPlotCanvas(FigureCanvasQTAgg):
    def __init__(self) -> None:
        self.figure = Figure(figsize=(6.8, 5.4), tight_layout=True)
        super().__init__(self.figure)
        configure_matplotlib_canvas(self)
        self.plot_placeholder("Select a projectile to preview its configured reaction processes.")

    def plot_placeholder(self, message: str) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.axis("off")
        axis.text(0.5, 0.5, message, ha="center", va="center", fontsize=11)
        self.draw_idle()

    def plot_reaction(
        self,
        label: str,
        energy_ev: np.ndarray,
        sigma_m2: np.ndarray,
        fit_energy_ev: np.ndarray,
        fit_sigma_m2: np.ndarray,
        fit_degree: int,
    ) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        axis.loglog(energy_ev, sigma_m2, "o", markersize=3.5, alpha=0.85, label="Sampled data")
        axis.loglog(fit_energy_ev, fit_sigma_m2, linewidth=2.0, label=f"{fit_degree}th-degree fit")
        axis.set_xlabel("Energy [eV]")
        axis.set_ylabel("Cross section [m^2]")
        axis.set_title(label)
        axis.grid(True, which="both", alpha=0.28)
        axis.legend()
        self.draw_idle()

    def plot_projectile_reactions(
        self,
        projectile_label: str,
        gas_species: str,
        reactions: list[dict[str, Any]],
    ) -> None:
        self.figure.clear()
        axis = self.figure.add_subplot(111)
        colors = list(plt_color for plt_color in ("#c24d2c", "#355070", "#6d597a", "#2a9d8f", "#bc6c25", "#577590"))

        for index, reaction in enumerate(reactions):
            color = colors[index % len(colors)]
            label = str(reaction["label"])
            if bool(reaction.get("selected")):
                label += " [selected]"
            axis.loglog(
                np.asarray(reaction["fitEnergyEV"], dtype=float),
                np.asarray(reaction["fitSigmaM2"], dtype=float),
                color=color,
                linewidth=2.4 if bool(reaction.get("selected")) else 1.8,
                alpha=1.0 if bool(reaction.get("selected")) else 0.85,
                label=label,
            )
            axis.loglog(
                np.asarray(reaction["energyEV"], dtype=float),
                np.asarray(reaction["sigmaM2"], dtype=float),
                "o",
                color=color,
                markersize=4.0 if bool(reaction.get("selected")) else 3.0,
                alpha=0.7 if bool(reaction.get("selected")) else 0.38,
            )

        axis.set_xlabel("Energy [eV]")
        axis.set_ylabel("Cross section [m^2]")
        axis.set_title(f"{projectile_label} + {gas_species} reaction processes")
        axis.grid(True, which="both", alpha=0.28)
        axis.legend(fontsize=8)
        self.draw_idle()


class _InteractionsWorkspaceWidget(QWidget):
    def __init__(self, window) -> None:
        super().__init__()
        self._window = window
        self._curve_cache: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        self._projectile_entries: list[dict[str, str]] = []
        self._reaction_entries: list[dict[str, Any]] = []
        self._visible_reaction_entries: list[dict[str, Any]] = []
        self._fallback_species = "H2"
        self._updating_ui = False
        self._mapped_particle_family = self._active_particle_family()
        self._parameter_toggles: list[ParameterBindingToggle] = []
        self._build_ui()

    def _current_reaction_entry_index(self) -> int:
        entry = self._current_reaction_entry()
        if entry is None:
            return 0
        process_id = str(entry.get("processId") or "")
        for index, candidate in enumerate(self._reaction_entries):
            if str(candidate.get("processId") or "") == process_id:
                return index
        return 0

    def _reaction_parameter_path(self, suffix: str) -> str:
        return f"gasInteractions.reactions[{self._current_reaction_entry_index()}].{suffix}"

    def _refresh_parameter_toggles(self) -> None:
        for toggle in self._parameter_toggles:
            toggle.refresh_binding_state()

    def showEvent(self, event) -> None:
        super().showEvent(event)
        self._sync_gas_species()
        self._sync_projectiles()

    def _build_ui(self) -> None:
        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        root.addWidget(splitter)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(6, 6, 6, 6)

        gas_box = QGroupBox("Gas density")
        gas_layout = QFormLayout(gas_box)
        self._profile_name = QLineEdit("background")
        self._gas_species = QLineEdit("H2")
        self._gas_species.setReadOnly(True)
        self._density_path = _PathSelector(
            self._window,
            "Select gas density profile",
            "Density profiles (*.dens *.txt);;All files (*)",
            _DEFAULT_DENSITY_PROFILE_PATH,
            self._on_general_changed,
        )
        species_hint = QLabel("Gas species follows the active Particles family (H -> H2, D -> D2).")
        species_hint.setWordWrap(True)
        gas_layout.addRow("Profile name", self._profile_name)
        gas_layout.addRow("Gas species", self._gas_species)
        gas_layout.addRow("Profile path", self._density_path)
        gas_layout.addRow(species_hint)

        stripping_box = QGroupBox("Stripping")
        stripping_layout = QFormLayout(stripping_box)
        self._stripping_enabled = QCheckBox("Enable stripping")
        self._generate_secondaries = QCheckBox("Generate secondaries")
        self._stripping_minimum_z = _double_spin(0.0, 10.0, 6, 0.001)
        stripping_layout.addRow(self._stripping_enabled)
        stripping_layout.addRow(self._generate_secondaries)
        stripping_layout.addRow(
            "Minimum z [m]",
            self._window._build_parameterized_editor(
                lambda: "gasInteractions.stripping.minimumZMeters",
                "Minimum z [m]",
                self._stripping_minimum_z,
            ),
        )

        projectiles_box = QGroupBox("Projectiles")
        projectiles_layout = QVBoxLayout(projectiles_box)
        projectiles_hint = QLabel("Projectiles follow the Available species in the Particles tab.")
        projectiles_hint.setWordWrap(True)
        self._projectile_list = QListWidget()
        self._projectile_list.setMinimumHeight(220)
        projectiles_layout.addWidget(projectiles_hint)
        projectiles_layout.addWidget(self._projectile_list, 1)

        left_layout.addWidget(gas_box)
        left_layout.addWidget(stripping_box)
        left_layout.addWidget(projectiles_box, 1)
        left_layout.addStretch(1)

        middle = QWidget()
        middle_layout = QVBoxLayout(middle)
        middle_layout.setContentsMargins(6, 6, 6, 6)

        middle_splitter = QSplitter(Qt.Orientation.Vertical)

        reactions_box = QGroupBox("Reaction processes")
        reactions_layout = QVBoxLayout(reactions_box)
        self._reaction_scope = QLabel("Select a projectile to inspect its reactions with the background gas.")
        self._reaction_scope.setWordWrap(True)
        reaction_actions = QHBoxLayout()
        self._add_process = QPushButton("Add process")
        self._remove_process = QPushButton("Remove process")
        reaction_actions.addWidget(self._add_process)
        reaction_actions.addWidget(self._remove_process)
        reaction_actions.addStretch(1)
        self._reaction_list = QListWidget()
        self._reaction_list.setMinimumHeight(180)
        reactions_layout.addWidget(self._reaction_scope)
        reactions_layout.addLayout(reaction_actions)
        reactions_layout.addWidget(self._reaction_list, 1)

        detail_box = QGroupBox("Selected reaction fit")
        detail_layout = QVBoxLayout(detail_box)
        detail_form = QFormLayout()
        self._reaction_title = QLabel("No reaction selected")
        self._reaction_title.setStyleSheet("font-weight: 600;")
        self._reaction_name = QLineEdit()
        self._reaction_source_path = _PathSelector(
            self._window,
            "Select cross-section file",
            "Cross-section files (*.cross);;All files (*)",
            "Cross-sections/example.cross",
            self._on_reaction_changed,
        )
        self._reaction_fit_degree = QSpinBox()
        self._reaction_fit_degree.setRange(0, 6)
        self._reaction_scale_by_mass = QCheckBox("Scale incident energy by ion mass")
        self._reaction_projectile_fate = QComboBox()
        self._reaction_projectile_fate.addItem("Consume projectile", "consume")
        self._reaction_projectile_fate.addItem("Projectile survives", "survive")
        self._reaction_minimum_energy_row, self._reaction_minimum_enabled, self._reaction_minimum_energy, minimum_toggle = (
            _optional_energy_widget(
                self._window,
                lambda: self._reaction_parameter_path("minimumEnergyEV"),
                "Minimum energy [eV]",
            )
        )
        self._reaction_maximum_energy_row, self._reaction_maximum_enabled, self._reaction_maximum_energy, maximum_toggle = (
            _optional_energy_widget(
                self._window,
                lambda: self._reaction_parameter_path("maximumEnergyEV"),
                "Maximum energy [eV]",
            )
        )
        if minimum_toggle is not None:
            self._parameter_toggles.append(minimum_toggle)
        if maximum_toggle is not None:
            self._parameter_toggles.append(maximum_toggle)
        fit_degree_editor = self._window._build_parameterized_editor(
            lambda: self._reaction_parameter_path("fitDegree"),
            "Fit degree",
            self._reaction_fit_degree,
        )
        self._parameter_toggles.append(fit_degree_editor.parameter_toggle())
        detail_form.addRow("Reaction type", self._reaction_title)
        detail_form.addRow("Name", self._reaction_name)
        detail_form.addRow("Source file", self._reaction_source_path)
        detail_form.addRow("Fit degree", fit_degree_editor)
        detail_form.addRow("Projectile fate", self._reaction_projectile_fate)
        detail_form.addRow("", self._reaction_scale_by_mass)
        detail_form.addRow("Minimum energy [eV]", self._reaction_minimum_energy_row)
        detail_form.addRow("Maximum energy [eV]", self._reaction_maximum_energy_row)
        detail_layout.addLayout(detail_form)

        self._reaction_error = QLabel("")
        self._reaction_error.setStyleSheet("color: #8f1d1d;")
        self._reaction_error.setWordWrap(True)
        detail_layout.addWidget(self._reaction_error)

        coefficients_label = QLabel("Fit coefficients (highest degree first)")
        products_label = QLabel("Outcome products")
        products_hint = QLabel(
            "One line per product. Use particleKind:speedClass or particleKind:count:speedClass, for example H0:fast or e-:2:slow."
        )
        products_hint.setWordWrap(True)
        self._reaction_products = QPlainTextEdit()
        self._reaction_products.setMinimumHeight(84)
        detail_layout.addWidget(products_label)
        detail_layout.addWidget(products_hint)
        detail_layout.addWidget(self._reaction_products)
        self._coefficients = QPlainTextEdit()
        self._coefficients.setReadOnly(True)
        self._coefficients.setMinimumHeight(100)
        detail_layout.addWidget(coefficients_label)
        detail_layout.addWidget(self._coefficients)

        middle_splitter.addWidget(reactions_box)
        middle_splitter.addWidget(detail_box)
        middle_splitter.setStretchFactor(0, 0)
        middle_splitter.setStretchFactor(1, 1)
        middle_splitter.setSizes([240, 760])

        middle_layout.addWidget(middle_splitter, 1)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(6, 6, 6, 6)

        projectile_plot_box = QGroupBox("Projectile reaction plot")
        projectile_plot_layout = QVBoxLayout(projectile_plot_box)
        self._projectile_plot_scope = QLabel("Select a projectile to preview all configured reaction processes.")
        self._projectile_plot_scope.setWordWrap(True)
        self._projectile_plot_error = QLabel("")
        self._projectile_plot_error.setStyleSheet("color: #8f1d1d;")
        self._projectile_plot_error.setWordWrap(True)
        self._projectile_canvas = CrossSectionPlotCanvas()
        self._projectile_toolbar = NavigationToolbar2QT(self._projectile_canvas, projectile_plot_box)
        configure_matplotlib_toolbar(self._projectile_toolbar)
        projectile_plot_layout.addWidget(self._projectile_plot_scope)
        projectile_plot_layout.addWidget(self._projectile_plot_error)
        projectile_plot_layout.addWidget(self._projectile_toolbar)
        projectile_plot_layout.addWidget(self._projectile_canvas, 1)
        right_layout.addWidget(projectile_plot_box, 1)

        splitter.addWidget(left)
        splitter.addWidget(middle)
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 2)
        splitter.setSizes([420, 520, 760])

        self._profile_name.textChanged.connect(self._on_general_changed)
        self._stripping_enabled.toggled.connect(self._on_general_changed)
        self._generate_secondaries.toggled.connect(self._on_general_changed)
        self._stripping_minimum_z.valueChanged.connect(self._on_general_changed)
        self._projectile_list.currentRowChanged.connect(self._on_projectile_selected)
        self._reaction_list.currentRowChanged.connect(self._on_reaction_selected)
        self._add_process.clicked.connect(self._on_add_process)
        self._remove_process.clicked.connect(self._on_remove_process)
        self._reaction_name.textChanged.connect(self._on_reaction_changed)
        self._reaction_fit_degree.valueChanged.connect(self._on_reaction_changed)
        self._reaction_projectile_fate.currentIndexChanged.connect(self._on_reaction_changed)
        self._reaction_scale_by_mass.toggled.connect(self._on_reaction_changed)
        self._reaction_minimum_enabled.toggled.connect(self._on_reaction_changed)
        self._reaction_minimum_energy.valueChanged.connect(self._on_reaction_changed)
        self._reaction_maximum_enabled.toggled.connect(self._on_reaction_changed)
        self._reaction_maximum_energy.valueChanged.connect(self._on_reaction_changed)
        self._reaction_products.textChanged.connect(self._on_reaction_changed)

    def _on_general_changed(self) -> None:
        if self._updating_ui:
            return
        self._sync_gas_species()
        self._refresh_reaction_scope()
        self._refresh_projectile_plot()
        self._window.schedule_preview_refresh()

    def _on_projectile_selected(self, row: int) -> None:
        if self._updating_ui:
            return
        self._refresh_visible_reactions()

    def _on_reaction_selected(self, row: int) -> None:
        if self._updating_ui:
            return
        self._load_selected_reaction(row)

    def _on_reaction_changed(self) -> None:
        if self._updating_ui:
            return
        if not self._sync_current_reaction_entry_from_controls():
            self._refresh_projectile_plot()
            return
        self._update_current_reaction_list_item()
        self._refresh_current_reaction_plot()
        self._window.schedule_preview_refresh()

    def _on_add_process(self) -> None:
        projectile = self._current_projectile_entry()
        if projectile is None:
            return

        projectile_kind = str(projectile.get("kind") or "").strip()
        if not projectile_kind:
            return

        entry = {
            "processId": self._generate_unique_process_id(projectile_kind),
            "defaultLabel": "Custom process",
            "label": "Custom process",
            "projectileKind": projectile_kind,
            "sourcePath": "",
            "fitDegree": 0,
            "coefficients": [],
            "scaleEnergyByIonMass": True,
            "minimumEnergyEV": None,
            "maximumEnergyEV": None,
            "projectileFate": "consume",
            "products": [],
        }
        self._reaction_entries.append(entry)
        self._refresh_visible_reactions()
        for index, visible_entry in enumerate(self._visible_reaction_entries):
            if visible_entry.get("processId") == entry["processId"]:
                self._reaction_list.setCurrentRow(index)
                break
        self._window.schedule_preview_refresh()

    def _on_remove_process(self) -> None:
        entry = self._current_reaction_entry()
        if entry is None:
            return

        process_id = str(entry.get("processId") or "")
        self._reaction_entries = [
            reaction_entry
            for reaction_entry in self._reaction_entries
            if str(reaction_entry.get("processId") or "") != process_id
        ]
        self._refresh_visible_reactions()
        self._window.schedule_preview_refresh()

    def _generate_unique_process_id(self, projectile_kind: str) -> str:
        existing = {
            str(entry.get("processId") or "")
            for entry in self._reaction_entries
            if str(entry.get("processId") or "")
        }
        base = f"{_process_id_suffix(projectile_kind)}_process"
        candidate = base
        index = 1
        while candidate in existing:
            index += 1
            candidate = f"{base}_{index}"
        return candidate

    def _allowed_product_kinds(self) -> set[str]:
        return {
            str(particle_type.get("kind") or "")
            for particle_type in build_family_particle_types(self._active_particle_family())
            if isinstance(particle_type, dict) and str(particle_type.get("kind") or "")
        }

    def _format_products(self, products: list[dict[str, Any]]) -> str:
        lines: list[str] = []
        for product in products:
            kind = str(product.get("particleKind") or "").strip()
            count = int(product.get("count") or 1)
            speed_class = str(product.get("speedClass") or "").strip().lower()
            if not kind:
                continue
            if speed_class in _SUPPORTED_PRODUCT_SPEED_CLASSES:
                lines.append(f"{kind}:{speed_class}" if count == 1 else f"{kind}:{count}:{speed_class}")
            else:
                lines.append(kind if count == 1 else f"{kind}:{count}")
        return "\n".join(lines)

    def _parse_products_text(self, raw_text: str) -> list[dict[str, Any]]:
        allowed_kinds = self._allowed_product_kinds()
        parsed_products: list[dict[str, Any]] = []
        for line_number, raw_line in enumerate(raw_text.splitlines(), start=1):
            line = raw_line.strip()
            if not line:
                continue
            parts = [part.strip() for part in line.split(":")]
            if len(parts) == 2:
                raw_kind, raw_speed_class = parts
                raw_count = "1"
            elif len(parts) == 3:
                raw_kind, raw_count, raw_speed_class = parts
            else:
                raise WorkflowError(
                    f"Outcome line {line_number} must use particleKind:speedClass or particleKind:count:speedClass"
                )

            kind = raw_kind.strip()
            if kind not in allowed_kinds:
                raise WorkflowError(
                    f"Unknown product particle kind '{kind}' on outcome line {line_number}"
                )
            try:
                count = int(raw_count.strip())
            except ValueError as exc:
                raise WorkflowError(
                    f"Outcome line {line_number} must use an integer count"
                ) from exc
            if count < 1:
                raise WorkflowError(
                    f"Outcome line {line_number} must use a count greater than or equal to 1"
                )
            speed_class = raw_speed_class.strip().lower()
            if speed_class not in _SUPPORTED_PRODUCT_SPEED_CLASSES:
                raise WorkflowError(
                    f"Outcome line {line_number} must use speedClass fast or slow"
                )
            parsed_products.append({"particleKind": kind, "count": count, "speedClass": speed_class})
        return parsed_products

    def _reaction_entry_template(self, definition: _ReactionDefinition, projectile_kind: str) -> dict[str, Any]:
        family = self._active_particle_family()
        return {
            "processId": _default_process_id(definition, projectile_kind),
            "defaultLabel": definition.label,
            "label": definition.label,
            "projectileKind": projectile_kind,
            "sourcePath": definition.default_source_path,
            "fitDegree": definition.default_fit_degree,
            "coefficients": [],
            "scaleEnergyByIonMass": definition.scale_energy_by_ion_mass,
            "projectileFate": definition.projectile_fate,
            "products": [
                {
                    "particleKind": map_particle_kind_to_family(kind, family),
                    "count": count,
                    "speedClass": speed_class,
                }
                for kind, count, speed_class in definition.product_templates
            ],
            "minimumEnergyEV": definition.minimum_energy_ev,
            "maximumEnergyEV": definition.maximum_energy_ev,
        }

    def _default_reaction_entries(self) -> list[dict[str, Any]]:
        family = self._active_particle_family()
        entries: list[dict[str, Any]] = []
        for definition in _REACTION_DEFINITIONS:
            for projectile_template in definition.projectile_templates:
                entries.append(
                    self._reaction_entry_template(
                        definition,
                        map_particle_kind_to_family(projectile_template, family),
                    )
                )
        return entries

    def _normalize_products(self, raw_products: Any) -> list[dict[str, Any]]:
        products: list[dict[str, Any]] = []
        if not isinstance(raw_products, list):
            return products
        for item in raw_products:
            if not isinstance(item, dict):
                continue
            kind = str(item.get("particleKind") or "").strip()
            if not kind:
                continue
            count = item.get("count", 1)
            if not isinstance(count, int) or count < 1:
                count = 1
            speed_class = str(item.get("speedClass") or "").strip().lower()
            products.append({"particleKind": kind, "count": count, "speedClass": speed_class})
        return products

    def _legacy_entries_from_source(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        definition = _REACTION_DEFINITIONS_BY_ID.get(str(source.get("reactionId") or ""))
        if definition is None:
            return []

        family = self._active_particle_family()
        entries: list[dict[str, Any]] = []
        for projectile_template in definition.projectile_templates:
            entry = self._reaction_entry_template(
                definition,
                map_particle_kind_to_family(projectile_template, family),
            )
            if isinstance(source.get("name"), str) and source["name"].strip():
                entry["label"] = source["name"].strip()
            if isinstance(source.get("sourcePath"), str) and source["sourcePath"].strip():
                entry["sourcePath"] = source["sourcePath"].strip()
            if isinstance(source.get("fitDegree"), int) and 0 <= source["fitDegree"] <= 6:
                entry["fitDegree"] = int(source["fitDegree"])
            if isinstance(source.get("scaleEnergyByIonMass"), bool):
                entry["scaleEnergyByIonMass"] = bool(source["scaleEnergyByIonMass"])
            if isinstance(source.get("coefficients"), list):
                entry["coefficients"] = [float(value) for value in source["coefficients"] if isinstance(value, (int, float))]
            if "minimumEnergyEV" in source and isinstance(source["minimumEnergyEV"], (int, float)):
                entry["minimumEnergyEV"] = float(source["minimumEnergyEV"])
            if "maximumEnergyEV" in source and isinstance(source["maximumEnergyEV"], (int, float)):
                entry["maximumEnergyEV"] = float(source["maximumEnergyEV"])
            entries.append(entry)
        return entries

    def _normalize_generic_reaction_entry(self, source: dict[str, Any]) -> dict[str, Any]:
        projectile_kind = str(
            source.get("projectileKind")
            or map_particle_kind_to_family("H-", self._active_particle_family())
        ).strip()
        process_id = str(source.get("processId") or "").strip() or self._generate_unique_process_id(projectile_kind)
        label = str(source.get("name") or process_id).strip() or process_id
        fit_degree = int(source.get("fitDegree")) if isinstance(source.get("fitDegree"), int) and 0 <= int(source.get("fitDegree")) <= 6 else 0
        projectile_fate = str(source.get("projectileFate") or "consume").strip().lower()
        if projectile_fate not in {"consume", "survive"}:
            projectile_fate = "consume"

        entry: dict[str, Any] = {
            "processId": process_id,
            "defaultLabel": label,
            "label": label,
            "projectileKind": projectile_kind,
            "sourcePath": str(source.get("sourcePath") or "").strip(),
            "fitDegree": fit_degree,
            "coefficients": [
                float(value)
                for value in source.get("coefficients", [])
                if isinstance(value, (int, float))
            ] if isinstance(source.get("coefficients"), list) else [],
            "scaleEnergyByIonMass": bool(source.get("scaleEnergyByIonMass", True)),
            "projectileFate": projectile_fate,
            "products": self._normalize_products(source.get("products", [])),
            "minimumEnergyEV": float(source["minimumEnergyEV"]) if isinstance(source.get("minimumEnergyEV"), (int, float)) else None,
            "maximumEnergyEV": float(source["maximumEnergyEV"]) if isinstance(source.get("maximumEnergyEV"), (int, float)) else None,
        }
        return entry

    def _merge_reaction_entries(self, raw_reactions: Any) -> list[dict[str, Any]]:
        if not isinstance(raw_reactions, list) or not raw_reactions:
            return self._default_reaction_entries()

        merged: list[dict[str, Any]] = []
        seen_process_ids: set[str] = set()
        for item in raw_reactions:
            if not isinstance(item, dict):
                continue

            entries: list[dict[str, Any]] = []
            if isinstance(item.get("processId"), str) and item["processId"].strip():
                entries = [self._normalize_generic_reaction_entry(item)]
            elif isinstance(item.get("reactionId"), str) and item["reactionId"].strip():
                entries = self._legacy_entries_from_source(item)

            for entry in entries:
                process_id = str(entry.get("processId") or "").strip()
                if not process_id or process_id in seen_process_ids:
                    continue
                seen_process_ids.add(process_id)
                merged.append(entry)

        return merged or self._default_reaction_entries()

    def _refresh_reaction_labels(self) -> None:
        current_reaction_id = None
        current_entry = self._current_reaction_entry()
        if current_entry is not None:
            current_reaction_id = str(current_entry.get("processId", ""))

        self._updating_ui = True
        try:
            self._reaction_list.clear()
            for entry in self._visible_reaction_entries:
                self._reaction_list.addItem(entry["label"] + f" (deg {int(entry['fitDegree'])})")
            if self._visible_reaction_entries:
                current_row = 0
                if current_reaction_id:
                    for index, entry in enumerate(self._visible_reaction_entries):
                        if entry.get("processId") == current_reaction_id:
                            current_row = index
                            break
                self._reaction_list.setCurrentRow(current_row)
            else:
                self._reaction_list.setCurrentRow(-1)
        finally:
            self._updating_ui = False

        if self._visible_reaction_entries:
            self._load_selected_reaction(self._reaction_list.currentRow())
        else:
            self._load_selected_reaction(-1)

    def _update_current_reaction_list_item(self) -> None:
        row = self._reaction_list.currentRow()
        if row < 0 or row >= len(self._visible_reaction_entries):
            return
        item = self._reaction_list.item(row)
        if item is None:
            return
        entry = self._visible_reaction_entries[row]
        item.setText(entry["label"] + f" (deg {int(entry['fitDegree'])})")

    def _current_projectile_entry(self) -> dict[str, str] | None:
        row = self._projectile_list.currentRow()
        if row < 0 or row >= len(self._projectile_entries):
            return None
        return self._projectile_entries[row]

    def _active_particle_family(self) -> str:
        particles_editor = self._window.widgets.get("particles.editor")
        if particles_editor is not None and hasattr(particles_editor, "get_particle_family"):
            try:
                family = particles_editor.get_particle_family()
            except Exception:
                family = None
            if family in _GAS_SPECIES_BY_FAMILY:
                return str(family)
        return "D" if self._fallback_species == "D2" else "H"

    def _available_projectile_entries(self) -> list[dict[str, str]]:
        particle_types = None
        particles_editor = self._window.widgets.get("particles.editor")
        if particles_editor is not None and hasattr(particles_editor, "get_particle_types"):
            try:
                particle_types = particles_editor.get_particle_types()
            except Exception:
                particle_types = None
        if not isinstance(particle_types, list) or not particle_types:
            particle_types = build_family_particle_types(self._active_particle_family())

        projectiles: list[dict[str, str]] = []
        seen_kinds: set[str] = set()
        for particle_type in particle_types:
            if not isinstance(particle_type, dict):
                continue
            kind = str(particle_type.get("kind") or "").strip()
            if not kind or kind in seen_kinds:
                continue
            label = str(particle_type.get("name") or kind).strip() or kind
            projectiles.append({"kind": kind, "label": label})
            seen_kinds.add(kind)
        return projectiles

    def _reactions_for_projectile(self, projectile_kind: str) -> list[dict[str, Any]]:
        return [
            entry
            for entry in self._reaction_entries
            if str(entry.get("projectileKind") or "") == projectile_kind
        ]

    def _remap_reaction_entries_to_family(self) -> None:
        family = self._active_particle_family()
        if family == self._mapped_particle_family:
            return

        for entry in self._reaction_entries:
            projectile_kind = str(entry.get("projectileKind") or "")
            if projectile_kind:
                entry["projectileKind"] = map_particle_kind_to_family(projectile_kind, family)
            raw_products = entry.get("products", [])
            if not isinstance(raw_products, list):
                continue
            for product in raw_products:
                if not isinstance(product, dict):
                    continue
                product_kind = str(product.get("particleKind") or "")
                if product_kind:
                    product["particleKind"] = map_particle_kind_to_family(product_kind, family)

        self._mapped_particle_family = family

    def _sync_projectiles(self) -> None:
        self._remap_reaction_entries_to_family()
        preferred_kind = None
        current_projectile = self._current_projectile_entry()
        if current_projectile is not None:
            preferred_kind = current_projectile.get("kind")

        self._projectile_entries = self._available_projectile_entries()

        self._updating_ui = True
        try:
            self._projectile_list.clear()
            for entry in self._projectile_entries:
                self._projectile_list.addItem(entry["label"])
            if self._projectile_entries:
                current_row = 0
                if preferred_kind:
                    for index, entry in enumerate(self._projectile_entries):
                        if entry.get("kind") == preferred_kind:
                            current_row = index
                            break
                self._projectile_list.setCurrentRow(current_row)
            else:
                self._projectile_list.setCurrentRow(-1)
        finally:
            self._updating_ui = False

        self._refresh_visible_reactions()

    def _refresh_visible_reactions(self) -> None:
        projectile = self._current_projectile_entry()
        projectile_kind = str(projectile.get("kind") or "") if projectile is not None else ""
        self._visible_reaction_entries = self._reactions_for_projectile(projectile_kind)
        self._refresh_reaction_scope()
        self._refresh_reaction_labels()

    def _refresh_reaction_scope(self) -> None:
        projectile = self._current_projectile_entry()
        gas_species = self._gas_species.text().strip() or self._fallback_species or "background gas"
        if projectile is None:
            self._reaction_scope.setText("Select a projectile to inspect its reactions with the background gas.")
            return
        if self._visible_reaction_entries:
            self._reaction_scope.setText(
                f"Reaction processes for {projectile['label']} + {gas_species}."
            )
            return
        self._reaction_scope.setText(
            f"No configured gas reactions for {projectile['label']} + {gas_species}."
        )

    def _refresh_projectile_plot(self) -> None:
        projectile = self._current_projectile_entry()
        gas_species = self._gas_species.text().strip() or self._fallback_species or "background gas"

        if projectile is None:
            self._projectile_plot_scope.setText("Select a projectile to preview all configured reaction processes.")
            self._projectile_plot_error.setText("")
            self._projectile_canvas.plot_placeholder(
                "Select a projectile to preview its configured reaction processes."
            )
            return

        if not self._visible_reaction_entries:
            self._projectile_plot_scope.setText(
                f"No configured gas reactions for {projectile['label']} + {gas_species}."
            )
            self._projectile_plot_error.setText("")
            self._projectile_canvas.plot_placeholder(
                f"No configured reaction processes for {projectile['label']} + {gas_species}."
            )
            return

        self._projectile_plot_scope.setText(
            f"All configured reaction processes for {projectile['label']} + {gas_species}."
        )
        selected_entry = self._current_reaction_entry()
        selected_reaction_id = str(selected_entry.get("processId") or "") if selected_entry is not None else ""
        reactions_to_plot: list[dict[str, Any]] = []
        plot_errors: list[str] = []

        for entry in self._visible_reaction_entries:
            try:
                fit_energy, fit_sigma, coefficients, fit_curve_energy, fit_curve_sigma = self._fit_reaction_curve(entry)
                entry["coefficients"] = coefficients
                reactions_to_plot.append(
                    {
                        "label": entry["label"],
                        "energyEV": fit_energy,
                        "sigmaM2": fit_sigma,
                        "fitEnergyEV": fit_curve_energy,
                        "fitSigmaM2": fit_curve_sigma,
                        "selected": entry.get("processId") == selected_reaction_id,
                    }
                )
            except WorkflowError as exc:
                plot_errors.append(f"{entry['label']}: {exc}")

        self._projectile_plot_error.setText("\n".join(plot_errors))
        if reactions_to_plot:
            self._projectile_canvas.plot_projectile_reactions(
                projectile["label"],
                gas_species,
                reactions_to_plot,
            )
            return

        self._projectile_canvas.plot_placeholder(
            "Plot unavailable due to invalid cross-section data or fit settings."
        )

    def _current_reaction_entry(self) -> dict[str, Any] | None:
        row = self._reaction_list.currentRow()
        if row < 0 or row >= len(self._visible_reaction_entries):
            return None
        return self._visible_reaction_entries[row]

    def _sync_current_reaction_entry_from_controls(self) -> bool:
        entry = self._current_reaction_entry()
        if entry is None:
            return True

        minimum_energy = float(self._reaction_minimum_energy.value()) if self._reaction_minimum_enabled.isChecked() else None
        maximum_energy = float(self._reaction_maximum_energy.value()) if self._reaction_maximum_enabled.isChecked() else None
        if minimum_energy is not None and maximum_energy is not None and maximum_energy < minimum_energy:
            self._reaction_error.setText("Maximum energy must be greater than or equal to minimum energy.")
            return False

        try:
            products = self._parse_products_text(self._reaction_products.toPlainText())
        except WorkflowError as exc:
            self._reaction_error.setText(str(exc))
            return False

        entry["label"] = self._reaction_name.text().strip() or str(entry.get("defaultLabel") or entry.get("label") or "Reaction")
        entry["sourcePath"] = self._reaction_source_path.text().strip()
        entry["fitDegree"] = int(self._reaction_fit_degree.value())
        entry["scaleEnergyByIonMass"] = bool(self._reaction_scale_by_mass.isChecked())
        entry["projectileFate"] = str(self._reaction_projectile_fate.currentData() or "consume")
        entry["products"] = products
        entry["minimumEnergyEV"] = minimum_energy
        entry["maximumEnergyEV"] = maximum_energy
        return True

    def _resolve_data_path(self, raw_path: str) -> Path:
        candidate = Path(raw_path).expanduser()
        if candidate.is_absolute():
            return candidate
        roots: list[Path] = []
        if self._window.authoring_path is not None:
            roots.append(self._window.authoring_path.parent.resolve())
        roots.append(REPO_ROOT.resolve())
        for root in roots:
            resolved = (root / candidate).resolve()
            if resolved.exists():
                return resolved
        return (roots[0] / candidate).resolve()

    def _load_cross_section_curve(self, source_path: str) -> tuple[np.ndarray, np.ndarray]:
        if not source_path:
            raise WorkflowError("Cross-section sourcePath is required")
        resolved_path = self._resolve_data_path(source_path)
        cache_key = str(resolved_path)
        if cache_key in self._curve_cache:
            return self._curve_cache[cache_key]
        if not resolved_path.exists():
            raise WorkflowError(f"Cross-section file not found: {resolved_path}")
        try:
            data = np.loadtxt(str(resolved_path), comments="#")
        except Exception as exc:
            raise WorkflowError(f"Failed to load cross-section file: {resolved_path}") from exc
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] < 2:
            raise WorkflowError(f"Cross-section file must contain at least two numeric columns: {resolved_path}")

        energy = np.asarray(data[:, 0], dtype=float)
        sigma = np.asarray(data[:, 1], dtype=float)
        finite_mask = np.isfinite(energy) & np.isfinite(sigma) & (energy > 0.0)
        energy = energy[finite_mask]
        sigma = sigma[finite_mask]
        if energy.size == 0:
            raise WorkflowError(f"Cross-section file contains no positive finite energy samples: {resolved_path}")

        order = np.argsort(energy)
        curve = (energy[order], sigma[order])
        self._curve_cache[cache_key] = curve
        return curve

    def _fit_reaction_curve(self, entry: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, list[float], np.ndarray, np.ndarray]:
        energy, sigma = self._load_cross_section_curve(str(entry.get("sourcePath", "")).strip())
        positive_mask = sigma > 0.0
        if int(entry["fitDegree"]) + 1 > int(np.count_nonzero(positive_mask)):
            raise WorkflowError(
                f"{entry['label']} needs at least fitDegree + 1 positive samples to compute the polynomial fit"
            )
        fit_energy = energy[positive_mask]
        fit_sigma = sigma[positive_mask]
        coefficients = np.polyfit(np.log10(fit_energy), np.log10(fit_sigma), int(entry["fitDegree"]))
        fit_curve_energy = np.geomspace(fit_energy.min(), fit_energy.max(), 240)
        fit_curve_sigma = np.power(10.0, np.polyval(coefficients, np.log10(fit_curve_energy)))
        return fit_energy, fit_sigma, [float(value) for value in coefficients.tolist()], fit_curve_energy, fit_curve_sigma

    def _format_coefficients(self, coefficients: list[float]) -> str:
        lines = []
        highest_degree = len(coefficients) - 1
        for index, value in enumerate(coefficients):
            lines.append(f"a{highest_degree - index} = {value:.10e}")
        return "\n".join(lines)

    def _refresh_current_reaction_plot(self) -> None:
        entry = self._current_reaction_entry()
        if entry is None:
            self._reaction_error.setText("")
            self._coefficients.setPlainText("")
            self._refresh_projectile_plot()
            return
        try:
            fit_energy, fit_sigma, coefficients, fit_curve_energy, fit_curve_sigma = self._fit_reaction_curve(entry)
            entry["coefficients"] = coefficients
            self._reaction_error.setText("")
            self._coefficients.setPlainText(self._format_coefficients(coefficients))
        except WorkflowError as exc:
            self._reaction_error.setText(str(exc))
            self._coefficients.setPlainText("")
        self._refresh_projectile_plot()

    def _load_selected_reaction(self, row: int) -> None:
        entry = self._visible_reaction_entries[row] if 0 <= row < len(self._visible_reaction_entries) else None
        self._updating_ui = True
        try:
            if entry is None:
                self._reaction_title.setText("No reaction selected")
                self._reaction_name.setText("")
                self._reaction_source_path.setText("")
                self._reaction_fit_degree.setValue(0)
                self._reaction_projectile_fate.setCurrentIndex(0)
                self._reaction_scale_by_mass.setChecked(True)
                self._reaction_minimum_enabled.setChecked(False)
                self._reaction_minimum_energy.setValue(0.0)
                self._reaction_maximum_enabled.setChecked(False)
                self._reaction_maximum_energy.setValue(0.0)
                self._reaction_products.setPlainText("")
            else:
                self._reaction_title.setText(str(entry.get("label") or entry.get("defaultLabel") or "Reaction"))
                self._reaction_name.setText(str(entry.get("label") or ""))
                self._reaction_source_path.setText(str(entry.get("sourcePath", "")))
                self._reaction_fit_degree.setValue(int(entry.get("fitDegree", 0)))
                projectile_fate = str(entry.get("projectileFate") or "consume")
                self._reaction_projectile_fate.setCurrentIndex(0 if projectile_fate == "consume" else 1)
                self._reaction_scale_by_mass.setChecked(bool(entry.get("scaleEnergyByIonMass", True)))
                minimum_energy = entry.get("minimumEnergyEV")
                maximum_energy = entry.get("maximumEnergyEV")
                self._reaction_minimum_enabled.setChecked(minimum_energy is not None)
                self._reaction_minimum_energy.setValue(float(minimum_energy) if minimum_energy is not None else 0.0)
                self._reaction_maximum_enabled.setChecked(maximum_energy is not None)
                self._reaction_maximum_energy.setValue(float(maximum_energy) if maximum_energy is not None else 0.0)
                self._reaction_products.setPlainText(self._format_products(entry.get("products", [])))
        finally:
            self._updating_ui = False
        self._refresh_parameter_toggles()
        self._refresh_current_reaction_plot()

    def _sync_gas_species(self) -> str:
        family = None
        particles_editor = self._window.widgets.get("particles.editor")
        if particles_editor is not None and hasattr(particles_editor, "get_particle_family"):
            try:
                family = particles_editor.get_particle_family()
            except Exception:
                family = None
        species = _GAS_SPECIES_BY_FAMILY.get(str(family), self._fallback_species or "H2")
        self._gas_species.setText(species)
        return species

    def populate(self, spec: dict[str, object]) -> None:
        self._updating_ui = True
        try:
            density_profile = nested_get(spec, "gasInteractions", "densityProfile", default={})
            if not isinstance(density_profile, dict):
                density_profile = {}
            self._profile_name.setText(str(density_profile.get("name") or "background"))
            self._density_path.setText(str(density_profile.get("path") or _DEFAULT_DENSITY_PROFILE_PATH))
            self._fallback_species = str(density_profile.get("species") or "H2")

            stripping = nested_get(spec, "gasInteractions", "stripping", default={})
            if not isinstance(stripping, dict):
                stripping = {}
            self._stripping_enabled.setChecked(bool(stripping.get("enabled", False)))
            self._generate_secondaries.setChecked(bool(stripping.get("generateSecondaries", False)))
            self._stripping_minimum_z.setValue(float(stripping.get("minimumZMeters") or 0.0))

            self._reaction_entries = self._merge_reaction_entries(
                nested_get(spec, "gasInteractions", "reactions", default=[])
            )
        finally:
            self._updating_ui = False

        self._sync_gas_species()
        self._sync_projectiles()

    def _collect_reaction_entries(self) -> list[dict[str, Any]]:
        if not self._sync_current_reaction_entry_from_controls():
            raise WorkflowError(self._reaction_error.text() or "Invalid reaction settings")

        reactions: list[dict[str, Any]] = []
        for entry in self._reaction_entries:
            source_path = str(entry.get("sourcePath", "")).strip()
            if not source_path:
                raise WorkflowError(f"{entry['label']} is missing a cross-section file path")
            _, _, coefficients, _, _ = self._fit_reaction_curve(entry)
            entry["coefficients"] = coefficients
            reaction: dict[str, Any] = {
                "processId": str(entry.get("processId") or self._generate_unique_process_id(str(entry.get("projectileKind") or "H-"))),
                "name": str(entry.get("label") or entry.get("defaultLabel") or entry.get("processId") or "process"),
                "projectileKind": str(entry.get("projectileKind") or ""),
                "projectileFate": str(entry.get("projectileFate") or "consume"),
                "products": [
                    {
                        "particleKind": str(product.get("particleKind") or ""),
                        "count": int(product.get("count") or 1),
                        "speedClass": str(product.get("speedClass") or "").strip().lower(),
                    }
                    for product in entry.get("products", [])
                    if isinstance(product, dict)
                    and str(product.get("particleKind") or "")
                    and str(product.get("speedClass") or "").strip().lower() in _SUPPORTED_PRODUCT_SPEED_CLASSES
                ],
                "sourcePath": source_path,
                "fitDegree": int(entry["fitDegree"]),
                "coefficients": coefficients,
                "scaleEnergyByIonMass": bool(entry.get("scaleEnergyByIonMass", True)),
            }
            if len(reaction["products"]) != len(entry.get("products", [])):
                raise WorkflowError(
                    f"{reaction['name']} has outcome products missing a required fast/slow speedClass"
                )
            if entry.get("minimumEnergyEV") is not None:
                reaction["minimumEnergyEV"] = float(entry["minimumEnergyEV"])
            if entry.get("maximumEnergyEV") is not None:
                reaction["maximumEnergyEV"] = float(entry["maximumEnergyEV"])
            reactions.append(reaction)
        return reactions

    def collect(self) -> dict[str, object]:
        profile_name = self._profile_name.text().strip() or "background"
        density_path = self._density_path.text().strip()
        if not density_path:
            raise WorkflowError("gasInteractions.densityProfile.path must not be empty")

        gas_interactions: dict[str, object] = {
            "densityProfile": {
                "name": profile_name,
                "species": self._sync_gas_species(),
                "path": density_path,
            },
            "stripping": {
                "enabled": self._stripping_enabled.isChecked(),
                "generateSecondaries": self._generate_secondaries.isChecked(),
                "minimumZMeters": float(self._stripping_minimum_z.value()),
            },
            "reactions": self._collect_reaction_entries(),
        }
        return gas_interactions


def build_workspace(window) -> QWidget:
    workspace = _InteractionsWorkspaceWidget(window)
    window.widgets["interactions.workspace"] = workspace
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setWidget(workspace)
    return scroll


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    layout.addRow(QLabel("Gas interactions section is rendered as a dedicated workspace."))
    return layout


def populate(window, spec: dict[str, object]) -> None:
    workspace = window.widgets.get("interactions.workspace")
    if isinstance(workspace, _InteractionsWorkspaceWidget):
        workspace.populate(spec)


def collect(window, spec: dict[str, object]) -> None:
    workspace = window.widgets.get("interactions.workspace")
    if not isinstance(workspace, _InteractionsWorkspaceWidget):
        return
    spec["gasInteractions"] = workspace.collect()

