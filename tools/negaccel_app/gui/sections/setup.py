"""Run-setup form section for the NegAccel GUI."""

from __future__ import annotations

from negaccel_app.workflow.domains.setup import multigrid_level_support

from ..common import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QSpinBox,
    QVBoxLayout,
    QWidget,
    nested_get,
    set_combo_value,
)


def _solver_settings_group(title: str, rows: list[tuple[str, QWidget]]) -> QGroupBox:
    box = QGroupBox(title)
    layout = QFormLayout(box)
    for label, widget in rows:
        layout.addRow(label, widget)
    return box


def _refresh_multigrid_hint(window) -> None:
    hint = getattr(window, "_multigrid_level_hint", None)
    if hint is None:
        return

    try:
        geometry_workspace = window.widgets.get("geometry.solidsEditor")
        geometry_document = geometry_workspace.current_geometry_document() if geometry_workspace is not None else None
        if not isinstance(geometry_document, dict):
            hint.setText("MG level limit unavailable until a valid geometry is loaded.")
            return

        domain = geometry_document.get("domain")
        if not isinstance(domain, dict):
            hint.setText("MG level limit unavailable until geometry.domain is valid.")
            return

        support = multigrid_level_support(
            float(geometry_document.get("meshSizeMeters", 0.0)),
            float(domain.get("xSizeMeters", 0.0)),
            float(domain.get("ySizeMeters", 0.0)),
            float(domain.get("zSizeMeters", 0.0)),
        )
        counts = support["nodeCounts"]
        max_levels = int(support["maxLevels"])
        requested_levels = int(window.widgets["run.solver.multigrid.levels"].value())
        if requested_levels > max_levels:
            status = f"Requested {requested_levels} is too high."
        else:
            status = f"Requested {requested_levels} is allowed."
        hint.setText(
            f"Mesh nodes: {counts[0]} x {counts[1]} x {counts[2]}. "
            f"Max allowed MG levels: {max_levels}. {status} "
            "More levels are not always better; stop once the coarse-grid hierarchy is deep enough to help convergence without adding unnecessary overhead."
        )
    except Exception:
        hint.setText("MG level limit unavailable until geometry and mesh values are valid.")


def _sync_solver_visibility(window) -> None:
    solver_type = window.widgets["run.solverType"].currentText().strip() or "bicgstab"
    is_bicgstab = solver_type == "bicgstab"
    window._bicgstab_solver_box.setVisible(is_bicgstab)
    window._multigrid_solver_box.setVisible(not is_bicgstab)
    _refresh_multigrid_hint(window)
    window.schedule_preview_refresh()


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    window.widgets["run.iterations"] = QSpinBox()
    window.widgets["run.iterations"].setRange(1, 10000)
    solver_type = QComboBox()
    solver_type.setEditable(True)
    solver_type.addItems(["bicgstab", "multigrid"])
    window.widgets["run.solverType"] = solver_type
    window.widgets["run.solver.bicgstab.eps"] = window._double_spin(1.0e-9, 1.0, 8, 1.0e-4)
    window.widgets["run.solver.bicgstab.maxIterations"] = QSpinBox()
    window.widgets["run.solver.bicgstab.maxIterations"].setRange(1, 1000000)
    window.widgets["run.solver.bicgstab.newtonEps"] = window._double_spin(1.0e-9, 1.0, 8, 1.0e-4)
    window.widgets["run.solver.bicgstab.newtonMaxIterations"] = QSpinBox()
    window.widgets["run.solver.bicgstab.newtonMaxIterations"].setRange(1, 1000000)
    window.widgets["run.solver.bicgstab.globallyConvergentNewton"] = QCheckBox("Enable globally convergent Newton")
    window.widgets["run.solver.multigrid.levels"] = QSpinBox()
    window.widgets["run.solver.multigrid.levels"].setRange(1, 32)
    window.widgets["run.solver.multigrid.mgTolerance"] = window._double_spin(1.0e-9, 1.0, 8, 1.0e-4)
    window.widgets["run.solver.multigrid.maxCycles"] = QSpinBox()
    window.widgets["run.solver.multigrid.maxCycles"].setRange(1, 1000000)
    window.widgets["run.solver.multigrid.gamma"] = QSpinBox()
    window.widgets["run.solver.multigrid.gamma"].setRange(1, 16)
    window.widgets["run.solver.multigrid.preSmooth"] = QSpinBox()
    window.widgets["run.solver.multigrid.preSmooth"].setRange(1, 1000)
    window.widgets["run.solver.multigrid.postSmooth"] = QSpinBox()
    window.widgets["run.solver.multigrid.postSmooth"].setRange(1, 1000)
    window.widgets["run.solver.multigrid.coarseRelaxation"] = window._double_spin(0.1, 3.0, 4, 1.7)
    window.widgets["run.solver.multigrid.coarseMaxIterations"] = QSpinBox()
    window.widgets["run.solver.multigrid.coarseMaxIterations"].setRange(1, 1000000)
    window.widgets["run.solver.multigrid.localPlasmaMaxIterations"] = QSpinBox()
    window.widgets["run.solver.multigrid.localPlasmaMaxIterations"].setRange(1, 1000000)
    window.widgets["execution.serverCommand"] = QLineEdit()
    window.widgets["execution.serverCommand"].setPlaceholderText("srun --pty /bin/bash")
    window.widgets["execution.serverCommand"].setToolTip(
        "Optional shell-launch command used to open a remote or cluster shell before running NegAccelExec. "
        "The GUI appends '-lc <simulation command>' to this shell launcher."
    )
    window._multigrid_level_hint = QLabel()
    window._multigrid_level_hint.setWordWrap(True)
    window.widgets["run.alphaCoeff"] = window._double_spin(0.0, 10.0, 4, 0.05)
    window.widgets["run.currentDensityTolerance"] = window._double_spin(0.0, 1000.0, 4, 0.1)
    window.widgets["run.splitDomain"] = QCheckBox("Enable domain split")
    window._bicgstab_solver_box = _solver_settings_group(
        "BiCGSTAB settings",
        [
            ("Residual tolerance", window.widgets["run.solver.bicgstab.eps"]),
            ("Max iterations", window.widgets["run.solver.bicgstab.maxIterations"]),
            ("Newton tolerance", window.widgets["run.solver.bicgstab.newtonEps"]),
            ("Newton max iterations", window.widgets["run.solver.bicgstab.newtonMaxIterations"]),
            ("", window.widgets["run.solver.bicgstab.globallyConvergentNewton"]),
        ],
    )
    window._multigrid_solver_box = _solver_settings_group(
        "Multigrid settings",
        [
            ("Levels", window.widgets["run.solver.multigrid.levels"]),
            ("MG tolerance", window.widgets["run.solver.multigrid.mgTolerance"]),
            ("Max cycles", window.widgets["run.solver.multigrid.maxCycles"]),
            ("Gamma", window.widgets["run.solver.multigrid.gamma"]),
            ("Pre-smoothing rounds", window.widgets["run.solver.multigrid.preSmooth"]),
            ("Post-smoothing rounds", window.widgets["run.solver.multigrid.postSmooth"]),
            ("Coarse relaxation", window.widgets["run.solver.multigrid.coarseRelaxation"]),
            ("Coarse max iterations", window.widgets["run.solver.multigrid.coarseMaxIterations"]),
            ("Local plasma max iterations", window.widgets["run.solver.multigrid.localPlasmaMaxIterations"]),
            ("Level feasibility", window._multigrid_level_hint),
        ],
    )
    layout.addRow("Iterations", window.widgets["run.iterations"])
    layout.addRow("Solver type", window.widgets["run.solverType"])
    layout.addRow(window._bicgstab_solver_box)
    layout.addRow(window._multigrid_solver_box)
    layout.addRow("Alpha coefficient", window.widgets["run.alphaCoeff"])
    layout.addRow("Current density tolerance", window.widgets["run.currentDensityTolerance"])
    layout.addRow(window.widgets["run.splitDomain"])
    solver_type.currentTextChanged.connect(lambda _text: _sync_solver_visibility(window))
    window.widgets["run.solver.multigrid.levels"].valueChanged.connect(lambda _value: _refresh_multigrid_hint(window))
    _sync_solver_visibility(window)

    runtime_help = QLabel(
        "Run Simulation first generates a runtime JSON here from the current authoring form, then launches NegAccelExec with that file."
    )
    runtime_help.setWordWrap(True)
    runtime_wrapper = QWidget()
    runtime_wrapper_layout = QVBoxLayout(runtime_wrapper)
    runtime_wrapper_layout.setContentsMargins(0, 0, 0, 0)
    runtime_wrapper_layout.setSpacing(4)
    runtime_wrapper_layout.addWidget(window.runtime_path_edit)
    runtime_wrapper_layout.addWidget(runtime_help)
    layout.addRow("Generated runtime JSON", runtime_wrapper)
    layout.addRow("Command to run on server", window.widgets["execution.serverCommand"])

    materialize_button_row = QHBoxLayout()
    materialize_button_row.addWidget(window.materialize_button)
    materialize_button_row.addStretch(1)
    materialize_button_wrapper = QWidget()
    materialize_button_wrapper.setLayout(materialize_button_row)
    layout.addRow("", materialize_button_wrapper)
    return layout


def populate(window, spec: dict[str, object]) -> None:
    window.widgets["run.iterations"].setValue(int(nested_get(spec, "run", "iterations", default=4)))
    set_combo_value(window.widgets["run.solverType"], str(nested_get(spec, "run", "solver", "type", default="bicgstab")))
    window.widgets["run.solver.bicgstab.eps"].setValue(
        float(nested_get(spec, "run", "solver", "bicgstab", "eps", default=1.0e-4))
    )
    window.widgets["run.solver.bicgstab.maxIterations"].setValue(
        int(nested_get(spec, "run", "solver", "bicgstab", "maxIterations", default=10000))
    )
    window.widgets["run.solver.bicgstab.newtonEps"].setValue(
        float(nested_get(spec, "run", "solver", "bicgstab", "newtonEps", default=1.0e-4))
    )
    window.widgets["run.solver.bicgstab.newtonMaxIterations"].setValue(
        int(nested_get(spec, "run", "solver", "bicgstab", "newtonMaxIterations", default=10))
    )
    window.widgets["run.solver.bicgstab.globallyConvergentNewton"].setChecked(
        bool(nested_get(spec, "run", "solver", "bicgstab", "globallyConvergentNewton", default=True))
    )
    window.widgets["run.solver.multigrid.levels"].setValue(
        int(nested_get(spec, "run", "solver", "multigrid", "levels", default=1))
    )
    window.widgets["run.solver.multigrid.mgTolerance"].setValue(
        float(nested_get(spec, "run", "solver", "multigrid", "mgTolerance", default=1.0e-4))
    )
    window.widgets["run.solver.multigrid.maxCycles"].setValue(
        int(nested_get(spec, "run", "solver", "multigrid", "maxCycles", default=100))
    )
    window.widgets["run.solver.multigrid.gamma"].setValue(
        int(nested_get(spec, "run", "solver", "multigrid", "gamma", default=1))
    )
    window.widgets["run.solver.multigrid.preSmooth"].setValue(
        int(nested_get(spec, "run", "solver", "multigrid", "preSmooth", default=5))
    )
    window.widgets["run.solver.multigrid.postSmooth"].setValue(
        int(nested_get(spec, "run", "solver", "multigrid", "postSmooth", default=5))
    )
    window.widgets["run.solver.multigrid.coarseRelaxation"].setValue(
        float(nested_get(spec, "run", "solver", "multigrid", "coarseRelaxation", default=1.7))
    )
    window.widgets["run.solver.multigrid.coarseMaxIterations"].setValue(
        int(nested_get(spec, "run", "solver", "multigrid", "coarseMaxIterations", default=10000))
    )
    window.widgets["run.solver.multigrid.localPlasmaMaxIterations"].setValue(
        int(nested_get(spec, "run", "solver", "multigrid", "localPlasmaMaxIterations", default=1))
    )
    window.widgets["run.alphaCoeff"].setValue(float(nested_get(spec, "run", "spaceCharge", "alphaCoeff", default=0.5)))
    window.widgets["run.currentDensityTolerance"].setValue(
        float(nested_get(spec, "run", "convergence", "currentDensityTolerance", default=1.0))
    )
    window.widgets["run.splitDomain"].setChecked(
        bool(nested_get(spec, "run", "domainDecomposition", "splitDomain", default=False))
    )
    window.widgets["execution.serverCommand"].setText(
        str(nested_get(spec, "execution", "serverCommand", default=""))
    )
    _refresh_multigrid_hint(window)
    _sync_solver_visibility(window)


def collect(window, spec: dict[str, object]) -> None:
    run = spec.setdefault("run", {})
    run["iterations"] = int(window.widgets["run.iterations"].value())
    solver = run.setdefault("solver", {})
    solver["type"] = window.widgets["run.solverType"].currentText().strip() or "bicgstab"
    solver["bicgstab"] = {
        "eps": float(window.widgets["run.solver.bicgstab.eps"].value()),
        "maxIterations": int(window.widgets["run.solver.bicgstab.maxIterations"].value()),
        "newtonEps": float(window.widgets["run.solver.bicgstab.newtonEps"].value()),
        "newtonMaxIterations": int(window.widgets["run.solver.bicgstab.newtonMaxIterations"].value()),
        "globallyConvergentNewton": window.widgets["run.solver.bicgstab.globallyConvergentNewton"].isChecked(),
    }
    solver["multigrid"] = {
        "levels": int(window.widgets["run.solver.multigrid.levels"].value()),
        "mgTolerance": float(window.widgets["run.solver.multigrid.mgTolerance"].value()),
        "maxCycles": int(window.widgets["run.solver.multigrid.maxCycles"].value()),
        "gamma": int(window.widgets["run.solver.multigrid.gamma"].value()),
        "preSmooth": int(window.widgets["run.solver.multigrid.preSmooth"].value()),
        "postSmooth": int(window.widgets["run.solver.multigrid.postSmooth"].value()),
        "coarseRelaxation": float(window.widgets["run.solver.multigrid.coarseRelaxation"].value()),
        "coarseMaxIterations": int(window.widgets["run.solver.multigrid.coarseMaxIterations"].value()),
        "localPlasmaMaxIterations": int(window.widgets["run.solver.multigrid.localPlasmaMaxIterations"].value()),
    }
    space_charge = run.setdefault("spaceCharge", {})
    space_charge["alphaCoeff"] = float(window.widgets["run.alphaCoeff"].value())
    convergence = run.setdefault("convergence", {})
    convergence["currentDensityTolerance"] = float(window.widgets["run.currentDensityTolerance"].value())
    domain_decomposition = run.setdefault("domainDecomposition", {})
    domain_decomposition["splitDomain"] = window.widgets["run.splitDomain"].isChecked()

    server_command = window.widgets["execution.serverCommand"].text().strip()
    execution = spec.get("execution")
    if server_command:
        if not isinstance(execution, dict):
            execution = {}
            spec["execution"] = execution
        execution["serverCommand"] = server_command
    else:
        if isinstance(execution, dict):
            execution.pop("serverCommand", None)
            if not execution:
                spec.pop("execution", None)

    _refresh_multigrid_hint(window)
