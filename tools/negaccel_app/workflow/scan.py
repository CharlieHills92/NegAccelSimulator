"""Scan expansion and execution helpers for NegAccel workflow tooling."""

from __future__ import annotations

import copy
import subprocess
import sys
from pathlib import Path
from typing import Any

from .common import (
    WorkflowError,
    apply_override_pairs,
    ensure_case_metadata,
    load_json,
    sanitize_case_tag,
    set_path,
    write_json,
)
from .runtime import authored_to_runtime_case


def build_case_tag(scan_tag: str, base_tag: str, index: int, label: str, template: str | None) -> str:
    if template:
        case_tag = template
        case_tag = case_tag.replace("{scan}", scan_tag)
        case_tag = case_tag.replace("{base}", base_tag)
        case_tag = case_tag.replace("{index}", str(index))
        case_tag = case_tag.replace("{label}", label)
    else:
        case_tag = f"{scan_tag}_{label}"
    return sanitize_case_tag(case_tag)


def resolve_scan_output_root(scan_spec_path: Path, scan_spec: dict[str, Any], cli_override: str | None) -> Path:
    configured = cli_override
    if configured is None:
        configured = scan_spec.get("outputs", {}).get("directory")
    if not isinstance(configured, str) or not configured:
        raise WorkflowError("Scan spec requires outputs.directory or --output-dir")
    root = Path(configured)
    if not root.is_absolute():
        root = (scan_spec_path.parent / root).resolve()
    return root


def resolve_scan_template(scan_spec_path: Path, scan_spec: dict[str, Any]) -> tuple[Path, str, dict[str, Any]]:
    template_block = scan_spec.get("template")
    if not isinstance(template_block, dict):
        raise WorkflowError("Scan spec requires a template object")

    template_value = template_block.get("path")
    if not isinstance(template_value, str) or not template_value:
        raise WorkflowError("Scan spec requires template.path")

    template_kind = template_block.get("kind", "runtime-case")
    if template_kind not in {"runtime-case", "authoring-case"}:
        raise WorkflowError("template.kind must be either runtime-case or authoring-case")

    template_path = (scan_spec_path.parent / template_value).resolve()
    template_document = load_json(template_path)
    if not isinstance(template_document, dict):
        raise WorkflowError("Scan template must be a JSON object")

    return template_path, template_kind, template_document


def build_scan_case(
    template_kind: str,
    template_document: dict[str, Any],
    common_overrides: list[Any],
    parameter_path: str,
    value: Any,
    case_tag: str,
) -> dict[str, Any]:
    working_document = copy.deepcopy(template_document)
    apply_override_pairs(working_document, common_overrides)
    set_path(working_document, parameter_path, value)

    if template_kind == "authoring-case":
        return authored_to_runtime_case(working_document, case_tag)

    return working_document


def expand_scan(scan_spec_path: Path, output_dir_override: str | None = None) -> dict[str, Any]:
    scan_spec = load_json(scan_spec_path)
    if not isinstance(scan_spec, dict):
        raise WorkflowError("Scan specification must be a JSON object")

    template_path, template_kind, template_document = resolve_scan_template(scan_spec_path, scan_spec)

    scan_metadata = scan_spec.get("metadata", {})
    scan_tag = scan_metadata.get("scanTag")
    if not isinstance(scan_tag, str) or not scan_tag:
        raise WorkflowError("Scan spec requires metadata.scanTag")
    scan_tag = sanitize_case_tag(scan_tag)

    scan_block = scan_spec.get("scan")
    if not isinstance(scan_block, dict):
        raise WorkflowError("Scan spec requires a scan object")

    parameter_path = scan_block.get("parameterPath")
    values = scan_block.get("values")
    if not isinstance(parameter_path, str) or not parameter_path:
        raise WorkflowError("Scan spec requires scan.parameterPath")
    if not isinstance(values, list) or not values:
        raise WorkflowError("Scan spec requires a non-empty scan.values array")

    labels = scan_block.get("labels")
    if labels is not None:
        if not isinstance(labels, list) or len(labels) != len(values):
            raise WorkflowError("scan.labels must be an array with the same length as scan.values")
        for label in labels:
            if not isinstance(label, str) or not label:
                raise WorkflowError("scan.labels must only contain non-empty strings")

    output_root = resolve_scan_output_root(scan_spec_path, scan_spec, output_dir_override)
    output_root.mkdir(parents=True, exist_ok=True)

    common_overrides = scan_spec.get("overrides", [])
    if not isinstance(common_overrides, list):
        raise WorkflowError("overrides must be an array")

    base_case_tag = template_document.get("metadata", {}).get("caseTag", template_path.stem)
    manifest_cases = []

    for index, value in enumerate(values):
        label = labels[index] if labels is not None else str(index)
        case_tag = build_case_tag(
            scan_tag=scan_tag,
            base_tag=str(base_case_tag),
            index=index,
            label=label,
            template=scan_block.get("caseTagTemplate"),
        )

        case_dir = output_root / case_tag
        case_path = case_dir / f"{case_tag}.json"
        case_config = build_scan_case(
            template_kind,
            template_document,
            common_overrides,
            parameter_path,
            value,
            case_tag,
        )
        ensure_case_metadata(case_config, case_tag, case_dir)
        write_json(case_path, case_config)

        manifest_cases.append(
            {
                "index": index,
                "caseTag": case_tag,
                "label": label,
                "value": value,
                "configPath": case_path.as_posix(),
            }
        )

    manifest_file = scan_spec.get("outputs", {}).get("manifestFile", "scan-manifest.json")
    if not isinstance(manifest_file, str) or not manifest_file:
        raise WorkflowError("outputs.manifestFile must be a non-empty string when provided")

    manifest_path = output_root / manifest_file
    manifest = {
        "metadata": {
            "schemaVersion": scan_metadata.get("schemaVersion", "1.0.0"),
            "scanTag": scan_tag,
            "templateKind": template_kind,
            "templatePath": template_path.as_posix(),
            "scanSpecPath": scan_spec_path.as_posix(),
        },
        "execution": scan_spec.get("execution", {}),
        "cases": manifest_cases,
    }
    write_json(manifest_path, manifest)
    return manifest


def run_scan(
    scan_spec_path: Path,
    simulator: str | None,
    output_dir_override: str | None,
    load_existing_override: bool | None,
) -> int:
    manifest = expand_scan(scan_spec_path, output_dir_override)
    execution = manifest.get("execution", {})

    simulator_value = simulator or execution.get("simulator")
    if not isinstance(simulator_value, str) or not simulator_value:
        raise WorkflowError("run-scan requires a simulator path via --simulator or execution.simulator")

    simulator_path = Path(simulator_value)
    if not simulator_path.is_absolute():
        simulator_path = (scan_spec_path.parent / simulator_path).resolve()

    load_existing = load_existing_override
    if load_existing is None:
        load_existing = bool(execution.get("loadExisting", False))
    stop_on_error = bool(execution.get("stopOnError", True))

    for case in manifest.get("cases", []):
        config_path = Path(case["configPath"])
        command = [simulator_path.as_posix(), config_path.as_posix()]
        if load_existing:
            command.append("1")

        print("Running", " ".join(command))
        result = subprocess.run(command, check=False)
        if result.returncode != 0:
            message = f"Simulator failed for {config_path} with exit code {result.returncode}"
            if stop_on_error:
                raise WorkflowError(message)
            print(message, file=sys.stderr)

    return 0
