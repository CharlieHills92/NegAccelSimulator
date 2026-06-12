"""Shared helpers for NegAccel workflow tooling."""

from __future__ import annotations

import ast
import copy
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence, Union


JsonToken = Union[str, int]

REPO_ROOT = Path(__file__).resolve().parents[3]

DEFAULT_OUTPUTS = {
    "rootDirectory": "",
    "summary": {
        "enabled": True,
        "directory": "Summary",
    },
    "plots": {
        "enabled": True,
        "directory": "Plots",
    },
    "data": {
        "enabled": True,
        "directory": "Data",
    },
    "vtk": {
        "enabled": True,
        "directory": "VTK",
        "exportGeometry": True,
        "exportSimulationState": True,
        "exportTracedParticles": True,
    },
    "logging": {
        "consoleLevel": "info",
        "fileLevel": "debug",
        "captureStdout": True,
        "writeDebugArtifacts": False,
        "structuredLogFile": "run.log",
    },
}


class WorkflowError(RuntimeError):
    """Raised for invalid workflow inputs."""


def normalize_symbol_name(name: str) -> str:
    normalized: list[str] = []
    previous_was_separator = True
    for char in name.strip().upper():
        if char.isalnum():
            normalized.append(char)
            previous_was_separator = False
            continue
        if not previous_was_separator:
            normalized.append("_")
        previous_was_separator = True

    return "".join(normalized).strip("_")


def normalize_expression_symbols(symbols: Mapping[str, float], context: str) -> dict[str, float]:
    normalized: dict[str, float] = {}
    for raw_name, raw_value in symbols.items():
        normalized_name = normalize_symbol_name(str(raw_name))
        if not normalized_name:
            raise WorkflowError(f"{context} has an empty expression symbol name")
        if isinstance(raw_value, bool) or not isinstance(raw_value, (int, float)):
            raise WorkflowError(f"{context}.{raw_name} must resolve to a numeric value")
        if normalized_name in normalized:
            raise WorkflowError(f"{context} has duplicate expression symbol '{normalized_name}'")
        normalized[normalized_name] = float(raw_value)
    return normalized


def evaluate_numeric_expression(
    expression: str,
    symbols: Mapping[str, float] | None = None,
    context: str = "expression",
) -> float:
    source = expression.strip()
    if not source:
        raise WorkflowError(f"{context} must not be empty")

    try:
        parsed = ast.parse(source, mode="eval")
    except SyntaxError as exc:
        raise WorkflowError(f"Invalid {context}: {expression}") from exc

    normalized_symbols = normalize_expression_symbols(symbols or {}, context)

    def _evaluate(node: ast.AST) -> float:
        if isinstance(node, ast.Expression):
            return _evaluate(node.body)

        if isinstance(node, ast.Constant):
            if isinstance(node.value, bool) or not isinstance(node.value, (int, float)):
                raise WorkflowError(f"{context} only supports numeric literals")
            return float(node.value)

        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            operand = _evaluate(node.operand)
            return operand if isinstance(node.op, ast.UAdd) else -operand

        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            left = _evaluate(node.left)
            right = _evaluate(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            return left / right

        if isinstance(node, ast.Name):
            normalized_name = normalize_symbol_name(node.id)
            if normalized_name not in normalized_symbols:
                raise WorkflowError(f"Unknown name '{node.id}' in {context}")
            return normalized_symbols[normalized_name]

        raise WorkflowError(f"{context} only supports numbers, names, parentheses, and + - * /")

    return float(_evaluate(parsed))


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise WorkflowError(f"JSON file not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise WorkflowError(f"Invalid JSON in {path}: {exc}") from exc


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def parse_value(raw: str) -> Any:
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        return raw


def parse_path(path: str) -> list[JsonToken]:
    if not path:
        raise WorkflowError("Configuration path must not be empty")

    if path.startswith("/"):
        tokens: list[JsonToken] = []
        for raw_token in path.split("/")[1:]:
            token = raw_token.replace("~1", "/").replace("~0", "~")
            if token.isdigit():
                tokens.append(int(token))
            else:
                tokens.append(token)
        return tokens

    tokens: list[JsonToken] = []
    token: list[str] = []
    index = 0
    while index < len(path):
        char = path[index]
        if char == ".":
            if token:
                tokens.append("".join(token))
                token = []
            index += 1
            continue

        if char == "[":
            if token:
                tokens.append("".join(token))
                token = []
            end = path.find("]", index)
            if end == -1:
                raise WorkflowError(f"Unclosed bracket in path: {path}")
            index_token = path[index + 1:end]
            if not index_token.isdigit():
                raise WorkflowError(f"Only numeric list indexes are supported in path: {path}")
            tokens.append(int(index_token))
            index = end + 1
            continue

        token.append(char)
        index += 1

    if token:
        tokens.append("".join(token))

    return tokens


def set_path(root: Any, path: str, value: Any) -> None:
    tokens = parse_path(path)
    current = root
    for token in tokens[:-1]:
        if isinstance(token, int):
            if not isinstance(current, list):
                raise WorkflowError(f"Path token {token} expects a list in {path}")
            if token < 0 or token >= len(current):
                raise WorkflowError(f"List index {token} out of range in {path}")
            current = current[token]
        else:
            if not isinstance(current, dict):
                raise WorkflowError(f"Path token '{token}' expects an object in {path}")
            if token not in current:
                raise WorkflowError(f"Missing object key '{token}' in {path}")
            current = current[token]

    final_token = tokens[-1]
    if isinstance(final_token, int):
        if not isinstance(current, list):
            raise WorkflowError(f"Path token {final_token} expects a list in {path}")
        if final_token < 0 or final_token >= len(current):
            raise WorkflowError(f"List index {final_token} out of range in {path}")
        current[final_token] = value
        return

    if not isinstance(current, dict):
        raise WorkflowError(f"Path token '{final_token}' expects an object in {path}")
    current[final_token] = value


def sanitize_case_tag(tag: str) -> str:
    sanitized = []
    for char in tag:
        if char.isalnum() or char in {"_", "-", "."}:
            sanitized.append(char)
        else:
            sanitized.append("_")
    return "".join(sanitized)


def merge_objects(target: dict[str, Any], overrides: dict[str, Any]) -> None:
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            merge_objects(target[key], value)
        else:
            target[key] = copy.deepcopy(value)


def require_object(document: dict[str, Any], key: str, context: str) -> dict[str, Any]:
    value = document.get(key)
    if not isinstance(value, dict):
        raise WorkflowError(f"{context} requires an object '{key}'")
    return value


def require_number(document: dict[str, Any], key: str, context: str) -> float:
    value = document.get(key)
    if not isinstance(value, (int, float)):
        raise WorkflowError(f"{context} requires numeric '{key}'")
    return float(value)


def require_string(document: dict[str, Any], key: str, context: str) -> str:
    value = document.get(key)
    if not isinstance(value, str) or not value:
        raise WorkflowError(f"{context} requires non-empty string '{key}'")
    return value


def apply_override_pairs(document: Any, overrides: Iterable[dict[str, Any]]) -> None:
    for override in overrides:
        path = override.get("path")
        if not isinstance(path, str) or not path:
            raise WorkflowError("Each override requires a non-empty string 'path'")
        set_path(document, path, override.get("value"))


def apply_cli_overrides(document: Any, overrides: Sequence[str]) -> None:
    for item in overrides:
        if "=" not in item:
            raise WorkflowError(f"Override must use PATH=VALUE syntax: {item}")
        path, raw_value = item.split("=", 1)
        set_path(document, path, parse_value(raw_value))


def ensure_case_metadata(case_config: dict[str, Any], case_tag: str, output_dir: Path) -> None:
    metadata = case_config.setdefault("metadata", {})
    metadata["caseTag"] = case_tag

    outputs = case_config.setdefault("outputs", {})
    outputs["rootDirectory"] = output_dir.as_posix()


def resolve_default_case_output(case_tag: str) -> Path:
    return Path(case_tag) / f"{case_tag}.json"
