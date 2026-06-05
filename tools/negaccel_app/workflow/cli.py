"""Command-line entry points for NegAccel workflow tooling."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

from .common import WorkflowError
from .runtime import materialize_authored_case, materialize_case
from .scan import expand_scan, run_scan


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Python-side authoring and orchestration for NegAccelSimulator canonical case JSON files.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    case_parser = subparsers.add_parser(
        "build-case",
        help="Materialize one resolved case JSON from a template plus optional overrides.",
    )
    case_parser.add_argument("template", help="Path to the canonical case JSON template.")
    case_parser.add_argument("--output", help="Output path for the resolved case JSON.")
    case_parser.add_argument("--case-tag", help="Override metadata.caseTag for the output case.")
    case_parser.add_argument(
        "--set",
        dest="overrides",
        action="append",
        default=[],
        help="Override using PATH=JSON_VALUE syntax. Supports dot paths, array indexes, and JSON Pointer.",
    )

    author_parser = subparsers.add_parser(
        "author-case",
        help="Materialize one resolved runtime case JSON from the organized Python-side authoring specification.",
    )
    author_parser.add_argument("authoring_spec", help="Path to the authoring specification JSON file.")
    author_parser.add_argument("--output", help="Output path for the resolved case JSON.")
    author_parser.add_argument("--case-tag", help="Override metadata.caseTag for the output case.")
    author_parser.add_argument(
        "--set",
        dest="overrides",
        action="append",
        default=[],
        help="Override authoring fields using PATH=JSON_VALUE syntax before materialization.",
    )

    scan_parser = subparsers.add_parser(
        "expand-scan",
        help="Expand a Python-side scan specification into per-case JSON files and a manifest.",
    )
    scan_parser.add_argument("scan_spec", help="Path to the scan specification JSON file.")
    scan_parser.add_argument("--output-dir", help="Override outputs.directory from the scan specification.")

    run_parser = subparsers.add_parser(
        "run-scan",
        help="Expand a scan specification and run the simulator for each generated case.",
    )
    run_parser.add_argument("scan_spec", help="Path to the scan specification JSON file.")
    run_parser.add_argument("--output-dir", help="Override outputs.directory from the scan specification.")
    run_parser.add_argument("--simulator", help="Path to the simulator executable, such as ./runtest_new_v2.")
    run_parser.add_argument(
        "--load-existing",
        action="store_true",
        default=None,
        help="Pass the legacy load-existing flag to the simulator instead of running a fresh solve.",
    )

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.command == "build-case":
            template_path = Path(args.template).resolve()
            output_path = Path(args.output).resolve() if args.output else None
            case_path = materialize_case(template_path, output_path, args.case_tag, args.overrides)
            print(case_path.as_posix())
            return 0

        if args.command == "author-case":
            authoring_path = Path(args.authoring_spec).resolve()
            output_path = Path(args.output).resolve() if args.output else None
            case_path = materialize_authored_case(authoring_path, output_path, args.case_tag, args.overrides)
            print(case_path.as_posix())
            return 0

        if args.command == "expand-scan":
            manifest = expand_scan(Path(args.scan_spec).resolve(), args.output_dir)
            print(json.dumps(manifest, indent=2))
            return 0

        if args.command == "run-scan":
            return run_scan(
                Path(args.scan_spec).resolve(),
                args.simulator,
                args.output_dir,
                args.load_existing,
            )

        raise WorkflowError(f"Unsupported command: {args.command}")
    except WorkflowError as exc:
        print(exc, file=sys.stderr)
        return 2
