"""CLI entrypoint for the NegAccel GUI."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .common import QApplication
from .window import NegAccelMainWindow


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="PySide6 GUI for NegAccel authoring, execution, and visualization."
    )
    parser.add_argument("--authoring", type=Path, help="Initial authoring JSON to open.")
    parser.add_argument("--runtime", type=Path, help="Initial runtime JSON output path.")
    parser.add_argument("--simulator", type=Path, help="Simulator executable path.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    app = QApplication(sys.argv if argv is None else [sys.argv[0], *argv])
    window = NegAccelMainWindow(args.authoring, args.runtime, args.simulator)
    window.show()
    return app.exec()