"""CLI entrypoint for the NegAccel GUI."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def _configure_qt_platform() -> None:
    if not sys.platform.startswith("linux"):
        return
    if os.environ.get("QT_QPA_PLATFORM"):
        return
    session_type = os.environ.get("XDG_SESSION_TYPE", "").strip().lower()
    wayland_display = os.environ.get("WAYLAND_DISPLAY", "").strip()
    if session_type == "wayland" or wayland_display:
        os.environ["QT_QPA_PLATFORM"] = "xcb"


_configure_qt_platform()

from .common import QApplication
from .window import NegAccelMainWindow


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="PySide6 GUI for NegAccel authoring, execution, and visualization."
    )
    parser.add_argument("--authoring", type=Path, help="Initial authoring JSON to open.")
    parser.add_argument("--runtime", type=Path, help="Initial runtime JSON output path.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    app = QApplication(sys.argv if argv is None else [sys.argv[0], *argv])
    window = NegAccelMainWindow(args.authoring, args.runtime)
    window.show()
    return app.exec()