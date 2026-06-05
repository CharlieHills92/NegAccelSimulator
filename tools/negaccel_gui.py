#!/usr/bin/env python3
"""Compatibility wrapper for the NegAccel GUI package."""

from __future__ import annotations

from negaccel_app.gui import NegAccelMainWindow, build_parser, main

__all__ = ["NegAccelMainWindow", "build_parser", "main"]


if __name__ == "__main__":
    raise SystemExit(main())