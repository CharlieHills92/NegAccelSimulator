#!/usr/bin/env python3
"""Compatibility wrapper for NegAccel workflow tooling."""

from __future__ import annotations

from negaccel_app.workflow import *  # noqa: F401,F403
from negaccel_app.workflow import __all__ as _WORKFLOW_ALL
from negaccel_app.workflow import main

__all__ = _WORKFLOW_ALL


if __name__ == "__main__":
    raise SystemExit(main())