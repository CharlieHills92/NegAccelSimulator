"""Domain-specific GUI form sections."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

from . import boundary, fields, geometry, interactions, metadata, outputs, setup, source


@dataclass(frozen=True)
class GuiSection:
    title: str
    build_form: Callable[[Any], Any]
    populate: Callable[[Any, dict[str, Any]], None]
    collect: Callable[[Any, dict[str, Any]], None]


FORM_SECTIONS = [
    GuiSection("Metadata", metadata.build_form, metadata.populate, metadata.collect),
    GuiSection("Geometry", geometry.build_form, geometry.populate, geometry.collect),
    GuiSection("Boundary", boundary.build_form, boundary.populate, boundary.collect),
    GuiSection("Source", source.build_form, source.populate, source.collect),
    GuiSection("Magnetic", fields.build_form, fields.populate, fields.collect),
    GuiSection("Interactions", interactions.build_form, interactions.populate, interactions.collect),
    GuiSection("Run", setup.build_form, setup.populate, setup.collect),
    GuiSection("Outputs", outputs.build_form, outputs.populate, outputs.collect),
]
