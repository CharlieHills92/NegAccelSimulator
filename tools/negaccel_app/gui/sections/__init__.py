"""Domain-specific GUI form sections."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

from . import boundary, diagnostics, fields, geometry, interactions, metadata, outputs, setup, source, surface_interactions


@dataclass(frozen=True)
class GuiSection:
    title: str
    build_form: Callable[[Any], Any]
    populate: Callable[[Any, dict[str, Any]], None]
    collect: Callable[[Any, dict[str, Any]], None]
    build_workspace: Callable[[Any], Any] | None = None


FORM_SECTIONS = [
    GuiSection("Metadata", metadata.build_form, metadata.populate, metadata.collect),
    GuiSection("Geometry", geometry.build_form, geometry.populate, geometry.collect, geometry.build_workspace),
    GuiSection("Boundary", boundary.build_form, boundary.populate, boundary.collect),
    GuiSection("Particles", source.build_form, source.populate, source.collect),
    GuiSection("Magnetic", fields.build_form, fields.populate, fields.collect, fields.build_workspace),
    GuiSection("Gas interactions", interactions.build_form, interactions.populate, interactions.collect, interactions.build_workspace),
    GuiSection("Surface interactions", surface_interactions.build_form, surface_interactions.populate, surface_interactions.collect),
    GuiSection("Run", setup.build_form, setup.populate, setup.collect),
    GuiSection("Diagnostics", diagnostics.build_form, diagnostics.populate, diagnostics.collect),
    GuiSection("Outputs", outputs.build_form, outputs.populate, outputs.collect),
]
