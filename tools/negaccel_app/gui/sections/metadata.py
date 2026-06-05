"""Metadata form section for the NegAccel GUI."""

from __future__ import annotations

from ..common import QFormLayout, QLineEdit, QTextEdit, nested_get


def build_form(window) -> QFormLayout:
    layout = QFormLayout()
    window.widgets["caseTag"] = QLineEdit()
    window.widgets["caseTag"].textEdited.connect(window.on_case_tag_edited)
    window.widgets["title"] = QLineEdit()
    window.widgets["author"] = QLineEdit()
    window.widgets["description"] = QTextEdit()
    window.widgets["description"].setMinimumHeight(100)
    layout.addRow("Case tag", window.widgets["caseTag"])
    layout.addRow("Title", window.widgets["title"])
    layout.addRow("Author", window.widgets["author"])
    layout.addRow("Description", window.widgets["description"])
    return layout


def populate(window, spec: dict[str, object]) -> None:
    window.widgets["caseTag"].setText(str(nested_get(spec, "metadata", "caseTag", default="negaccel_case")))
    window.widgets["title"].setText(str(nested_get(spec, "metadata", "title", default="")))
    window.widgets["author"].setText(str(nested_get(spec, "metadata", "author", default="")))
    window.widgets["description"].setPlainText(str(nested_get(spec, "metadata", "description", default="")))


def collect(window, spec: dict[str, object]) -> None:
    metadata = spec.setdefault("metadata", {})
    metadata["schemaVersion"] = str(metadata.get("schemaVersion", "1.0.0"))
    metadata["caseTag"] = window.widgets["caseTag"].text().strip() or "negaccel_case"
    window._set_optional_text(metadata, "title", window.widgets["title"].text())
    window._set_optional_text(metadata, "author", window.widgets["author"].text())
    window._set_optional_text(metadata, "description", window.widgets["description"].toPlainText())
