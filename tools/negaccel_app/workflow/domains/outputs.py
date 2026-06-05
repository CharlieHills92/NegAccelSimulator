"""Output builders for workflow materialization."""

from __future__ import annotations

import copy
from typing import Any

from ..common import DEFAULT_OUTPUTS, WorkflowError, merge_objects


def build_outputs(authoring_outputs: dict[str, Any] | None) -> dict[str, Any]:
    outputs = copy.deepcopy(DEFAULT_OUTPUTS)
    if authoring_outputs is None:
        return outputs
    if not isinstance(authoring_outputs, dict):
        raise WorkflowError("outputs must be an object when provided")
    merge_objects(outputs, authoring_outputs)
    return outputs
