"""Pure selection and session-state rules for the Streamlit application."""

from __future__ import annotations

from collections.abc import MutableMapping
from typing import Any

from .models import Selection

SCIENTIFIC_STATE_KEY = "_scientific_selection_key"
REGION_STATE_KEYS = ("region_query", "region_result")
FILTER_STATE_KEYS = ("filter_recipe", "filter_result")


def invalidate_derived_state(
    state: MutableMapping[str, Any],
    previous: Selection | None,
    current: Selection,
) -> bool:
    """Clear submitted region/filter results after a scientific selection change.

    Renderer, line/bar mode, and point visibility are intentionally ignored because they do not
    alter :attr:`Selection.data_key`.

    Returns
    -------
    bool
        Whether region/filter state was cleared.
    """

    changed = previous is not None and previous.data_key != current.data_key
    if changed:
        for key in (*REGION_STATE_KEYS, *FILTER_STATE_KEYS):
            state.pop(key, None)
    state[SCIENTIFIC_STATE_KEY] = current.data_key
    return changed


def previous_data_key(state: MutableMapping[str, Any]) -> tuple[str, str, int, int] | None:
    """Return the last submitted scientific key when present."""

    value = state.get(SCIENTIFIC_STATE_KEY)
    if isinstance(value, tuple) and len(value) == 4:
        return value
    return None
