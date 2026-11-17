"""Consistent, null-safe application formatting helpers."""

from __future__ import annotations

import math
from typing import Any

import pandas as pd


def is_null(value: Any) -> bool:
    """Return whether a scalar should be displayed as unavailable."""

    result = pd.isna(value)
    return bool(result) if not hasattr(result, "__len__") else False


def format_count(value: Any) -> str:
    """Format an integer count or return an em dash for null."""

    return "—" if is_null(value) else f"{int(value):,}"


def format_fraction(value: Any, *, decimals: int = 1) -> str:
    """Format a zero-to-one proportion as a percentage."""

    if is_null(value):
        return "—"
    return f"{float(value):.{decimals}%}"


def format_number(value: Any, *, decimals: int = 3) -> str:
    """Format a finite numeric value without converting null to zero."""

    if is_null(value):
        return "—"
    number = float(value)
    if not math.isfinite(number):
        return "—"
    return f"{number:,.{decimals}f}"


def format_coordinate(start: int, end: int | None = None) -> str:
    """Format one-based inclusive AA coordinates."""

    return f"AA {start:,}" if end is None or end == start else f"AA {start:,}–{end:,}"


def format_bytes(value: int) -> str:
    """Format an artifact size using binary units."""

    amount = float(value)
    for unit in ("B", "KiB", "MiB", "GiB"):
        if amount < 1024 or unit == "GiB":
            return f"{amount:.1f} {unit}"
        amount /= 1024
    return f"{amount:.1f} GiB"
