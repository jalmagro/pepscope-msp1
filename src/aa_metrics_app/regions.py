"""Inclusive coordinate-region selection over one cached window frame."""

from __future__ import annotations

import pandas as pd

from .models import RegionMode, RegionQuery, RegionResult


def select_region(
    frame: pd.DataFrame,
    query: RegionQuery,
    *,
    block_start: int,
    block_end: int,
    window_size: int,
) -> RegionResult:
    """Select windows overlapping or fully containing an inclusive interval.

    Parameters
    ----------
    frame
        One selected block/window-size frame with ``aa_start`` and ``aa_end``.
    query
        One-based inclusive coordinate interval and match mode.
    block_start, block_end
        Inclusive bounds used to reject coordinates outside the selected block.
    window_size
        Fixed peptide length represented by every row.
    """

    required = {"aa_start", "aa_end", "window_size"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Window frame is missing columns: {sorted(missing)}")
    if query.start < block_start or query.end > block_end:
        raise ValueError(
            f"Region {query.start}–{query.end} lies outside block {block_start}–{block_end}"
        )
    if not frame.window_size.eq(window_size).all():
        raise ValueError("Window frame contains a different peptide length")
    warnings: list[str] = []
    if query.mode is RegionMode.CONTAIN and query.end - query.start + 1 > window_size:
        warnings.append(
            "The query interval is longer than the selected peptide size, so no peptide can "
            "contain it completely."
        )
        return RegionResult(query=query, rows=frame.iloc[0:0].copy(), warnings=tuple(warnings))
    if query.mode is RegionMode.OVERLAP:
        mask = frame.aa_start.le(query.end) & frame.aa_end.ge(query.start)
    else:
        mask = frame.aa_start.le(query.start) & frame.aa_end.ge(query.end)
    rows = frame.loc[mask].copy().sort_values(["aa_start", "aa_end"])
    if not rows.empty:
        rows.insert(
            min(6, len(rows.columns)),
            "overlap_length",
            (
                rows.aa_end.clip(upper=query.end)
                - rows.aa_start.clip(lower=query.start)
                + 1
            ).astype(int),
        )
        rows.insert(
            min(7, len(rows.columns)),
            "query_relation",
            "overlaps" if query.mode is RegionMode.OVERLAP else "contains",
        )
    return RegionResult(query=query, rows=rows, warnings=tuple(warnings))
