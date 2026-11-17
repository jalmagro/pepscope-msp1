"""Coordinate context diagrams for peptide-region exploration."""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go

from .base import LIGHT_GRID, NEUTRAL_GREY, REFERENCE_ORANGE, STATIC_DPI

_INTERVAL_COLORS = {
    "block": NEUTRAL_GREY,
    "query": REFERENCE_ORANGE,
}
_INTERVAL_LABELS = {
    "block": "Block",
    "query": "Submitted query",
}


def region_interval_frame(
    rows: pd.DataFrame,
    *,
    block_start: int,
    block_end: int,
    query_start: int,
    query_end: int,
) -> pd.DataFrame:
    """Build renderer-neutral inclusive intervals for the block and submitted query.

    ``rows`` is retained for API compatibility with the region result workflow. Matching peptides
    belong in the result table and inspectors rather than in this compact orientation chart.
    """

    if block_start > block_end or query_start > query_end:
        raise ValueError("Region context coordinates must be increasing")
    if query_start < block_start or query_end > block_end:
        raise ValueError("Region query must lie within the selected block")
    del rows

    records = [
        {
            "interval_type": "block",
            "display_key": "Block",
            "aa_start": block_start,
            "aa_end": block_end,
            "lane": 0,
            "lane_label": "Block",
        },
        {
            "interval_type": "query",
            "display_key": "Submitted query",
            "aa_start": query_start,
            "aa_end": query_end,
            "lane": 1,
            "lane_label": "Query",
        },
    ]
    frame = pd.DataFrame.from_records(records)
    frame["aa_length"] = frame.aa_end - frame.aa_start + 1
    frame["x_left"] = frame.aa_start - 0.5
    frame["x_right"] = frame.aa_end + 0.5
    return frame


def _context_height() -> int:
    """Return a compact fixed height for the two context intervals."""

    return 260


def plotly_region_context(
    rows: pd.DataFrame,
    *,
    block_start: int,
    block_end: int,
    query_start: int,
    query_end: int,
) -> go.Figure:
    """Render the selected block and submitted interval interactively."""

    frame = region_interval_frame(
        rows,
        block_start=block_start,
        block_end=block_end,
        query_start=query_start,
        query_end=query_end,
    )
    figure = go.Figure()
    for interval_type in ("block", "query"):
        selected = frame.loc[frame.interval_type.eq(interval_type)]
        figure.add_bar(
            x=selected.aa_length,
            base=selected.x_left,
            y=selected.lane_label,
            orientation="h",
            width=0.72,
            marker={"color": _INTERVAL_COLORS[interval_type]},
            name=_INTERVAL_LABELS[interval_type],
            customdata=selected[
                ["display_key", "aa_start", "aa_end", "aa_length"]
            ].to_numpy(),
            hovertemplate=(
                "%{customdata[0]}<br>AA %{customdata[1]}–%{customdata[2]}"
                "<br>Length: %{customdata[3]} AA<extra></extra>"
            ),
        )
    figure.update_layout(
        title="Block context and submitted query",
        barmode="overlay",
        xaxis={
            "title": "MSP1 amino-acid coordinate",
            "range": [block_start - 0.5, block_end + 0.5],
            "gridcolor": LIGHT_GRID,
        },
        yaxis={
            "title": None,
            "categoryorder": "array",
            "categoryarray": frame.lane_label.tolist(),
            "autorange": "reversed",
            "tickfont": {"size": 9},
        },
        height=_context_height(),
        margin={"l": 115, "r": 20, "t": 82, "b": 50},
        paper_bgcolor="white",
        plot_bgcolor="white",
        legend={"orientation": "h", "x": 1.0, "xanchor": "right", "y": 1.14},
    )
    return figure


def matplotlib_region_context(
    rows: pd.DataFrame,
    *,
    block_start: int,
    block_end: int,
    query_start: int,
    query_end: int,
) -> plt.Figure:
    """Render the same block/query context with Matplotlib."""

    frame = region_interval_frame(
        rows,
        block_start=block_start,
        block_end=block_end,
        query_start=query_start,
        query_end=query_end,
    )
    figure_height = _context_height() / 100
    figure, axis = plt.subplots(
        figsize=(10.5, figure_height), constrained_layout=True, dpi=STATIC_DPI
    )
    for interval_type in ("block", "query"):
        selected = frame.loc[frame.interval_type.eq(interval_type)]
        axis.barh(
            selected.lane,
            selected.aa_length,
            left=selected.x_left,
            height=0.72,
            color=_INTERVAL_COLORS[interval_type],
            label=_INTERVAL_LABELS[interval_type],
        )
    axis.set(xlabel="MSP1 amino-acid coordinate", xlim=(block_start - 0.5, block_end + 0.5))
    axis.set_title(
        "Block context and submitted query",
        loc="left",
        fontsize=10,
        fontweight="semibold",
        pad=24,
    )
    axis.set_yticks(frame.lane, labels=frame.lane_label, fontsize=7)
    axis.invert_yaxis()
    axis.grid(axis="x", color=LIGHT_GRID)
    axis.spines[["right", "top"]].set_visible(False)
    axis.legend(
        loc="upper right",
        bbox_to_anchor=(1.0, 1.16),
        ncols=2,
        frameon=False,
        fontsize=8,
    )
    return figure
