"""Renderer-paired diagnostic charts for peptide filter results."""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go

from .base import NEUTRAL_GREY, PRIMARY_TOMATO, STATIC_DPI


def diagnostic_frame(
    all_rows: pd.DataFrame,
    retained_rows: pd.DataFrame,
    *,
    conservation_metric: str,
    comparison_metric: str,
) -> pd.DataFrame:
    """Create one renderer-neutral diagnostic table in coordinate order."""

    required = {
        "display_key",
        "aa_start",
        "aa_end",
        conservation_metric,
        comparison_metric,
    }
    missing = required.difference(all_rows.columns)
    if missing:
        raise ValueError(f"Diagnostic frame is missing columns: {sorted(missing)}")
    retained = set(retained_rows.display_key)
    frame = all_rows.loc[:, list(required)].copy()
    frame["retained"] = frame.display_key.isin(retained)
    return frame.sort_values(["aa_start", "aa_end", "display_key"]).reset_index(drop=True)


def plotly_coordinate_strip(frame: pd.DataFrame) -> go.Figure:
    """Show excluded and retained window spans across the selected block."""

    figure = go.Figure()
    for retained, color, label, y in (
        (False, NEUTRAL_GREY, "Excluded", 0),
        (True, PRIMARY_TOMATO, "Retained", 1),
    ):
        selected = frame.loc[frame.retained.eq(retained)]
        if selected.empty:
            figure.add_trace(
                go.Scatter(
                    x=[],
                    y=[],
                    mode="lines",
                    line={"color": color, "width": 5},
                    name=label,
                    legendgroup=label,
                    showlegend=True,
                    hoverinfo="skip",
                )
            )
            continue
        for trace_index, row in enumerate(selected.itertuples(index=False)):
            figure.add_trace(
                go.Scatter(
                    x=[row.aa_start, row.aa_end],
                    y=[y, y],
                    mode="lines",
                    line={"color": color, "width": 5},
                    name=label,
                    legendgroup=label,
                    showlegend=bool(trace_index == 0),
                    text=[row.display_key, row.display_key],
                    hovertemplate="%{text}<br>AA %{x}<extra></extra>",
                )
            )
    figure.update_layout(
        title="Window spans across the selected block",
        height=190,
        margin={"l": 40, "r": 20, "t": 45, "b": 45},
        xaxis_title="MSP1 amino-acid coordinate",
        yaxis={"tickvals": [0, 1], "ticktext": ["Excluded", "Retained"], "range": [-0.5, 1.5]},
    )
    return figure


def matplotlib_coordinate_strip(frame: pd.DataFrame) -> plt.Figure:
    """Static counterpart of :func:`plotly_coordinate_strip`."""

    figure, axis = plt.subplots(figsize=(10, 1.8), constrained_layout=True, dpi=STATIC_DPI)
    for retained, color, y in ((False, NEUTRAL_GREY, 0), (True, PRIMARY_TOMATO, 1)):
        for row in frame.loc[frame.retained.eq(retained)].itertuples():
            axis.plot([row.aa_start, row.aa_end], [y, y], color=color, linewidth=4)
    axis.set_yticks([0, 1], labels=["Excluded", "Retained"])
    axis.set_xlabel("MSP1 amino-acid coordinate")
    axis.set_title("Window spans across the selected block")
    axis.set_ylim(-0.5, 1.5)
    return figure


def plotly_diagnostic_scatter(
    frame: pd.DataFrame, *, conservation_metric: str, comparison_metric: str
) -> go.Figure:
    """Plot quality/diversity against the selected conservation measure."""

    figure = go.Figure()
    for retained, color, symbol, label in (
        (False, NEUTRAL_GREY, "circle", "Excluded"),
        (True, PRIMARY_TOMATO, "diamond", "Retained"),
    ):
        selected = frame.loc[frame.retained.eq(retained)]
        figure.add_scatter(
            x=selected[conservation_metric],
            y=selected[comparison_metric],
            mode="markers",
            name=label,
            marker={"color": color, "symbol": symbol, "size": 8},
            customdata=selected[["display_key", "aa_start", "aa_end"]],
            hovertemplate=(
                "%{customdata[0]} · AA %{customdata[1]}–%{customdata[2]}<br>"
                "%{x:.3f}, %{y:.3f}<extra></extra>"
            ),
        )
    figure.update_layout(
        title="Conservation versus peptide evidence",
        height=360,
        xaxis_title=conservation_metric,
        yaxis_title=comparison_metric,
        margin={"l": 60, "r": 20, "t": 50, "b": 55},
    )
    return figure


def matplotlib_diagnostic_scatter(
    frame: pd.DataFrame, *, conservation_metric: str, comparison_metric: str
) -> plt.Figure:
    """Static counterpart of :func:`plotly_diagnostic_scatter`."""

    figure, axis = plt.subplots(figsize=(7.5, 3.6), constrained_layout=True, dpi=STATIC_DPI)
    for retained, color, marker, label in (
        (False, NEUTRAL_GREY, "o", "Excluded"),
        (True, PRIMARY_TOMATO, "D", "Retained"),
    ):
        selected = frame.loc[frame.retained.eq(retained)]
        axis.scatter(
            selected[conservation_metric],
            selected[comparison_metric],
            c=color,
            marker=marker,
            label=label,
            s=28,
        )
    axis.set_xlabel(conservation_metric)
    axis.set_ylabel(comparison_metric)
    axis.set_title("Conservation versus peptide evidence")
    axis.legend(frameon=False)
    return figure
