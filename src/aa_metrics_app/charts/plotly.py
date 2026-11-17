"""Plotly implementation of shared chart specifications."""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go

from ..models import PlotStyle
from .base import (
    LIGHT_GRID,
    MISSING_CELL_GREY,
    RenderedTrack,
    TrackSpec,
    normalize_binary_values,
    normalized_hover_columns,
)


def render_track(
    frame: pd.DataFrame,
    spec: TrackSpec,
    *,
    show_points: bool = False,
) -> RenderedTrack:
    """Render a compact Plotly line, bar, or binary track."""

    spec.validate_frame(frame)
    hover_columns = normalized_hover_columns(spec)
    tidy = frame.loc[:, list(dict.fromkeys((spec.x, spec.y, *hover_columns)))].copy()
    custom = tidy.loc[:, list(hover_columns)].to_numpy() if hover_columns else None
    if spec.binary:
        values, labels = normalize_binary_values(tidy[spec.y])
        binary_custom = [
            [
                [labels[index], *[tidy.iloc[index][name] for name in hover_columns]]
                for index in range(len(tidy))
            ]
        ]
        trace: go.BaseTraceType = go.Heatmap(
            x=tidy[spec.x].tolist(),
            y=[spec.title],
            z=[values.tolist()],
            zmin=-1,
            zmax=1,
            colorscale=[
                [0.0, MISSING_CELL_GREY],
                [0.249, MISSING_CELL_GREY],
                [0.25, "rgba(255,255,255,0)"],
                [0.749, "rgba(255,255,255,0)"],
                [0.75, spec.color],
                [1.0, spec.color],
            ],
            showscale=False,
            xgap=1,
            ygap=1,
            customdata=binary_custom,
        )
    elif spec.style is PlotStyle.BAR:
        trace = go.Bar(
            x=tidy[spec.x],
            y=tidy[spec.y],
            marker={"color": spec.color},
            customdata=custom,
        )
    else:
        trace = go.Scatter(
            x=tidy[spec.x],
            y=tidy[spec.y],
            mode="lines+markers" if show_points else "lines",
            line={"color": spec.color, "width": 1.6},
            marker={"size": 5},
            connectgaps=False,
            customdata=custom,
        )
    if spec.binary:
        hover = f"{spec.x}: %{{x}}<br>{spec.y}: %{{customdata[0]}}"
        custom_offset = 1
    else:
        hover = f"{spec.x}: %{{x}}<br>{spec.y}: %{{y:.4g}}"
        custom_offset = 0
    for index, column in enumerate(hover_columns, start=custom_offset):
        hover += f"<br>{column}: %{{customdata[{index}]}}"
    trace.hovertemplate = hover + "<extra></extra>"
    figure = go.Figure(trace)
    figure.update_layout(
        title={"text": spec.title, "font": {"size": 15}},
        height=spec.height,
        margin={"l": 62, "r": 20, "t": 56 if spec.height >= 200 else 40, "b": 46},
        paper_bgcolor="white",
        plot_bgcolor="white",
        showlegend=False,
        hovermode="x unified" if not spec.binary else "closest",
    )
    figure.update_xaxes(title=spec.x_label, gridcolor=LIGHT_GRID, zeroline=False)
    if spec.binary:
        figure.update_yaxes(title=None, fixedrange=True, showticklabels=False, gridcolor=LIGHT_GRID)
    else:
        figure.update_yaxes(
            title=spec.y_label,
            autorange=True,
            gridcolor=LIGHT_GRID,
            zeroline=False,
        )
    return RenderedTrack(figure=figure, data=tidy, backend="plotly")
