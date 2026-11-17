"""Matplotlib implementation of shared chart specifications."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap

from ..models import PlotStyle
from .base import (
    LIGHT_GRID,
    MISSING_CELL_GREY,
    STATIC_DPI,
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
    """Render a compact static line, bar, or binary track."""

    spec.validate_frame(frame)
    hover_columns = normalized_hover_columns(spec)
    tidy = frame.loc[:, list(dict.fromkeys((spec.x, spec.y, *hover_columns)))].copy()
    figure, axis = plt.subplots(
        figsize=(10.5, spec.height / 100), constrained_layout=True, dpi=STATIC_DPI
    )
    if spec.binary:
        values, _ = normalize_binary_values(tidy[spec.y])
        coordinates = tidy[spec.x].to_numpy(dtype=float)
        if len(coordinates) == 1:
            edges = np.asarray([coordinates[0] - 0.5, coordinates[0] + 0.5])
        else:
            midpoints = (coordinates[:-1] + coordinates[1:]) / 2
            edges = np.concatenate(
                ([coordinates[0] - (midpoints[0] - coordinates[0])], midpoints,
                 [coordinates[-1] + (coordinates[-1] - midpoints[-1])])
            )
        colormap = ListedColormap([MISSING_CELL_GREY, (1, 1, 1, 0), spec.color])
        norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], colormap.N)
        axis.pcolormesh(
            edges,
            [-0.5, 0.5],
            values.reshape(1, -1),
            cmap=colormap,
            norm=norm,
            shading="flat",
            edgecolors="white",
            linewidth=0.35,
        )
    elif spec.style is PlotStyle.BAR:
        axis.bar(tidy[spec.x], tidy[spec.y], color=spec.color, width=0.82)
    else:
        axis.plot(
            tidy[spec.x],
            tidy[spec.y],
            color=spec.color,
            linewidth=1.35,
            marker="o" if show_points else None,
            markersize=3.2,
        )
    axis.set_title(spec.title, loc="left", fontsize=10, fontweight="semibold")
    axis.set_xlabel(spec.x_label, fontsize=8)
    axis.set_ylabel("" if spec.binary else spec.y_label, fontsize=8)
    if spec.binary:
        axis.set_ylim(-0.5, 0.5)
        axis.set_yticks([])
    elif spec.y_range is not None:
        axis.set_ylim(*spec.y_range)
    axis.grid(axis="y", color=LIGHT_GRID, linewidth=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(labelsize=7)
    return RenderedTrack(figure=figure, data=tidy, backend="matplotlib")
