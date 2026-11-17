"""Renderer-neutral chart specifications."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from ..models import PlotStyle

PRIMARY_TOMATO = "#FF4B4B"
CORAL = "#E05D4D"
REFERENCE_ORANGE = "#D48632"
NEUTRAL_GREY = "#637083"
LIGHT_GRID = "#E8EDF2"
STATIC_DPI = 180
PLOTLY_CONFIG = {
    "displaylogo": False,
    "scrollZoom": True,
    "modeBarButtonsToRemove": ["lasso2d", "select2d"],
}
MISSING_CELL_GREY = "#D9DEE5"

# Compatibility aliases for callers of the initial chart API. Both now resolve to the app's
# tomato primary rather than preserving the superseded blue/teal palette.
PRIMARY_BLUE = PRIMARY_TOMATO
SELECTED_TEAL = PRIMARY_TOMATO


def normalized_hover_columns(spec: TrackSpec) -> tuple[str, ...]:
    """Return unique hover fields that do not duplicate the plotted x/y fields."""

    return tuple(dict.fromkeys(name for name in spec.hover_columns if name not in {spec.x, spec.y}))


def normalize_binary_values(series: pd.Series) -> tuple[np.ndarray, tuple[str, ...]]:
    """Normalize Boolean-like values to ``-1`` missing, ``0`` false, and ``1`` true."""

    values: list[int] = []
    labels: list[str] = []
    for value in series.astype(object):
        if pd.isna(value):
            values.append(-1)
            labels.append("Missing")
        elif isinstance(value, (bool, np.bool_)):
            normalized = bool(value)
            values.append(int(normalized))
            labels.append(str(normalized))
        elif isinstance(value, (int, float, np.integer, np.floating)) and value in {0, 1}:
            normalized = bool(value)
            values.append(int(normalized))
            labels.append(str(normalized))
        else:
            raise ValueError(f"Binary track contains unsupported value: {value!r}")
    return np.asarray(values, dtype=np.int8), tuple(labels)


@dataclass(frozen=True, slots=True)
class TrackSpec:
    """Describe one compact quantitative track independently of renderer.

    Parameters
    ----------
    title, x_label, y_label
        Human-readable chart labels.
    x, y
        Column names in the supplied tidy frame.
    style
        Line or bar rendering for quantitative metrics.
    color
        Accessible CSS color used by both renderers.
    y_range
        Optional explicit y-axis range.
    height
        Target pixel height for the compact track.
    binary
        Render as categorical markers instead of an interpolated continuous line.
    """

    title: str
    x_label: str
    y_label: str
    x: str
    y: str
    style: PlotStyle = PlotStyle.LINE
    color: str = PRIMARY_TOMATO
    y_range: tuple[float, float] | None = None
    height: int = 165
    binary: bool = False
    hover_columns: tuple[str, ...] = ()
    description: str = ""
    unit: str = ""

    def validate_frame(self, frame: pd.DataFrame) -> None:
        """Validate required columns and coordinate order."""

        required = {self.x, self.y, *self.hover_columns}
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"Track frame is missing columns: {sorted(missing)}")
        if not frame[self.x].is_monotonic_increasing:
            raise ValueError(f"Track coordinate {self.x!r} must be ordered")
        if self.height < 100:
            raise ValueError("Track height must be at least 100 pixels")
        if self.y_range is not None and self.y_range[0] >= self.y_range[1]:
            raise ValueError("Track y_range must be increasing")


@dataclass(frozen=True, slots=True)
class RenderedTrack:
    """Renderer output plus normalized source data for parity tests."""

    figure: object
    data: pd.DataFrame
    backend: Literal["plotly", "matplotlib"]
