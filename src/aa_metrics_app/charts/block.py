"""Block-level distribution charts in Plotly and Matplotlib."""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go

from .base import CORAL, LIGHT_GRID, PRIMARY_TOMATO, REFERENCE_ORANGE, STATIC_DPI


def plotly_callability_histogram(frame: pd.DataFrame, *, mean: float | None) -> go.Figure:
    """Plot stored block callability bins with exact interval hover."""

    centers = (frame.bin_left + frame.bin_right) / 2
    widths = frame.bin_right - frame.bin_left
    intervals = [
        f"{left:.2f}–{right:.2f}{']' if inclusive else ')'}"
        for left, right, inclusive in zip(
            frame.bin_left, frame.bin_right, frame.right_inclusive, strict=True
        )
    ]
    figure = go.Figure(
        go.Bar(
            x=centers,
            y=frame.sample_fraction,
            width=widths * 0.92,
            marker={"color": PRIMARY_TOMATO},
            customdata=pd.DataFrame(
                {"interval": intervals, "count": frame.sample_count, "n": frame.n_samples}
            ).to_numpy(),
            hovertemplate=(
                "Callability %{customdata[0]}<br>Samples: %{customdata[1]:,} / "
                "%{customdata[2]:,}<br>Fraction: %{y:.1%}<extra></extra>"
            ),
        )
    )
    if mean is not None and pd.notna(mean):
        figure.add_vline(
            x=float(mean),
            line_color=REFERENCE_ORANGE,
            line_dash="dash",
            annotation_text=f"Mean {float(mean):.1%}",
        )
    figure.update_layout(
        title="Sample callability distribution",
        xaxis_title="Callable fraction",
        yaxis_title="Fraction of samples",
        height=310,
        margin={"l": 60, "r": 20, "t": 50, "b": 55},
        paper_bgcolor="white",
        plot_bgcolor="white",
        showlegend=False,
    )
    figure.update_xaxes(autorange=True, gridcolor=LIGHT_GRID)
    figure.update_yaxes(autorange=True, tickformat=".0%", gridcolor=LIGHT_GRID)
    return figure


def matplotlib_callability_histogram(frame: pd.DataFrame, *, mean: float | None) -> plt.Figure:
    """Render the stored callability distribution as a static chart."""

    centers = (frame.bin_left + frame.bin_right) / 2
    widths = (frame.bin_right - frame.bin_left) * 0.92
    figure, axis = plt.subplots(figsize=(7.2, 3.1), constrained_layout=True, dpi=STATIC_DPI)
    axis.bar(centers, frame.sample_fraction, width=widths, color=PRIMARY_TOMATO)
    if mean is not None and pd.notna(mean):
        axis.axvline(float(mean), color=REFERENCE_ORANGE, linestyle="--", linewidth=1.3)
        axis.text(float(mean), axis.get_ylim()[1], f" Mean {float(mean):.1%}", va="top")
    axis.set(title="Sample callability distribution", xlabel="Callable fraction", ylabel="Samples")
    axis.set_xlim(0, 1)
    axis.yaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    axis.grid(axis="y", color=LIGHT_GRID)
    axis.spines[["top", "right"]].set_visible(False)
    return figure


def plotly_rank_distribution(frame: pd.DataFrame, *, top_n: int) -> go.Figure:
    """Plot complete haplotype frequency and cumulative coverage by rank."""

    shown = frame.head(top_n)
    colors = [
        REFERENCE_ORANGE if bool(value) else PRIMARY_TOMATO
        for value in shown.is_reference_sequence
    ]
    figure = go.Figure()
    figure.add_bar(
        x=shown.sequence_rank,
        y=shown.frequency,
        marker={"color": colors},
        name="Complete-haplotype frequency",
        customdata=shown[["observation_count", "observation_denominator", "sequence"]].to_numpy(),
        hovertemplate=(
            "Rank %{x}<br>Frequency: %{y:.2%}<br>Count: %{customdata[0]:,} / "
            "%{customdata[1]:,}<br>Sequence: %{customdata[2]}<extra></extra>"
        ),
    )
    figure.add_scatter(
        x=shown.sequence_rank,
        y=shown.cumulative_frequency,
        mode="lines+markers",
        line={"color": CORAL, "width": 1.5},
        marker={"color": CORAL, "size": 5},
        name="Cumulative",
        yaxis="y2",
    )
    figure.update_layout(
        title="Complete haplotypes by frequency rank",
        xaxis_title="Haplotype rank",
        yaxis={
            "title": "Complete-haplotype frequency",
            "tickformat": ".0%",
            "gridcolor": LIGHT_GRID,
        },
        yaxis2={
            "title": "Cumulative",
            "tickformat": ".0%",
            "autorange": True,
            "overlaying": "y",
            "side": "right",
        },
        height=310,
        margin={"l": 60, "r": 60, "t": 50, "b": 55},
        paper_bgcolor="white",
        plot_bgcolor="white",
        legend={"orientation": "h", "y": 1.13},
    )
    return figure


def matplotlib_rank_distribution(frame: pd.DataFrame, *, top_n: int) -> plt.Figure:
    """Render complete-haplotype frequency ranks and cumulative coverage statically."""

    shown = frame.head(top_n)
    colors = [
        REFERENCE_ORANGE if bool(value) else PRIMARY_TOMATO
        for value in shown.is_reference_sequence
    ]
    figure, axis = plt.subplots(figsize=(7.2, 3.1), constrained_layout=True, dpi=STATIC_DPI)
    axis.bar(shown.sequence_rank, shown.frequency, color=colors)
    cumulative = axis.twinx()
    cumulative.plot(
        shown.sequence_rank,
        shown.cumulative_frequency,
        color=CORAL,
        linewidth=1.4,
        marker="o",
        markersize=3,
    )
    axis.set(
        title="Complete haplotypes by frequency rank",
        xlabel="Rank",
        ylabel="Complete-haplotype frequency",
    )
    cumulative.set_ylabel("Cumulative")
    cumulative.set_ylim(0, 1)
    axis.yaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    cumulative.yaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    axis.grid(axis="y", color=LIGHT_GRID)
    axis.spines["top"].set_visible(False)
    cumulative.spines["top"].set_visible(False)
    return figure
