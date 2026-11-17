"""Sample-normalized, renderer-independent residue-composition logos."""

from __future__ import annotations

import math
from dataclasses import dataclass

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

AA_COLORS = {
    "A": "#4E79A7",
    "C": "#F28E2B",
    "D": "#E15759",
    "E": "#D45087",
    "F": "#7B61A8",
    "G": "#59A14F",
    "H": "#B07AA1",
    "I": "#76B7B2",
    "K": "#FF4B4B",
    "L": "#499894",
    "M": "#86BCB6",
    "N": "#EDC948",
    "P": "#9C755F",
    "Q": "#B6992D",
    "R": "#E76F51",
    "S": "#8CD17D",
    "T": "#B4E07A",
    "V": "#2F6B8A",
    "W": "#6F4E7C",
    "Y": "#AF7AA1",
    "Uncalled": "#E6E8EB",
}


@dataclass(frozen=True, slots=True)
class LogoData:
    """Plot-ready residue spectrum and exact wrapping metadata."""

    values: pd.DataFrame
    aa_start: int
    aa_end: int
    wrap: int
    n_rows: int


def prepare_logo_data(frame: pd.DataFrame, *, wrap: int = 50) -> LogoData:
    """Convert sparse allele weights into complete sample-normalized stacks.

    Canonical residue heights are fractional counts divided by all selected samples. Any mass not
    assigned to canonical residues is represented as ``Uncalled``; it includes missing,
    unresolved, and other non-canonical evidence. Mixed calls have already been split in the
    database and therefore never appear as a separate glyph.
    """

    if wrap < 1:
        raise ValueError("wrap must be positive")
    required = {
        "aa_position",
        "reference_residue",
        "n_samples",
        "missing_fraction",
        "unresolved_fraction",
        "residue",
        "sample_normalized_height",
        "allele_frequency_callable",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Logo input is missing columns: {sorted(missing)}")
    if frame.empty:
        raise ValueError("Logo input cannot be empty")
    aa_start = int(frame.aa_position.min())
    aa_end = int(frame.aa_position.max())
    records: list[dict[str, object]] = []
    for position in range(aa_start, aa_end + 1):
        locus = frame.loc[frame.aa_position.eq(position)]
        if locus.empty:
            raise ValueError(f"Logo input has no locus row for AA {position}")
        metadata = locus.iloc[0]
        alleles = locus.loc[locus.residue.notna()].copy()
        alleles["sample_normalized_height"] = pd.to_numeric(
            alleles.sample_normalized_height, errors="coerce"
        ).fillna(0.0)
        canonical_sum = float(alleles.sample_normalized_height.sum())
        if canonical_sum > 1 + 1e-9:
            raise ValueError(f"Canonical logo height exceeds one at AA {position}")
        if alleles.empty:
            consensus = None
        else:
            maximum = float(alleles.sample_normalized_height.max())
            tied = alleles.loc[
                alleles.sample_normalized_height.sub(maximum).abs().le(1e-12), "residue"
            ].astype(str)
            reference = str(metadata.reference_residue)
            consensus = reference if reference in set(tied) else sorted(tied)[0]
        states = [
            (
                str(row.residue),
                float(row.sample_normalized_height),
                float(row.allele_frequency_callable),
            )
            for row in alleles.sort_values(
                ["sample_normalized_height", "residue"], kind="stable"
            ).itertuples()
        ]
        empty_mass = max(0.0, 1.0 - canonical_sum)
        if empty_mass > 1e-12:
            states.insert(0, ("Uncalled", empty_mass, math.nan))
        bottom = 0.0
        for residue, height, callable_frequency in states:
            row_index = (position - aa_start) // wrap
            records.append(
                {
                    "aa_position": position,
                    "row_index": row_index,
                    "row_start": aa_start + row_index * wrap,
                    "row_end": min(aa_end, aa_start + (row_index + 1) * wrap - 1),
                    "reference_residue": str(metadata.reference_residue),
                    "residue": residue,
                    "height": height,
                    "bottom": bottom,
                    "top": bottom + height,
                    "allele_frequency_callable": callable_frequency,
                    "missing_fraction": float(metadata.missing_fraction),
                    "unresolved_fraction": float(metadata.unresolved_fraction),
                    "is_consensus": residue == consensus,
                }
            )
            bottom += height
    values = pd.DataFrame.from_records(records)
    n_rows = (aa_end - aa_start) // wrap + 1
    return LogoData(values, aa_start, aa_end, wrap, n_rows)


def plotly_logo(logo: LogoData, *, title: str) -> go.Figure:
    """Render a wrapped interactive stacked residue logo."""

    figure = make_subplots(
        rows=logo.n_rows,
        cols=1,
        shared_yaxes=True,
        vertical_spacing=min(0.16, 0.12 / max(1, logo.n_rows - 1)),
        subplot_titles=[
            f"AA {int(group.row_start.iloc[0])}–{int(group.row_end.iloc[0])}"
            for _, group in logo.values.groupby("row_index", sort=True)
        ],
    )
    for row_index, group in logo.values.groupby("row_index", sort=True):
        row_start = int(group.row_start.iloc[0])
        row_end = int(group.row_end.iloc[0])
        text = [
            "" if item.residue == "Uncalled" or item.height < 0.055 else item.residue
            for item in group.itertuples()
        ]
        custom = group.loc[
            :,
            [
                "residue",
                "reference_residue",
                "allele_frequency_callable",
                "missing_fraction",
                "unresolved_fraction",
                "is_consensus",
            ],
        ].to_numpy()
        figure.add_trace(
            go.Bar(
                x=group.aa_position,
                y=group.height,
                base=group.bottom,
                width=0.92,
                marker={"color": [AA_COLORS[str(value)] for value in group.residue]},
                text=text,
                textposition="inside",
                textfont={"size": 11, "color": "#111111"},
                customdata=custom,
                hovertemplate=(
                    "AA %{x}<br>State: %{customdata[0]}<br>3D7: %{customdata[1]}"
                    "<br>Frequency: %{customdata[2]:.3f}<br>Missing: %{customdata[3]:.3f}"
                    "<br>Unresolved: %{customdata[4]:.3f}<br>Consensus: %{customdata[5]}"
                    "<extra></extra>"
                ),
                showlegend=False,
            ),
            row=int(row_index) + 1,
            col=1,
        )
        figure.update_xaxes(
            range=[row_start - 0.6, row_end + 0.6],
            dtick=max(1, math.ceil((row_end - row_start) / 20)),
            title_text="3D7 AA coordinate" if row_index == logo.n_rows - 1 else None,
            row=int(row_index) + 1,
            col=1,
        )
        figure.update_yaxes(
            range=[0, 1],
            title_text="Sample fraction",
            row=int(row_index) + 1,
            col=1,
        )
    figure.update_layout(
        title=title,
        barmode="overlay",
        height=170 + 150 * logo.n_rows,
        margin={"l": 65, "r": 25, "t": 75, "b": 45},
        plot_bgcolor="white",
        paper_bgcolor="white",
        hovermode="closest",
    )
    return figure


def matplotlib_logo(logo: LogoData, *, title: str):
    """Render the same wrapped stacks with Matplotlib."""

    figure, axes = plt.subplots(
        logo.n_rows,
        1,
        figsize=(14, 1.8 * logo.n_rows + 0.8),
        squeeze=False,
    )
    for row_index, group in logo.values.groupby("row_index", sort=True):
        axis = axes[int(row_index), 0]
        for item in group.itertuples():
            axis.bar(
                item.aa_position,
                item.height,
                bottom=item.bottom,
                width=0.92,
                color=AA_COLORS[item.residue],
                linewidth=0,
            )
            if item.residue != "Uncalled" and item.height >= 0.055:
                axis.text(
                    item.aa_position,
                    item.bottom + item.height / 2,
                    item.residue,
                    ha="center",
                    va="center",
                    fontsize=7,
                    fontweight="bold" if item.is_consensus else "normal",
                )
        row_start = int(group.row_start.iloc[0])
        row_end = int(group.row_end.iloc[0])
        axis.set(xlim=(row_start - 0.6, row_end + 0.6), ylim=(0, 1), ylabel="Sample fraction")
        axis.set_title(f"AA {row_start}–{row_end}", loc="left", fontsize=9)
        axis.grid(axis="y", color="#ECECEC", linewidth=0.6)
        axis.spines[["top", "right"]].set_visible(False)
    axes[-1, 0].set_xlabel("3D7 AA coordinate")
    figure.suptitle(title, x=0.01, ha="left", fontsize=12, fontweight="bold")
    figure.tight_layout()
    return figure
