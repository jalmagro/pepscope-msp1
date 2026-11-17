"""Plot-only preparation for locus-track hover details."""

from __future__ import annotations

import pandas as pd

LOCUS_HOVER_COLUMNS = (
    "block_offset",
    "reference_residue",
    "consensus_residue",
    "second_residue",
    "n_samples",
    "denominator",
)


def prepare_locus_hover_data(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, tuple[str, ...]]:
    """Return plot-only locus data and concise hover fields.

    Parameters
    ----------
    frame
        Locus metric rows containing canonical repository column names.

    Returns
    -------
    tuple of pandas.DataFrame and tuple of str
        A frame with a plot-only ``denominator`` alias, when its canonical source exists, and the
        ordered hover fields available in that frame.

    Notes
    -----
    The canonical ``allele_denominator_fractional`` field is not renamed or removed from the
    input. The alias exists only so an interactive tooltip can use a shorter label.
    """

    plot_frame = frame
    if "allele_denominator_fractional" in frame.columns:
        plot_frame = frame.assign(denominator=frame["allele_denominator_fractional"])
    hover = tuple(column for column in LOCUS_HOVER_COLUMNS if column in plot_frame.columns)
    return plot_frame, hover
