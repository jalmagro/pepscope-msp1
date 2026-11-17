"""Offline structure assets, reviewed coordinate maps, and Plotly backbone views."""

from __future__ import annotations

import hashlib
import shlex
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def verify_structure_asset(path: str | Path, expected_sha256: str) -> Path:
    """Return a local asset path only when its pinned SHA-256 matches."""

    asset = Path(path)
    if not asset.is_file():
        raise FileNotFoundError(f"Local structure asset does not exist: {asset}")
    observed = hashlib.sha256(asset.read_bytes()).hexdigest()
    if observed != expected_sha256.lower():
        raise ValueError("Local structure asset checksum does not match the database")
    return asset


def _tokens(line: str) -> list[str]:
    lexer = shlex.shlex(line, posix=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    return list(lexer)


def read_mmcif_loop(path: str | Path, category: str) -> pd.DataFrame:
    """Read one row-oriented mmCIF loop without a third-party structure package."""

    lines = Path(path).read_text(encoding="utf-8").splitlines()
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue
        index += 1
        tags: list[str] = []
        while index < len(lines) and lines[index].lstrip().startswith("_"):
            tags.append(lines[index].strip())
            index += 1
        if not tags or not tags[0].startswith(category + "."):
            while index < len(lines) and lines[index].strip() != "#":
                index += 1
            continue
        pending: list[str] = []
        rows: list[list[str]] = []
        while index < len(lines):
            line = lines[index].strip()
            if line == "#":
                break
            if line and not line.startswith(";"):
                pending.extend(_tokens(line))
                while len(pending) >= len(tags):
                    rows.append(pending[: len(tags)])
                    pending = pending[len(tags) :]
            index += 1
        if pending:
            raise ValueError(f"Incomplete {category} loop row in {path}")
        return pd.DataFrame(rows, columns=[tag.split(".", 1)[1] for tag in tags])
    raise ValueError(f"mmCIF category {category!r} not found")


def top_complete_mismatches(
    reference_sequence: str,
    top_complete_sequence: str | None,
    *,
    aa_start: int,
) -> pd.DataFrame:
    """Return canonical top-complete substitutions against the interval reference."""

    columns = ["aa_position", "reference_residue", "alternate_residue", "mutation"]
    if top_complete_sequence is None or pd.isna(top_complete_sequence):
        return pd.DataFrame(columns=columns)
    top = str(top_complete_sequence)
    if len(reference_sequence) != len(top):
        raise ValueError("Top complete sequence length does not match its reference interval")
    if not set(reference_sequence).issubset(CANONICAL_AA) or not set(top).issubset(CANONICAL_AA):
        raise ValueError("Structure highlights require canonical complete sequences")
    records = []
    for offset, (reference, alternate) in enumerate(
        zip(reference_sequence, top, strict=True)
    ):
        if reference != alternate:
            position = aa_start + offset
            records.append(
                {
                    "aa_position": position,
                    "reference_residue": reference,
                    "alternate_residue": alternate,
                    "mutation": f"{reference}{position}{alternate}",
                }
            )
    return pd.DataFrame.from_records(records, columns=columns)


def backbone_from_mmcif(
    path: str | Path,
    chain_id: str,
    residue_map: pd.DataFrame,
) -> pd.DataFrame:
    """Return model-1 Cα coordinates with reviewed project coordinates where mapped."""

    atoms = read_mmcif_loop(path, "_atom_site")
    required = {
        "label_asym_id",
        "label_atom_id",
        "pdbx_PDB_model_num",
        "auth_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "label_comp_id",
    }
    missing = required.difference(atoms.columns)
    if missing:
        raise ValueError(f"mmCIF atom table is missing fields: {sorted(missing)}")
    current = atoms.loc[
        atoms.label_asym_id.eq(chain_id)
        & atoms.label_atom_id.eq("CA")
        & atoms.pdbx_PDB_model_num.eq("1")
    ].copy()
    if "label_alt_id" in current:
        current = current.loc[current.label_alt_id.isin([".", "?", "A"])]
    if current.empty:
        raise ValueError(f"No model-1 C-alpha atoms found for chain {chain_id}")
    current["structure_residue_number"] = pd.to_numeric(
        current.auth_seq_id, errors="raise"
    ).astype(int)
    current["insertion_code"] = current.pdbx_PDB_ins_code.replace({".": "", "?": ""})
    current["x"] = pd.to_numeric(current.Cartn_x, errors="raise")
    current["y"] = pd.to_numeric(current.Cartn_y, errors="raise")
    current["z"] = pd.to_numeric(current.Cartn_z, errors="raise")
    current = current.sort_values(
        ["structure_residue_number", "insertion_code", "label_alt_id"]
        if "label_alt_id" in current
        else ["structure_residue_number", "insertion_code"],
        kind="stable",
    ).drop_duplicates(["structure_residue_number", "insertion_code"])
    mapping_columns = [
        "structure_residue_number",
        "insertion_code",
        "aa_position",
        "project_reference_residue",
        "mapping_status",
        "confidence_value",
    ]
    return current.loc[
        :,
        [
            "structure_residue_number",
            "insertion_code",
            "label_comp_id",
            "x",
            "y",
            "z",
        ],
    ].merge(
        residue_map.loc[:, mapping_columns],
        on=["structure_residue_number", "insertion_code"],
        how="left",
        validate="one_to_one",
    )


def map_mismatches(
    mismatches: pd.DataFrame,
    residue_map: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split top-complete mismatches into mapped and explicitly unmapped rows."""

    if mismatches.empty:
        empty = mismatches.copy()
        return empty, empty
    merged = mismatches.merge(
        residue_map.loc[
            :,
            [
                "aa_position",
                "structure_residue_number",
                "insertion_code",
                "mapping_status",
                "confidence_value",
            ],
        ],
        on="aa_position",
        how="left",
        validate="one_to_one",
    )
    mapped = merged.loc[merged.structure_residue_number.notna()].reset_index(drop=True)
    unmapped = merged.loc[merged.structure_residue_number.isna()].reset_index(drop=True)
    return mapped, unmapped


def plotly_structure_context(
    backbone: pd.DataFrame,
    *,
    aa_start: int,
    aa_end: int,
    mapped_mismatches: pd.DataFrame,
    title: str,
) -> go.Figure:
    """Render muted Cα context, selected interval, and mapped substitutions."""

    figure = go.Figure()
    figure.add_trace(
        go.Scatter3d(
            x=backbone.x,
            y=backbone.y,
            z=backbone.z,
            mode="lines",
            line={"color": "#B8BEC6", "width": 4},
            name="Structure context",
            hoverinfo="skip",
        )
    )
    interval = backbone.loc[backbone.aa_position.between(aa_start, aa_end, inclusive="both")]
    if not interval.empty:
        figure.add_trace(
            go.Scatter3d(
                x=interval.x,
                y=interval.y,
                z=interval.z,
                mode="lines+markers",
                line={"color": "#FF4B4B", "width": 8},
                marker={"color": "#FF4B4B", "size": 3},
                customdata=interval[["aa_position", "confidence_value"]].to_numpy(),
                hovertemplate=(
                    "Project AA %{customdata[0]}<br>Confidence: %{customdata[1]}"
                    "<extra>Selected interval</extra>"
                ),
                name="Selected interval",
            )
        )
    if not mapped_mismatches.empty:
        mutations = mapped_mismatches.merge(
            backbone,
            on=["structure_residue_number", "insertion_code", "aa_position"],
            how="inner",
            validate="one_to_one",
        )
        if not mutations.empty:
            figure.add_trace(
                go.Scatter3d(
                    x=mutations.x,
                    y=mutations.y,
                    z=mutations.z,
                    mode="markers",
                    marker={
                        "color": "#FFC857",
                        "size": 7,
                        "line": {"color": "#4A3420", "width": 1},
                    },
                    text=mutations.mutation,
                    hovertemplate="%{text}<extra>Top complete mismatch</extra>",
                    name="Top complete mismatches",
                )
            )
    figure.update_layout(
        title=title,
        height=620,
        margin={"l": 0, "r": 0, "t": 55, "b": 0},
        paper_bgcolor="white",
        scene={
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "zaxis": {"visible": False},
            "aspectmode": "data",
            "bgcolor": "white",
        },
        legend={"orientation": "h", "y": 0.02, "x": 0.02},
    )
    return figure
