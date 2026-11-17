"""Token-aware sequence comparison and safe HTML presentation."""

from __future__ import annotations

import html
from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd


def tokenize_sequence(sequence: str) -> tuple[str, ...]:
    """Split a haplotype into one token per amino-acid coordinate.

    Bracketed mixed calls such as ``[C/V]`` or ``[C|V]`` remain one token. Missing ``-`` and
    unresolved ``!``/``*`` markers are ordinary one-coordinate tokens.
    """

    tokens: list[str] = []
    index = 0
    while index < len(sequence):
        character = sequence[index]
        if character == "[":
            end = sequence.find("]", index + 1)
            if end < 0:
                raise ValueError(f"Unclosed mixed-call token at character {index}")
            token = sequence[index : end + 1]
            if len(token) < 5 or not any(separator in token for separator in ("/", "|")):
                raise ValueError(f"Invalid mixed-call token: {token!r}")
            tokens.append(token)
            index = end + 1
        elif character == "]":
            raise ValueError(f"Unexpected closing bracket at character {index}")
        else:
            tokens.append(character)
            index += 1
    return tuple(tokens)


def normalized_token(token: str) -> str:
    """Normalize a token for order-insensitive display comparison."""

    if token.startswith("[") and token.endswith("]"):
        body = token[1:-1].replace("|", "/")
        return "[" + "/".join(sorted(body.split("/"))) + "]"
    return token


def available_sequence_rank_columns(
    frame: pd.DataFrame, requested: tuple[str, ...]
) -> list[str]:
    """Return stored rank-table fields, omitting wholly null optional distance values."""

    columns: list[str] = []
    for column in requested:
        if column not in frame.columns:
            continue
        if column == "distance_to_reference" and not frame[column].notna().any():
            continue
        columns.append(column)
    return columns


@dataclass(frozen=True, slots=True)
class AlignmentRow:
    """One labeled sequence aligned to a reference token vector."""

    label: str
    tokens: tuple[str, ...]
    differences: tuple[bool, ...]
    annotation: str | None = None


def align_sequences(
    reference: str,
    sequences: tuple[tuple[str, str | None], ...],
    *,
    annotations: Mapping[str, str | None] | None = None,
) -> tuple[AlignmentRow, ...]:
    """Align labeled sequences and attach optional safe display annotations."""

    annotations = {} if annotations is None else annotations
    reference_tokens = tokenize_sequence(reference)
    rows = [
        AlignmentRow(
            "3D7",
            reference_tokens,
            (False,) * len(reference_tokens),
            annotations.get("3D7"),
        )
    ]
    for label, sequence in sequences:
        if sequence is None:
            continue
        tokens = tokenize_sequence(sequence)
        if len(tokens) != len(reference_tokens):
            raise ValueError(
                f"{label} contains {len(tokens)} calls; expected {len(reference_tokens)}"
            )
        differences = tuple(
            normalized_token(observed) != normalized_token(expected)
            for observed, expected in zip(tokens, reference_tokens, strict=True)
        )
        rows.append(AlignmentRow(label, tokens, differences, annotations.get(label)))
    return tuple(rows)


def _token_class(token: str, different: bool) -> str:
    classes = ["seq-token"]
    if different:
        classes.append("seq-different")
    if token == "-":
        classes.append("seq-missing")
    elif token in {"!", "*"}:
        classes.append("seq-unresolved")
    elif token.startswith("["):
        classes.append("seq-mixed")
    return " ".join(classes)


def alignment_html(rows: tuple[AlignmentRow, ...], *, start: int, wrap: int = 50) -> str:
    """Render a safe, wrapped monospaced sequence alignment."""

    if not rows:
        return ""
    length = len(rows[0].tokens)
    show_annotations = any(row.annotation is not None for row in rows)
    section_class = "seq-section seq-with-annotations" if show_annotations else "seq-section"
    sections: list[str] = ["<div class='sequence-alignment'>"]
    for offset in range(0, length, wrap):
        end = min(offset + wrap, length)
        sections.append(f"<div class='{section_class}'>")
        header = (
            "<span></span>"
            f"<div class='seq-coordinates'>AA {start + offset}–{start + end - 1}</div>"
        )
        if show_annotations:
            header += "<div class='seq-frequency-header'>Freq</div>"
        sections.append(header)
        for row in rows:
            rendered = "".join(
                f"<span class='{_token_class(token, different)}'>{html.escape(token)}</span>"
                for token, different in zip(
                    row.tokens[offset:end], row.differences[offset:end], strict=True
                )
            )
            annotation = (
                f"<span class='seq-annotation'>{html.escape(row.annotation or '—')}</span>"
                if show_annotations
                else ""
            )
            sections.append(
                "<div class='seq-row'>"
                f"<span class='seq-label'>{html.escape(row.label)}</span>"
                f"<span class='seq-calls'>{rendered}</span>"
                f"{annotation}</div>"
            )
        sections.append("</div>")
    sections.append("</div>")
    return "".join(sections)
