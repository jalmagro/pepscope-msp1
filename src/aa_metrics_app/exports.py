"""In-memory, self-describing CSV and ZIP exports for application results."""

from __future__ import annotations

import io
import json
import re
import zipfile
from collections.abc import Mapping
from dataclasses import asdict, dataclass, is_dataclass
from typing import Any

import pandas as pd

SEQUENCE_COLUMNS = {
    "reference_sequence",
    "allelic_consensus_sequence",
    "all_call_consensus_sequence",
    "top_exact_peptide_sequence",
    "modal_all_call_sequence",
}
FORMULA_PREFIXES = ("=", "+", "-", "@", "\t", "\r")


@dataclass(frozen=True, slots=True)
class ExportArtifact:
    """One downloadable in-memory file."""

    filename: str
    content: bytes
    mime_type: str


def safe_slug(value: object) -> str:
    """Return a lowercase traversal-safe filename component."""

    normalized = re.sub(r"[^a-zA-Z0-9]+", "-", str(value).strip()).strip("-").lower()
    return normalized or "unknown"


def safe_csv(frame: pd.DataFrame) -> bytes:
    """Serialize a table while neutralizing formulas outside biological sequences."""

    export = frame.copy()
    for column in export.select_dtypes(include=["object", "string"]).columns:
        if column in SEQUENCE_COLUMNS:
            continue
        export[column] = export[column].map(
            lambda value: (
                "'" + value
                if isinstance(value, str) and value.startswith(FORMULA_PREFIXES)
                else value
            )
        )
    return export.to_csv(index=False).encode("utf-8")


def canonical_filter_rows(frame: pd.DataFrame, conservation_metric: str) -> pd.DataFrame:
    """Return the documented stable result ranking without changing columns."""

    if frame.empty:
        return frame.reset_index(drop=True)
    sort_columns = [
        name
        for name in (
            conservation_metric,
            "mean_sample_call_rate",
            "block_number",
            "aa_start",
            "aa_end",
            "window_id",
        )
        if name in frame
    ]
    descending = {conservation_metric, "mean_sample_call_rate"}
    return frame.sort_values(
        sort_columns,
        ascending=[name not in descending for name in sort_columns],
        kind="stable",
        na_position="last",
    ).reset_index(drop=True)


def _json_bytes(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(payload, indent=2, ensure_ascii=False, default=str).encode("utf-8")


def _identity_payload(identity: object) -> dict[str, Any]:
    if is_dataclass(identity) and not isinstance(identity, type):
        return asdict(identity)
    return dict(identity)  # type: ignore[arg-type]


def _definitions_for_export(definitions: pd.DataFrame, rows: pd.DataFrame) -> pd.DataFrame:
    name_column = "name" if "name" in definitions else "column_name"
    return definitions.loc[definitions[name_column].isin(rows.columns)].reset_index(drop=True)


def _readme(kind: str) -> str:
    return f"""MSP1 Peptide Conservation Explorer — {kind} export

Coordinates are one-based and inclusive. Sliding windows never cross block boundaries.
Complete peptide frequencies exclude missing, unresolved, mixed, and non-canonical calls.
Mixed calls contribute fractionally to allele-frequency/site-diversity metrics only.
Null metric values fail an active filter constraint. Missingness is a quality signal and is not
direct proof of an insertion or deletion. Presets are screening rules, not evidence of vaccine
efficacy. Read the JSON manifests and metric-definitions.csv before interpreting the table.
"""


def _zip(files: Mapping[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, content in files.items():
            if name.startswith(("/", "\\")) or ".." in name.split("/"):
                raise ValueError(f"Unsafe ZIP member path: {name}")
            archive.writestr(name, content)
    return buffer.getvalue()


def build_filter_bundle(
    rows: pd.DataFrame,
    definitions: pd.DataFrame,
    *,
    summary: Mapping[str, Any],
    database_identity: object,
    filename_context: Mapping[str, Any],
) -> ExportArtifact:
    """Build the canonical reproducible peptide-filter bundle."""

    files = {
        "selected-peptides.csv": safe_csv(rows),
        "filter-summary.json": _json_bytes(summary),
        "metric-definitions.csv": safe_csv(_definitions_for_export(definitions, rows)),
        "database-manifest.json": _json_bytes(_identity_payload(database_identity)),
        "README.txt": _readme("filtered peptide").encode("utf-8"),
    }
    filename = (
        f"msp1_{safe_slug(filename_context['stratum_type'])}-"
        f"{safe_slug(filename_context['stratum_value'])}_block-"
        f"{int(filename_context['block_number']):02d}_"
        f"{int(filename_context['window_size'])}aa_"
        f"{safe_slug(filename_context.get('preset_id', 'custom'))}_peptides.zip"
    )
    return ExportArtifact(filename, _zip(files), "application/zip")


def build_region_bundle(
    rows: pd.DataFrame,
    definitions: pd.DataFrame,
    *,
    query: Mapping[str, Any],
    database_identity: object,
    filename_context: Mapping[str, Any],
) -> ExportArtifact:
    """Build a region-query bundle using the same conventions as filter exports."""

    files = {
        "overlapping-peptides.csv": safe_csv(rows),
        "region-query.json": _json_bytes(query),
        "metric-definitions.csv": safe_csv(_definitions_for_export(definitions, rows)),
        "database-manifest.json": _json_bytes(_identity_payload(database_identity)),
        "README.txt": _readme("coordinate-region peptide").encode("utf-8"),
    }
    filename = (
        f"msp1_{safe_slug(filename_context['stratum_type'])}-"
        f"{safe_slug(filename_context['stratum_value'])}_block-"
        f"{int(filename_context['block_number']):02d}_"
        f"{int(filename_context['window_size'])}aa_region-"
        f"{int(query['start'])}-{int(query['end'])}.zip"
    )
    return ExportArtifact(filename, _zip(files), "application/zip")


def build_mutation_bundle(
    wide: pd.DataFrame,
    tidy: pd.DataFrame,
    *,
    context: Mapping[str, Any],
    database_identity: object,
) -> ExportArtifact:
    """Build a reproducible wide/tidy mutation comparison bundle."""

    manifest = {
        **dict(context),
        "allele_frequency_estimator": (
            "canonical alternate fractional count / callable canonical allele weight; mixed "
            "calls split equally among their possible canonical residues"
        ),
        "missing_fraction_estimator": "literal missing calls / all selected samples",
        "coordinates": "one-based, inclusive",
    }
    files = {
        "mutation-comparison-wide.csv": safe_csv(wide),
        "mutation-facts-tidy.csv": safe_csv(tidy),
        "mutation-query.json": _json_bytes(manifest),
        "database-manifest.json": _json_bytes(_identity_payload(database_identity)),
        "README.txt": (
            b"MSP1 mutation comparison export\n\n"
            b"Rows describe canonical reference-to-alternate amino-acid mutations. The wide "
            b"table places the selected group first. The tidy table retains estimator counts, "
            b"denominators, zero-frequency facts, missingness, unresolvedness, and all stored "
            b"groups. See mutation-query.json for interval and estimator semantics.\n"
        ),
    }
    filename = (
        f"msp1_{safe_slug(context['stratum_type'])}-{safe_slug(context['stratum_value'])}_"
        f"aa-{int(context['aa_start'])}-{int(context['aa_end'])}_mutations.zip"
    )
    return ExportArtifact(filename, _zip(files), "application/zip")
