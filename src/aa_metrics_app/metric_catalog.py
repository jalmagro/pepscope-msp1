"""Application metadata layered over the authoritative database metric registry."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from .models import MetricDefinition

GROUPS: dict[str, dict[str, str]] = {
    "locus": {
        "missing_rate": "quality",
        "call_rate": "quality",
        "unresolved_rate": "quality",
        "non_callable_rate": "quality",
        "mixed_fraction_callable": "quality",
        "mixed_fraction_all_samples": "quality",
        "consensus_frequency": "allele representation",
        "reference_residue_frequency": "allele representation",
        "second_residue_frequency": "allele representation",
        "contains_reference_residue": "allele representation",
        "consensus_tied": "allele representation",
        "n_observed_residues": "residue richness",
        "n_residues_frequency_ge_0_01": "residue richness",
        "n_residues_frequency_ge_0_05": "residue richness",
        "expected_heterozygosity": "diversity and evenness",
        "shannon_entropy": "diversity and evenness",
        "allele_frequency_gini": "diversity and evenness",
        "pielou_evenness": "diversity and evenness",
    },
    "window": {
        "mean_sample_call_rate": "callability and missingness",
        "mean_sample_missing_fraction": "callability and missingness",
        "mean_sample_unresolved_fraction": "callability and missingness",
        "mean_sample_non_callable_fraction": "callability and missingness",
        "min_locus_call_rate": "callability and missingness",
        "mean_sample_mixed_fraction_callable": "callability and missingness",
        "pooled_mixed_fraction_callable": "callability and missingness",
        "max_locus_mixed_fraction_callable": "callability and missingness",
        "top_exact_peptide_frequency": "frequency and coverage",
        "reference_peptide_frequency_exact": "frequency and coverage",
        "allelic_consensus_exact_frequency": "frequency and coverage",
        "modal_all_call_sequence_frequency": "frequency and coverage",
        "peptide_diversity_hd": "diversity and distance",
        "shannon_diversity": "diversity and distance",
        "effective_peptides_simpson": "diversity and distance",
        "effective_peptides_shannon": "diversity and distance",
        "mean_pairwise_aa_distance": "diversity and distance",
        "mean_normalized_pairwise_aa_distance": "diversity and distance",
    },
    "block": {},
}


def _fallback_group(level: str, name: str) -> str:
    if name.startswith(("fraction_samples_", "n_samples", "n_total", "n_callable")):
        return "sample evidence"
    if "polymorphic" in name or "singleton_residue" in name:
        return "site composition"
    if "diversity" in name or "effective" in name or "gini" in name or "evenness" in name:
        return "diversity and evenness"
    if "frequency" in name or "cumulative" in name:
        return "frequency and coverage"
    if name.startswith("n_") or "singleton" in name:
        return "observations and richness"
    return f"{level} metrics"


def _metric_domain(data_type: str, unit: str, decimals: int) -> tuple[float, float | None, float]:
    """Return theoretical input bounds and a practical widget step."""

    minimum = 1.0 if unit == "effective_count" else 0.0
    maximum = 1.0 if unit in {"fraction", "rate", "proportion", "index"} else None
    step = 1.0 if data_type == "INTEGER" else 10.0 ** -min(max(decimals, 1), 4)
    return minimum, maximum, step


def build_metric_catalog(frame: pd.DataFrame) -> tuple[MetricDefinition, ...]:
    """Build validated application definitions from registry rows."""

    required = {
        "entity_level",
        "column_name",
        "display_name",
        "description",
        "data_type",
        "unit",
        "higher_is",
        "filterable",
        "default_visible",
        "display_decimals",
        "numerator_definition",
        "denominator_definition",
        "null_definition",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Metric registry is missing columns: {sorted(missing)}")
    keys = frame[["entity_level", "column_name"]].astype(str).agg("/".join, axis=1)
    if keys.duplicated().any():
        raise ValueError("Metric registry contains duplicate level/name rows")
    definitions: list[MetricDefinition] = []
    for row in frame.to_dict("records"):
        level = str(row["entity_level"])
        name = str(row["column_name"])
        group = GROUPS.get(level, {}).get(name, _fallback_group(level, name))
        decimals = int(row["display_decimals"])
        valid_minimum, valid_maximum, input_step = _metric_domain(
            str(row["data_type"]), str(row["unit"]), decimals
        )
        definitions.append(
            MetricDefinition(
                entity_level=level,
                name=name,
                display_name=str(row["display_name"]),
                description=str(row["description"]),
                data_type=str(row["data_type"]),
                unit=str(row["unit"]),
                higher_is=str(row["higher_is"]),
                filterable=bool(row["filterable"]),
                default_visible=bool(row["default_visible"]),
                decimals=decimals,
                numerator=str(row["numerator_definition"]),
                denominator=str(row["denominator_definition"]),
                null_definition=str(row["null_definition"]),
                group=group,
                short_label=str(row["display_name"]),
                valid_minimum=valid_minimum,
                valid_maximum=valid_maximum,
                input_step=input_step,
            )
        )
    return tuple(definitions)


def index_catalog(
    definitions: Iterable[MetricDefinition],
) -> dict[tuple[str, str], MetricDefinition]:
    """Index definitions by entity level and stored column name."""

    materialized = tuple(definitions)
    indexed = {(item.entity_level, item.name): item for item in materialized}
    if len(indexed) != len(materialized):
        raise ValueError("Metric definitions are not unique")
    return indexed
