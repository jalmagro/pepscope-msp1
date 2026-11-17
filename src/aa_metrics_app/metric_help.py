"""Versioned, application-owned guidance for database metrics."""

from __future__ import annotations

import json
from collections.abc import Iterable
from html import escape
from pathlib import Path

import pandas as pd

from .models import MetricDefinition

METRIC_HELP_SCHEMA_VERSION = "1.1.0"
METRIC_HELP_LEVELS = frozenset({"block", "locus", "window"})
METRIC_HELP_FIELDS = (
    "entity_level",
    "key",
    "display_name",
    "group",
    "value_type",
    "unit",
    "plain_language",
    "technical_definition",
    "interpretation",
    "example",
    "example_interpretation",
    "help_tip",
)


def _expected_value_type(definition: MetricDefinition) -> str:
    if definition.unit == "boolean":
        return "boolean"
    return "integer" if definition.data_type.upper() == "INTEGER" else "number"


def load_metric_help(path: str | Path) -> pd.DataFrame:
    """Load and structurally validate a metric-help JSON artifact.

    Parameters
    ----------
    path
        JSON file produced by ``scripts/build_metric_help.py``.

    Returns
    -------
    pandas.DataFrame
        One ordered row per ``(entity_level, key)`` metric.
    """

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema_version") != METRIC_HELP_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported metric-help schema version: "
            f"{payload.get('schema_version')!r}; expected {METRIC_HELP_SCHEMA_VERSION!r}"
        )
    database_schema_version = payload.get("database_schema_version")
    if not isinstance(database_schema_version, str) or not database_schema_version.strip():
        raise ValueError("Metric-help JSON must declare database_schema_version")
    metrics = payload.get("metrics")
    if not isinstance(metrics, list) or not metrics:
        raise ValueError("Metric-help JSON must contain a non-empty metrics list")
    frame = pd.DataFrame.from_records(metrics)
    missing = set(METRIC_HELP_FIELDS).difference(frame.columns)
    if missing:
        raise ValueError(f"Metric-help JSON is missing fields: {sorted(missing)}")
    frame = frame.loc[:, METRIC_HELP_FIELDS].copy()
    keys = frame[["entity_level", "key"]].astype(str).agg("/".join, axis=1)
    if keys.duplicated().any():
        duplicates = sorted(keys.loc[keys.duplicated(keep=False)].unique())
        raise ValueError(f"Metric-help JSON contains duplicate keys: {duplicates}")
    unknown_levels = sorted(set(frame.entity_level).difference(METRIC_HELP_LEVELS))
    if unknown_levels:
        raise ValueError(f"Metric-help JSON contains unknown levels: {unknown_levels}")
    optional_text = {"technical_definition", "example", "example_interpretation"}
    required_text = [field for field in METRIC_HELP_FIELDS if field not in optional_text]
    for field in required_text:
        invalid = frame[field].isna() | frame[field].astype(str).str.strip().eq("")
        if invalid.any():
            bad_keys = keys.loc[invalid].tolist()
            raise ValueError(f"Metric-help field {field!r} is empty for: {bad_keys}")
    example_missing = frame["example"].isna() | frame["example"].astype(str).str.strip().eq("")
    interpretation_missing = frame["example_interpretation"].isna() | (
        frame["example_interpretation"].astype(str).str.strip().eq("")
    )
    if not example_missing.equals(interpretation_missing):
        bad_keys = keys.loc[example_missing.ne(interpretation_missing)].tolist()
        raise ValueError(
            "Metric-help examples and example interpretations must be provided together for: "
            f"{bad_keys}"
        )
    bookkeeping = frame["help_tip"].str.contains(
        r"\b(?:numerator|denominator)\s*:", case=False, regex=True, na=False
    )
    if bookkeeping.any():
        bad_keys = keys.loc[bookkeeping].tolist()
        raise ValueError(
            "Metric-help tips must explain metrics without numerator/denominator labels for: "
            f"{bad_keys}"
        )
    frame = frame.sort_values(
        ["entity_level", "group", "display_name"], kind="stable"
    ).reset_index(drop=True)
    frame.attrs["schema_version"] = METRIC_HELP_SCHEMA_VERSION
    frame.attrs["database_schema_version"] = database_schema_version
    return frame


def validate_metric_help(
    frame: pd.DataFrame, definitions: Iterable[MetricDefinition]
) -> None:
    """Require exact key, type, and unit agreement with a database registry.

    Display labels intentionally belong to the application catalog and may translate stable
    database terminology into clearer user-facing language.
    """

    expected = {(item.entity_level, item.name): item for item in definitions}
    actual = set(zip(frame.entity_level, frame.key, strict=True))
    expected_keys = set(expected)
    missing = sorted(expected_keys - actual)
    unknown = sorted(actual - expected_keys)
    if missing or unknown:
        raise ValueError(
            "Metric-help/database registry mismatch. "
            f"Missing help: {missing}; unknown help: {unknown}"
        )
    indexed = frame.set_index(["entity_level", "key"])
    mismatches: list[str] = []
    for key, definition in expected.items():
        row = indexed.loc[key]
        checks = {
            "unit": definition.unit,
            "value_type": _expected_value_type(definition),
        }
        for field, expected_value in checks.items():
            if str(row[field]) != str(expected_value):
                mismatches.append(
                    f"{key[0]}/{key[1]} {field}: {row[field]!r} != {expected_value!r}"
                )
    if mismatches:
        raise ValueError("Metric-help/database metadata mismatch: " + "; ".join(mismatches))


def metric_help_for_level(frame: pd.DataFrame, level: str) -> pd.DataFrame:
    """Return a copy of the guidance for one validated entity level."""

    if level not in METRIC_HELP_LEVELS:
        raise ValueError(f"Unknown metric-help level: {level!r}")
    selected = frame.loc[frame.entity_level.eq(level)].copy()
    if selected.empty:
        raise ValueError(f"Metric-help catalog contains no {level!r} metrics")
    return selected.reset_index(drop=True)


def _unit_context(value_type: str, unit: str) -> str:
    """Return approachable type/unit context for one narrative metric entry."""

    if value_type == "boolean" or unit == "boolean":
        return "True/false indicator"
    labels = {
        "fraction": "Fraction, usually displayed as a percentage",
        "count": "Count",
        "samples": "Sample count",
        "sites": "Amino-acid position count",
        "AA": "Amino-acid distance",
        "nats": "Diversity measured in nats",
        "effective_calls": "Effective callable-sample weight",
        "index": "Unitless index",
        "rate": "Rate",
        "proportion": "Proportion",
    }
    return labels.get(unit, f"{value_type.title()} value ({unit})")


def metric_guide_html(frame: pd.DataFrame) -> str:
    """Render validated metric guidance as escaped, grouped semantic HTML.

    Parameters
    ----------
    frame
        Rows returned by :func:`metric_help_for_level`.

    Returns
    -------
    str
        Semantic HTML containing one section per group and one article per metric. Every catalog
        value is HTML-escaped before interpolation.
    """

    required = set(METRIC_HELP_FIELDS)
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Metric guide is missing catalog fields: {missing}")
    sections: list[str] = ["<div class='metric-guide'>"]
    for group, group_frame in frame.groupby("group", sort=False, dropna=False):
        sections.append("<section class='metric-guide-group'>")
        sections.append(f"<h3>{escape(str(group).replace('_', ' ').title())}</h3>")
        for row in group_frame.itertuples(index=False):
            context = _unit_context(str(row.value_type), str(row.unit))
            sections.extend(
                (
                    "<article class='metric-guide-entry'>",
                    f"<h4>{escape(str(row.display_name))}</h4>",
                    f"<p class='metric-guide-meta'>{escape(context)}</p>",
                    f"<p>{escape(str(row.plain_language))}</p>",
                )
            )
            if pd.notna(row.technical_definition) and str(row.technical_definition).strip():
                sections.append(
                    "<p><strong>Technical note.</strong> "
                    f"{escape(str(row.technical_definition))}</p>"
                )
            sections.append(
                f"<p><strong>How to read it.</strong> {escape(str(row.interpretation))}</p>"
            )
            if pd.notna(row.example) and str(row.example).strip():
                sections.append("<div class='metric-guide-example'>")
                sections.append(
                    f"<p><strong>Example.</strong> {escape(str(row.example))}</p>"
                )
                sections.append(
                    "<p><strong>What the example tells us.</strong> "
                    f"{escape(str(row.example_interpretation))}</p>"
                )
                sections.append("</div>")
            sections.append("</article>")
        sections.append("</section>")
    sections.append("</div>")
    return "".join(sections)
