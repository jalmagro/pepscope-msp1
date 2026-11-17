"""Deterministic in-memory peptide filtering and recipe serialization."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import asdict
from datetime import datetime, timezone
from typing import Any

import pandas as pd

from .models import (
    FilterRecipe,
    FilterResult,
    MetricConstraint,
    MetricDefinition,
    validate_metric_constraints,
)


def apply_filters(
    frame: pd.DataFrame,
    recipe: FilterRecipe,
    *,
    filterable_metrics: Iterable[str] | None = None,
    metric_definitions: Iterable[MetricDefinition] | None = None,
) -> FilterResult:
    """Apply ordered inclusive constraints with explicit null-failure accounting."""

    definitions = tuple(metric_definitions or ())
    if definitions:
        validate_metric_constraints(
            list(recipe.constraints),
            list(definitions),
        )
    allowed = set(frame.columns if filterable_metrics is None else filterable_metrics)
    working = frame.copy()
    steps: list[dict[str, Any]] = []
    initial = len(working)
    for number, constraint in enumerate(recipe.constraints, start=1):
        if constraint.metric not in allowed or constraint.metric not in working.columns:
            raise ValueError(f"Unknown filterable window metric: {constraint.metric}")
        values = pd.to_numeric(working[constraint.metric], errors="coerce")
        passes = values.notna()
        if constraint.minimum is not None:
            passes &= values.ge(constraint.minimum)
        if constraint.maximum is not None:
            passes &= values.le(constraint.maximum)
        before = len(working)
        null_failures = int(values.isna().sum())
        working = working.loc[passes].copy()
        steps.append(
            {
                "step": number,
                "metric": constraint.metric,
                "minimum": constraint.minimum,
                "maximum": constraint.maximum,
                "input_count": before,
                "null_failures": null_failures,
                "removed_at_step": before - len(working),
                "remaining": len(working),
                "remaining_fraction": len(working) / initial if initial else 0.0,
            }
        )
    warnings = ()
    if not recipe.constraints:
        warnings = ("No active constraints: every selected peptide is retained.",)
    return FilterResult(
        rows=working.reset_index(drop=True),
        initial_count=initial,
        retained_count=len(working),
        steps=pd.DataFrame(steps),
        recipe=recipe,
        warnings=warnings,
    )


def recipe_manifest(
    recipe: FilterRecipe,
    *,
    selection: Mapping[str, Any],
    database_identity: Mapping[str, Any],
    application_version: str,
    created_at: str | None = None,
) -> dict[str, Any]:
    """Return a JSON-safe reproducibility record for one submitted recipe."""

    return {
        "created_at": created_at or datetime.now(timezone.utc).isoformat(),
        "application_version": application_version,
        "database_identity": dict(database_identity),
        "selection": dict(selection),
        "preset_id": recipe.preset_id,
        "preset_version": recipe.preset_version,
        "conservation_target": recipe.conservation_target.value,
        "null_policy": recipe.null_policy,
        "constraints": [asdict(item) for item in recipe.constraints],
    }


def recipe_from_manifest(payload: Mapping[str, Any]) -> FilterRecipe:
    """Reconstruct the scientific recipe from a serialized manifest."""

    from .models import ConservationTarget

    return FilterRecipe(
        constraints=tuple(MetricConstraint(**item) for item in payload["constraints"]),
        conservation_target=ConservationTarget(str(payload["conservation_target"])),
        preset_id=payload.get("preset_id"),
        preset_version=payload.get("preset_version"),
        null_policy=str(payload.get("null_policy", "fail_active_constraint")),
    )
