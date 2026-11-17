"""Versioned, validated conservation-filter presets."""

from __future__ import annotations

import json
import math
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .models import (
    ConservationTarget,
    FilterRecipe,
    MetricConstraint,
    MetricDefinition,
    validate_metric_constraints,
)

CONSERVATION_PLACEHOLDER = "$conservation_frequency"
SUPPORTED_PRESET_SCHEMA = "1.0.0"


@dataclass(frozen=True, slots=True)
class PresetConstraint:
    """One directed threshold before conservation-target resolution."""

    metric: str
    direction: str
    value: float


@dataclass(frozen=True, slots=True)
class FilterPreset:
    """One named, versioned set of ordered filter thresholds."""

    identifier: str
    label: str
    version: str
    description: str
    constraints: tuple[PresetConstraint, ...]

    def recipe(self, target: ConservationTarget) -> FilterRecipe:
        """Resolve the semantic conservation field and create a typed recipe."""

        resolved = []
        for item in self.constraints:
            metric = target.value if item.metric == CONSERVATION_PLACEHOLDER else item.metric
            resolved.append(
                MetricConstraint(
                    metric=metric,
                    minimum=item.value if item.direction == "min" else None,
                    maximum=item.value if item.direction == "max" else None,
                )
            )
        return FilterRecipe(
            constraints=tuple(resolved),
            conservation_target=target,
            preset_id=self.identifier,
            preset_version=self.version,
        )


@dataclass(frozen=True, slots=True)
class PresetCatalog:
    """Validated preset collection in scientific strictness order."""

    schema_version: str
    database_compatibility: str
    scientific_caution: str
    presets: tuple[FilterPreset, ...]

    def get(self, identifier: str) -> FilterPreset:
        """Return a preset by stable identifier."""

        for preset in self.presets:
            if preset.identifier == identifier:
                return preset
        raise KeyError(f"Unknown filter preset: {identifier}")


def _parse_constraint(raw: dict[str, Any]) -> PresetConstraint:
    required = {"metric", "direction", "value"}
    if set(raw) != required:
        raise ValueError(f"Preset constraint fields must be exactly {sorted(required)}")
    direction = str(raw["direction"])
    if direction not in {"min", "max"}:
        raise ValueError("Preset constraint direction must be min or max")
    value = float(raw["value"])
    if not math.isfinite(value):
        raise ValueError("Preset constraint value must be finite")
    return PresetConstraint(str(raw["metric"]), direction, value)


def validate_monotonicity(catalog: PresetCatalog) -> None:
    """Ensure every successive preset is equally or less restrictive."""

    if len(catalog.presets) < 2:
        return
    for stricter, looser in zip(catalog.presets, catalog.presets[1:], strict=False):
        strict_items = {(item.metric, item.direction): item.value for item in stricter.constraints}
        loose_items = {(item.metric, item.direction): item.value for item in looser.constraints}
        if strict_items.keys() != loose_items.keys():
            raise ValueError("All ordered presets must constrain the same metric/direction pairs")
        for key, strict_value in strict_items.items():
            loose_value = loose_items[key]
            direction = key[1]
            valid = (
                strict_value >= loose_value if direction == "min" else strict_value <= loose_value
            )
            if not valid:
                raise ValueError(f"Preset strictness is not monotonic for {key[0]}")


def load_presets(
    path: str | Path,
    filterable_metrics: Iterable[str | MetricDefinition],
) -> PresetCatalog:
    """Load and fully validate a preset JSON file.

    Parameters
    ----------
    path
        JSON document using ``SUPPORTED_PRESET_SCHEMA``.
    filterable_metrics
        Authoritative filterable window definitions or field names from the database registry.
        Full definitions also validate every preset value against its metric domain.
    """

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("preset_schema_version") != SUPPORTED_PRESET_SCHEMA:
        raise ValueError("Unsupported filter preset schema version")
    raw_presets = payload.get("presets")
    order = payload.get("order")
    if not isinstance(raw_presets, dict) or not isinstance(order, list) or not order:
        raise ValueError("Preset JSON requires non-empty presets and order")
    if set(order) != set(raw_presets) or len(order) != len(set(order)):
        raise ValueError("Preset order must contain every preset exactly once")
    supplied_metrics = tuple(filterable_metrics)
    definitions = tuple(item for item in supplied_metrics if isinstance(item, MetricDefinition))
    known = {
        item.name if isinstance(item, MetricDefinition) else str(item) for item in supplied_metrics
    }
    presets = []
    for identifier in order:
        raw = raw_presets[identifier]
        constraints = tuple(_parse_constraint(item) for item in raw.get("constraints", ()))
        metrics = [item.metric for item in constraints]
        if not constraints or len(metrics) != len(set(metrics)):
            raise ValueError(f"Preset {identifier} has empty or duplicate constraints")
        for item in constraints:
            if item.metric != CONSERVATION_PLACEHOLDER and item.metric not in known:
                raise ValueError(f"Unknown filterable preset metric: {item.metric}")
        presets.append(
            FilterPreset(
                identifier=identifier,
                label=str(raw["label"]),
                version=str(raw["version"]),
                description=str(raw["description"]),
                constraints=constraints,
            )
        )
    catalog = PresetCatalog(
        schema_version=str(payload["preset_schema_version"]),
        database_compatibility=str(payload["database_schema_compatibility"]),
        scientific_caution=str(payload["scientific_caution"]),
        presets=tuple(presets),
    )
    validate_monotonicity(catalog)
    if definitions:
        for preset in catalog.presets:
            for target in ConservationTarget:
                validate_metric_constraints(
                    list(preset.recipe(target).constraints),
                    list(definitions),
                )
    return catalog
