"""Typed contracts shared by the metrics application services and user interface."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import pandas as pd


class StratumType(str, Enum):
    """Independent sample groupings stored in the metrics database."""

    GLOBAL = "global"
    CONTINENT = "continent"
    COUNTRY = "country"
    POPULATION = "population"


class Renderer(str, Enum):
    """Supported chart backends."""

    PLOTLY = "plotly"
    MATPLOTLIB = "matplotlib"


class PlotStyle(str, Enum):
    """Supported quantitative track styles."""

    LINE = "line"
    BAR = "bar"


class RegionMode(str, Enum):
    """How peptide windows relate to a submitted coordinate interval."""

    OVERLAP = "overlap"
    CONTAIN = "contain"


class ConservationTarget(str, Enum):
    """Exact-frequency field used to define peptide conservation."""

    TOP_EXACT = "top_exact_peptide_frequency"
    REFERENCE = "reference_peptide_frequency_exact"
    ALLELIC_CONSENSUS = "allelic_consensus_exact_frequency"


@dataclass(frozen=True, slots=True)
class Selection:
    """Describe one complete scientific and presentation selection.

    Parameters
    ----------
    stratum_type, stratum_value
        One independent stored stratum. Global must use the value ``all``.
    block_number
        Positive numeric block identifier.
    window_size
        Positive peptide length stored for the selected block.
    renderer, plot_style
        Presentation choices. They do not alter scientific query results.
    show_points
        Whether line tracks display individual observations.
    """

    stratum_type: StratumType
    stratum_value: str
    block_number: int
    window_size: int
    renderer: Renderer = Renderer.PLOTLY
    plot_style: PlotStyle = PlotStyle.LINE
    show_points: bool = True

    def __post_init__(self) -> None:
        if not self.stratum_value.strip():
            raise ValueError("stratum_value cannot be empty")
        if self.stratum_type is StratumType.GLOBAL and self.stratum_value != "all":
            raise ValueError("The global stratum must use stratum_value='all'")
        if self.block_number < 1:
            raise ValueError("block_number must be positive")
        if self.window_size < 1:
            raise ValueError("window_size must be positive")

    @property
    def data_key(self) -> tuple[str, str, int, int]:
        """Return the cache key that affects scientific data retrieval."""

        return (
            self.stratum_type.value,
            self.stratum_value,
            self.block_number,
            self.window_size,
        )


@dataclass(frozen=True, slots=True)
class StratumOption:
    """One selectable database stratum."""

    stratum_id: int
    stratum_type: StratumType
    value: str
    label: str
    n_samples: int
    small_sample_warning: bool
    is_default: bool = False
    sort_order: int = 0


@dataclass(frozen=True, slots=True)
class BlockOption:
    """Coordinates and reference information for one block."""

    block_id: int
    block_number: int
    name: str
    aa_start: int
    aa_end: int
    length: int
    reference_sequence: str

    def __post_init__(self) -> None:
        if self.aa_end - self.aa_start + 1 != self.length:
            raise ValueError("Block coordinates and length disagree")

    @property
    def label(self) -> str:
        """Return a compact selector label."""

        return f"Block {self.block_number:02d} · AA {self.aa_start}–{self.aa_end}"


@dataclass(frozen=True, slots=True)
class WindowSizeOption:
    """One peptide length available for a selected block."""

    window_size: int
    n_windows: int


@dataclass(frozen=True, slots=True)
class MetricDefinition:
    """Authoritative metric semantics plus application presentation metadata."""

    entity_level: str
    name: str
    display_name: str
    description: str
    data_type: str
    unit: str
    higher_is: str
    filterable: bool
    default_visible: bool
    decimals: int
    numerator: str
    denominator: str
    null_definition: str
    group: str
    short_label: str
    valid_minimum: float | None
    valid_maximum: float | None
    input_step: float

    @property
    def help_text(self) -> str:
        """Return concise UI help without losing denominator and null semantics."""

        return (
            f"{self.description} Denominator: {self.denominator}. "
            f"Unavailable when: {self.null_definition}."
        )


@dataclass(frozen=True, slots=True)
class MetricConstraint:
    """One inclusive numeric filter applied to a stored window metric."""

    metric: str
    minimum: float | None = None
    maximum: float | None = None

    def __post_init__(self) -> None:
        if not self.metric.strip():
            raise ValueError("metric cannot be empty")
        if self.minimum is None and self.maximum is None:
            raise ValueError("A constraint requires a minimum, maximum, or both")
        for label, value in (("minimum", self.minimum), ("maximum", self.maximum)):
            if value is not None and not math.isfinite(float(value)):
                raise ValueError(f"{label} must be finite")
        if self.minimum is not None and self.maximum is not None:
            if self.minimum > self.maximum:
                raise ValueError("minimum cannot exceed maximum")


def validate_metric_constraints(
    constraints: Sequence[MetricConstraint],
    definitions: Sequence[MetricDefinition],
) -> None:
    """Validate constraints against authoritative metric domains.

    Bounds are determined by the metric datatype and unit when the catalog is built.
    """

    catalog = {item.name: item for item in definitions if item.filterable}
    for constraint in constraints:
        try:
            definition = catalog[constraint.metric]
        except KeyError as error:
            raise ValueError(f"Unknown filterable window metric: {constraint.metric}") from error
        maximum = definition.valid_maximum
        for label, value in (
            ("minimum", constraint.minimum),
            ("maximum", constraint.maximum),
        ):
            if value is None:
                continue
            if definition.valid_minimum is not None and value < definition.valid_minimum:
                raise ValueError(
                    f"{definition.display_name} {label} must be at least "
                    f"{definition.valid_minimum:g}"
                )
            if maximum is not None and value > maximum:
                raise ValueError(f"{definition.display_name} {label} must be at most {maximum:g}")


@dataclass(frozen=True, slots=True)
class FilterRecipe:
    """Serializable, ordered peptide-filter specification."""

    constraints: tuple[MetricConstraint, ...] = ()
    conservation_target: ConservationTarget = ConservationTarget.TOP_EXACT
    preset_id: str | None = None
    preset_version: str | None = None
    null_policy: str = "fail_active_constraint"

    def __post_init__(self) -> None:
        names = [item.metric for item in self.constraints]
        if len(names) != len(set(names)):
            raise ValueError("A filter recipe cannot contain duplicate metric constraints")
        if self.null_policy != "fail_active_constraint":
            raise ValueError("Unsupported null policy")


@dataclass(slots=True)
class FilterResult:
    """Rows and accounting generated by a submitted filter recipe."""

    rows: pd.DataFrame
    initial_count: int
    retained_count: int
    steps: pd.DataFrame
    recipe: FilterRecipe
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class RegionQuery:
    """One inclusive coordinate query inside the selected block."""

    start: int
    end: int
    mode: RegionMode = RegionMode.OVERLAP

    def __post_init__(self) -> None:
        if self.start < 1 or self.end < 1:
            raise ValueError("Region coordinates must be positive")
        if self.start > self.end:
            raise ValueError("Region start cannot exceed region end")


@dataclass(slots=True)
class RegionResult:
    """Rows selected by a submitted coordinate query."""

    query: RegionQuery
    rows: pd.DataFrame
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class DatabaseIdentity:
    """Immutable identity and provenance for one database artifact."""

    path: str
    size_bytes: int
    schema_version: str
    sqlite_user_version: int | None
    build_timestamp: str
    build_status: str
    package_version: str
    gene_identifier: str
    gene_name: str | None
    artifact_sha256: str | None = None
    source_files: tuple[dict[str, Any], ...] = field(default_factory=tuple)
    storage_backend: str = "sqlite"
    artifact_root: str | None = None
    cache_identity: str | None = None
