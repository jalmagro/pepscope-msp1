"""MSP1 Peptide Conservation Explorer Streamlit entry point."""

from __future__ import annotations

import json
import logging
import math
import os
import sys
from dataclasses import replace
from html import escape
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

APP_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = APP_ROOT.parent
SOURCE_ROOT = PROJECT_ROOT / "src"
if SOURCE_ROOT.is_dir() and str(SOURCE_ROOT) not in sys.path:
    # Streamlit adds ``app/`` rather than the repository root to sys.path when the package has not
    # been installed. Keep the entry point runnable from a fresh checkout without PYTHONPATH.
    sys.path.insert(0, str(SOURCE_ROOT))

from aa_metrics_app import (  # noqa: E402
    AppContent,
    BlockOption,
    ConservationTarget,
    DataBackend,
    DatabaseCompatibilityError,
    DatabaseIdentity,
    DataSourceConfig,
    DuckDBRuntimeConfig,
    FilterRecipe,
    MetricConstraint,
    MetricDefinition,
    MetricsRepository,
    PlotStyle,
    RegionMode,
    RegionQuery,
    Renderer,
    Selection,
    StratumOption,
    StratumType,
    WindowSizeOption,
    apply_filters,
    build_filter_bundle,
    build_mutation_bundle,
    build_region_bundle,
    canonical_filter_rows,
    data_source_from_environment,
    load_app_content,
    load_data_source_config,
    load_metric_help,
    load_presets,
    metric_guide_html,
    metric_help_for_level,
    page_kicker_html,
    page_lead_html,
    page_title_html,
    recipe_manifest,
    safe_csv,
    section_intro_html,
    select_region,
    validate_metric_constraints,
    validate_metric_help,
)
from aa_metrics_app.charts import (  # noqa: E402
    PLOTLY_CONFIG,
    PRIMARY_TOMATO,
    TrackSpec,
    prepare_locus_hover_data,
)
from aa_metrics_app.formatting import (  # noqa: E402
    format_bytes,
    format_count,
    format_fraction,
    format_number,
)
from aa_metrics_app.sequences import (  # noqa: E402
    align_sequences,
    alignment_html,
    available_sequence_rank_columns,
)

LOGGER = logging.getLogger(__name__)
DEFAULT_DATABASE = PROJECT_ROOT / "results" / "database" / "msp1-africa-metrics-v2.sqlite"
DATABASE_ENV = "AA_METRICS_DB"
DATABASE_SHA_ENV = "AA_METRICS_DB_SHA256"
DATABASE_SIZE_ENV = "AA_METRICS_DB_SIZE"
SOURCE_CONFIG_ENV = "AA_METRICS_SOURCE_CONFIG"
DATA_SOURCE_CONFIG_PATH = APP_ROOT / "config" / "data-source.json"
PRESET_PATH = APP_ROOT / "config" / "filter-presets.json"
METRIC_HELP_PATH = APP_ROOT / "config" / "metric-help.json"
APP_CONTENT_PATH = APP_ROOT / "config" / "app-content.json"
MUTATION_HELP_PATH = APP_ROOT / "config" / "mutation-help.json"
LOGO_HELP_PATH = APP_ROOT / "config" / "logo-help.json"
STRUCTURE_HELP_PATH = APP_ROOT / "config" / "structure-help.json"
PAGE_HEADER_ICON_WIDTH = 72
PAGE_HEADER_COLUMN_WEIGHTS = (0.75, 11.25)
PAGE_HEADER_COLUMN_GAP = None

NAVIGATION = (
    "Block overview",
    "Locus tracks",
    "Window tracks",
    "Explore a region",
    "Filter & export",
    "About & limitations",
)

FRIENDLY_COLUMN_LABELS = {
    "reference_sequence": "3D7 sequence",
    "allelic_consensus_sequence": "Canonical consensus sequence",
    "all_call_consensus_sequence": "All-call consensus sequence",
    "top_exact_sequence": "Top complete haplotype",
    "top_exact_peptide_sequence": "Top complete peptide",
    "modal_all_call_sequence": "Top sequence (all)",
    "n_exact_haplotype_observations": "Complete haplotype observations",
    "n_exact_peptide_observations": "Complete peptide observations",
    "n_non_exact_peptide_observations": "Incomplete peptide observations",
    "n_unique_exact_haplotypes": "Unique complete haplotypes",
    "n_unique_exact_peptides": "Unique complete peptides",
    "top_exact_sequence_frequency": "Top complete haplotype frequency",
    "top_exact_peptide_frequency": "Top complete peptide frequency",
    "modal_all_call_sequence_frequency": "Top sequence frequency (all)",
    "reference_haplotype_frequency_exact": "3D7 haplotype frequency (complete)",
    "reference_peptide_frequency_exact": "3D7 peptide frequency (complete)",
    "allelic_consensus_exact_frequency": "Consensus frequency (complete)",
    "all_call_consensus_sequence_frequency": "All-call consensus frequency (all)",
}

# Compact labels are restricted to fixed-width sequence boxes. Full scientific names remain in
# cards, tabs, controls, prose, and exports.
TOP_HAP_LABEL = "Top hap"
TOP_HAP_ALL_LABEL = "Top hap (all)"
TOP_COMPLETE_PEP_LABEL = "Top complete pep"
TOP_PEPTIDE_ALL_LABEL = "Top peptide (all)"


def _plotly_chart(container: Any, figure: Any, *, key: str | None = None) -> Any:
    """Render Plotly with arguments supported by Streamlit 1.50 and newer.

    Both supported Streamlit versions already default to full container width. Passing the newer
    ``width`` argument to Streamlit 1.50 is interpreted as a deprecated arbitrary Plotly keyword,
    so this boundary deliberately sends only the shared Plotly ``config`` and an optional key.
    """

    options: dict[str, Any] = {"config": PLOTLY_CONFIG}
    if key is not None:
        options["key"] = key
    return container.plotly_chart(figure, **options)


def _column_config(columns: Any, registry: pd.DataFrame | None = None) -> dict[str, str]:
    """Return friendly display labels while preserving stored DataFrame column names."""

    labels = {name: label for name, label in FRIENDLY_COLUMN_LABELS.items() if name in columns}
    if registry is not None:
        labels.update(
            {
                str(row["name"]): str(row["display_name"])
                for row in registry.to_dict("records")
                if row["name"] in columns
            }
        )
    return labels


def _source_config_from_token(token: str) -> DataSourceConfig:
    """Decode one app-owned, backend-neutral cache token."""

    payload = json.loads(token)
    duckdb = payload["duckdb"]
    return DataSourceConfig(
        backend=DataBackend(payload["backend"]),
        artifact=Path(payload["artifact"]),
        artifact_root=(
            Path(payload["artifact_root"]) if payload.get("artifact_root") is not None else None
        ),
        expected_size_bytes=payload.get("expected_size_bytes"),
        expected_sha256=payload.get("expected_sha256"),
        expected_total_size_bytes=payload.get("expected_total_size_bytes"),
        expected_dataset_sha256=payload.get("expected_dataset_sha256"),
        duckdb=DuckDBRuntimeConfig(**duckdb),
    )


def _repository(token: str, *, validate: bool = False) -> MetricsRepository:
    """Open a short-lived repository from the deployment source token."""

    return MetricsRepository(_source_config_from_token(token), validate=validate)


def _data_source_token() -> str:
    """Resolve deployment configuration and include backend-aware cache identity."""

    if os.environ.get(SOURCE_CONFIG_ENV):
        source = load_data_source_config(os.environ[SOURCE_CONFIG_ENV])
    elif any(os.environ.get(name) for name in (DATABASE_ENV, DATABASE_SIZE_ENV, DATABASE_SHA_ENV)):
        source = data_source_from_environment(DEFAULT_DATABASE)
    elif DATA_SOURCE_CONFIG_PATH.is_file():
        source = load_data_source_config(DATA_SOURCE_CONFIG_PATH)
    else:
        source = data_source_from_environment(DEFAULT_DATABASE)
    cache_identity = MetricsRepository(source, validate=False).cache_identity
    payload = {
        "backend": source.backend.value,
        "artifact": str(source.artifact),
        "artifact_root": (str(source.artifact_root) if source.artifact_root is not None else None),
        "expected_size_bytes": source.expected_size_bytes,
        "expected_sha256": source.expected_sha256,
        "expected_total_size_bytes": source.expected_total_size_bytes,
        "expected_dataset_sha256": source.expected_dataset_sha256,
        "duckdb": {
            "memory_limit": source.duckdb.memory_limit,
            "threads": source.duckdb.threads,
            "temp_directory": source.duckdb.temp_directory,
            "max_temp_directory_size": source.duckdb.max_temp_directory_size,
        },
        "cache_identity": cache_identity,
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


@st.cache_data(show_spinner=False, max_entries=16)
def load_identity(source_token: str) -> DatabaseIdentity:
    """Validate the selected source once per backend-aware cache identity."""

    source = _source_config_from_token(source_token)
    repository = _repository(source_token, validate=True)
    identity = repository.identity(
        include_sha256=source.backend is DataBackend.PARQUET_DUCKDB
    )
    if source.expected_sha256 is not None:
        identity = replace(identity, artifact_sha256=source.expected_sha256.lower())
    return identity


@st.cache_data(show_spinner=False, max_entries=16)
def load_strata(path: str) -> tuple[StratumOption, ...]:
    """Load all independent stratum choices."""

    return _repository(path).strata()


@st.cache_data(show_spinner=False, max_entries=16)
def load_blocks(path: str) -> tuple[BlockOption, ...]:
    """Load ordered block choices."""

    return _repository(path).blocks()


@st.cache_data(show_spinner=False, max_entries=16)
def load_gene_metadata(path: str) -> pd.DataFrame:
    """Load decoded gene metadata for the About page."""

    return _repository(path).genes()


@st.cache_data(show_spinner=False, max_entries=16)
def load_database_statistics(path: str) -> pd.DataFrame:
    """Load aggregate database inventory counts."""

    return _repository(path).database_statistics()


@st.cache_data(show_spinner=False, max_entries=16)
def load_peptide_rows_by_stratum(path: str) -> pd.DataFrame:
    """Load stratum-specific peptide metric row counts."""

    return _repository(path).peptide_rows_by_stratum()


@st.cache_data(show_spinner=False, max_entries=16)
def load_window_sizes(path: str, block_number: int) -> tuple[WindowSizeOption, ...]:
    """Load peptide lengths available for one selected block."""

    return _repository(path).window_sizes(block_number)


@st.cache_data(show_spinner=False, max_entries=16)
def load_block_summary(
    path: str, stratum_type: str, stratum_value: str, block_number: int
) -> pd.DataFrame:
    """Load one selected block metric row."""

    return _repository(path).block_summary(stratum_type, stratum_value, block_number)


@st.cache_data(show_spinner=False, max_entries=16)
def load_block_histogram(
    path: str, stratum_type: str, stratum_value: str, block_number: int
) -> pd.DataFrame:
    """Load stored block callability bins."""

    return _repository(path).block_callability_histogram(stratum_type, stratum_value, block_number)


@st.cache_data(show_spinner=False, max_entries=16)
def load_block_sequences(
    path: str,
    stratum_type: str,
    stratum_value: str,
    block_number: int,
    scope: str,
) -> pd.DataFrame:
    """Load complete or all-sample ranked block sequences."""

    return _repository(path).block_sequence_frequencies(
        stratum_type, stratum_value, block_number, scope=scope
    )


@st.cache_data(show_spinner=False, max_entries=16)
def load_metric_definitions(path: str, level: str) -> tuple[MetricDefinition, ...]:
    """Load typed metric definitions for one entity level."""

    return _repository(path).metric_definitions(level)


@st.cache_data(show_spinner=False, max_entries=16)
def load_metric_registry(path: str, level: str) -> pd.DataFrame:
    """Load serializable metric definitions for one entity level."""

    definitions = load_metric_definitions(path, level)
    registry = pd.DataFrame(
        [
            {
                "name": item.name,
                "display_name": item.display_name,
                "description": item.description,
                "unit": item.unit,
                "decimals": item.decimals,
                "numerator": item.numerator,
                "denominator": item.denominator,
                "null_definition": item.null_definition,
                "group": item.group,
                "filterable": item.filterable,
                "higher_is": item.higher_is,
                "short_label": item.short_label,
                "data_type": item.data_type,
                "valid_minimum": item.valid_minimum,
                "valid_maximum": item.valid_maximum,
                "input_step": item.input_step,
            }
            for item in definitions
        ]
    )
    guidance_columns = [
        "key",
        "display_name",
        "plain_language",
        "technical_definition",
        "interpretation",
        "example",
        "example_interpretation",
        "help_tip",
    ]
    guidance = metric_help_for_level(
        load_metric_help_catalog(str(METRIC_HELP_PATH), path), level
    ).loc[:, guidance_columns]
    merged = (
        registry.drop(columns="display_name")
        .merge(guidance, left_on="name", right_on="key", validate="one_to_one")
        .drop(columns="key")
    )
    # Plot hover text and controls share the same approachable JSON-backed explanation.
    merged["description"] = merged["plain_language"]
    return merged


@st.cache_data(show_spinner=False, max_entries=16)
def load_metric_help_catalog(help_path: str, database_path: str) -> pd.DataFrame:
    """Load application guidance and require exact agreement with the database registry."""

    frame = load_metric_help(help_path)
    repository = _repository(database_path)
    identity = repository.identity()
    if frame.attrs.get("database_schema_version") != identity.schema_version:
        raise ValueError("Metric-help database schema does not match the selected metrics database")
    validate_metric_help(frame, repository.metric_definitions())
    return frame


@st.cache_data(show_spinner=False, max_entries=16)
def load_app_content_catalog(content_path: str) -> AppContent:
    """Load validated, editable page and section copy."""

    return load_app_content(content_path)


@st.cache_data(show_spinner=False, max_entries=16)
def load_locus_track(
    path: str, stratum_type: str, stratum_value: str, block_number: int
) -> pd.DataFrame:
    """Load one block-bounded locus metric frame."""

    return _repository(path).locus_track(stratum_type, stratum_value, block_number)


@st.cache_data(show_spinner=False, max_entries=16)
def load_window_track(
    path: str,
    stratum_type: str,
    stratum_value: str,
    block_number: int,
    window_size: int,
) -> pd.DataFrame:
    """Load one selected stratum/block/window-size metric frame."""

    return _repository(path).window_track(stratum_type, stratum_value, block_number, window_size)


@st.cache_data(show_spinner=False, max_entries=16)
def load_mutation_facts(path: str, aa_start: int, aa_end: int) -> pd.DataFrame:
    """Load tidy dense mutation facts for an inclusive interval."""

    return _repository(path).mutations(aa_start, aa_end)


@st.cache_data(show_spinner=False, max_entries=16)
def load_mutation_comparison(
    path: str,
    aa_start: int,
    aa_end: int,
    stratum_type: str,
    stratum_value: str,
    include_missing: bool,
) -> pd.DataFrame:
    """Load a selected-group-first wide mutation comparison."""

    return _repository(path).mutation_comparison(
        aa_start,
        aa_end,
        stratum_type,
        stratum_value,
        include_missing=include_missing,
    )


@st.cache_data(show_spinner=False, max_entries=16)
def load_mutation_help(path: str) -> dict[str, str]:
    """Load the versioned mutation explanation catalogue."""

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    required = {
        "schema_version",
        "title",
        "intro",
        "frequency",
        "missing",
        "canonical_id",
        "warning",
    }
    if set(payload) != required or payload["schema_version"] != "1.0.0":
        raise ValueError("Mutation-help JSON does not match schema 1.0.0")
    if any(not isinstance(value, str) or not value.strip() for value in payload.values()):
        raise ValueError("Mutation-help JSON values must be non-empty text")
    return payload


@st.cache_data(show_spinner=False, max_entries=16)
def load_allele_spectrum(
    path: str,
    stratum_type: str,
    stratum_value: str,
    aa_start: int,
    aa_end: int,
) -> pd.DataFrame:
    """Load sample-normalized residue weights for a logo interval."""

    return _repository(path).locus_allele_spectrum_interval(
        stratum_type, stratum_value, aa_start, aa_end
    )


@st.cache_data(show_spinner=False, max_entries=16)
def load_logo_help(path: str) -> dict[str, str]:
    """Load and minimally validate the logo explanation catalogue."""

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    required = {
        "schema_version",
        "title",
        "intro",
        "mixed",
        "uncalled",
        "consensus",
        "interpretation",
    }
    if set(payload) != required or payload["schema_version"] != "1.0.0":
        raise ValueError("Logo-help JSON does not match schema 1.0.0")
    if any(not isinstance(value, str) or not value.strip() for value in payload.values()):
        raise ValueError("Logo-help JSON values must be non-empty text")
    return payload


@st.cache_data(show_spinner=False, max_entries=16)
def load_structure_candidates(path: str, aa_start: int, aa_end: int) -> pd.DataFrame:
    """Load deterministic reviewed structure candidates for an interval."""

    return _repository(path).structures_for_interval(aa_start, aa_end)


@st.cache_data(show_spinner=False, max_entries=16)
def load_structure_map(
    path: str,
    structure_chain_id: int,
    aa_start: int,
    aa_end: int,
) -> pd.DataFrame:
    """Load one reviewed project-to-structure interval map."""

    return _repository(path).structure_residue_map(structure_chain_id, aa_start, aa_end)


@st.cache_data(show_spinner=False, max_entries=16)
def load_structure_backbone(
    asset_path: str,
    expected_sha256: str,
    chain_id: str,
    residue_map: pd.DataFrame,
) -> pd.DataFrame:
    """Verify and parse a local model-1 C-alpha backbone."""

    from aa_metrics_app.structures import backbone_from_mmcif, verify_structure_asset

    asset = verify_structure_asset(asset_path, expected_sha256)
    return backbone_from_mmcif(asset, chain_id, residue_map)


@st.cache_data(show_spinner=False, max_entries=16)
def load_structure_help(path: str) -> dict[str, str]:
    """Load and minimally validate structure-view guidance."""

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    required = {
        "schema_version",
        "title",
        "intro",
        "selection",
        "interpretation",
        "confidence",
        "fallback",
    }
    if set(payload) != required or payload["schema_version"] != "1.0.0":
        raise ValueError("Structure-help JSON does not match schema 1.0.0")
    if any(not isinstance(value, str) or not value.strip() for value in payload.values()):
        raise ValueError("Structure-help JSON values must be non-empty text")
    return payload


@st.cache_data(show_spinner=False, max_entries=16)
def load_filter_presets(path: str):
    """Load presets after validating them against the database metric registry."""

    definitions = load_metric_definitions(path, "window")
    return load_presets(PRESET_PATH, (item for item in definitions if item.filterable))


def _load_css() -> None:
    css_path = APP_ROOT / "assets" / "app.css"
    st.markdown(f"<style>{css_path.read_text(encoding='utf-8')}</style>", unsafe_allow_html=True)


def _stratum_selection(
    strata: tuple[StratumOption, ...],
) -> tuple[StratumType, StratumOption]:
    labels = {
        StratumType.GLOBAL: "All samples",
        StratumType.CONTINENT: "Continent",
        StratumType.POPULATION: "Population",
        StratumType.COUNTRY: "Country",
    }
    available_types = tuple(dict.fromkeys(item.stratum_type for item in strata))
    default = next(item for item in strata if item.is_default)
    type_labels = [labels[item] for item in available_types]
    label = st.sidebar.selectbox(
        "Cohort level",
        type_labels,
        index=available_types.index(default.stratum_type),
        key="cohort_level",
        help=(
            "Choose a set of samples. Country and Population are alternative views; "
            "samples from a country can belong to more than one population."
        ),
    )
    selected_type = next(item for item, text in labels.items() if text == label)
    options = [item for item in strata if item.stratum_type is selected_type]
    selected_label = st.sidebar.selectbox(
        label,
        [item.label for item in options],
        key=f"{selected_type.value}_value",
        help=f"Metrics will use only QC-pass samples in the selected {label.lower()}.",
    )
    selected = next(item for item in options if item.label == selected_label)
    return selected_type, selected


def _reset_display_state() -> None:
    """Restore display-only widget state before Streamlit creates the next page run."""

    st.session_state.interactive_plots = True
    st.session_state.track_style = "Line"
    st.session_state.show_points = True


def _sidebar(
    path: str,
    strata: tuple[StratumOption, ...],
    blocks: tuple[BlockOption, ...],
    content: AppContent,
) -> tuple[Selection, StratumOption, BlockOption, str]:
    st.sidebar.caption(content.sidebar_intro)

    st.sidebar.markdown("### Navigation")
    navigation = st.sidebar.selectbox(
        "Go to",
        NAVIGATION,
        key="navigation",
        label_visibility="collapsed",
        help="Choose the analysis or exploration page to display.",
    )

    st.sidebar.markdown("### Cohort")
    stratum_type, stratum = _stratum_selection(strata)
    st.sidebar.caption(f"{format_count(stratum.n_samples)} QC-pass samples")
    if stratum.small_sample_warning:
        st.sidebar.warning(
            "Small stratum: estimates may be unstable. Check the number of complete observations.",
            icon="⚠️",
        )

    st.sidebar.markdown("### Region")
    block_number = st.sidebar.selectbox(
        "Block",
        [item.block_number for item in blocks],
        format_func=lambda value: next(item.label for item in blocks if item.block_number == value),
        key="block_number",
        help="MSP1 blocks are one-based inclusive coordinate intervals.",
    )
    block = next(item for item in blocks if item.block_number == block_number)
    size_options = load_window_sizes(path, block_number)
    if not size_options:
        raise DatabaseCompatibilityError("The selected block has no stored peptide windows")
    size = st.sidebar.selectbox(
        "Peptide length",
        [item.window_size for item in size_options],
        index=min(8, len(size_options) - 1),
        format_func=lambda value: f"{value} AA",
        key="window_size",
        help="Sliding windows advance one amino-acid coordinate and never cross block boundaries.",
    )

    st.sidebar.markdown("### Display")
    interactive = st.sidebar.toggle(
        "Interactive plots (Plotly)",
        value=True,
        key="interactive_plots",
        help="Turn off to render static Matplotlib plots with the same values and scales.",
    )
    style_label = st.sidebar.radio(
        "Track style",
        ("Line", "Bar"),
        horizontal=True,
        key="track_style",
        help="Binary tracks remain categorical markers in either mode.",
    )
    show_points = st.sidebar.toggle(
        "Show points on lines",
        value=True,
        key="show_points",
        disabled=style_label == "Bar",
        help="Display individual coordinate observations on line tracks.",
    )
    st.sidebar.button(
        "Reset display",
        width="stretch",
        on_click=_reset_display_state,
    )

    selection = Selection(
        stratum_type=stratum_type,
        stratum_value=stratum.value,
        block_number=block_number,
        window_size=size,
        renderer=Renderer.PLOTLY if interactive else Renderer.MATPLOTLIB,
        plot_style=PlotStyle.LINE if style_label == "Line" else PlotStyle.BAR,
        show_points=show_points,
    )
    previous_key = st.session_state.get("_scientific_selection_key")
    if previous_key is not None and tuple(previous_key) != selection.data_key:
        for key in ("region_query", "region_result", "filter_recipe", "filter_result"):
            st.session_state.pop(key, None)
    st.session_state._scientific_selection_key = selection.data_key
    return selection, stratum, block, navigation


def _header(
    content: AppContent,
    navigation: str,
    identity: DatabaseIdentity,
    stratum: StratumOption,
    block: BlockOption,
    size: int,
) -> None:
    page = content.page(navigation)
    st.markdown(page_kicker_html(content), unsafe_allow_html=True)
    icon_column, title_column = st.columns(
        PAGE_HEADER_COLUMN_WEIGHTS,
        gap=PAGE_HEADER_COLUMN_GAP,
        vertical_alignment="bottom",
    )
    with icon_column:
        # Streamlit 1.58 does not yet expose an ``alt`` argument on ``st.image``. The adjacent
        # semantic page heading provides the same page name; retain ``icon_alt`` in configuration
        # so it can be passed natively when the deployed Streamlit version supports it.
        # Pass bytes instead of a filesystem path so rendering is independent of the deployment
        # process's working directory and Streamlit's handling of local media paths.
        st.image(page.icon_path.read_bytes(), width=PAGE_HEADER_ICON_WIDTH)
    with title_column:
        st.markdown(page_title_html(content, navigation), unsafe_allow_html=True)
    st.markdown(page_lead_html(content, navigation), unsafe_allow_html=True)
    chips = (
        f"<span class='context-chip'>{escape(stratum.label)} · "
        f"{stratum.n_samples:,} samples</span>"
        f"<span class='context-chip'>Block {block.block_number:02d} · "
        f"AA {block.aa_start}–{block.aa_end}</span>"
        f"<span class='context-chip'>{size} AA windows</span>"
        f"<span class='context-chip'>DB {escape(identity.schema_version)}</span>"
    )
    ribbon = (
        "<div class='context-ribbon'>"
        f"<div class='context-label'>{escape(content.current_configuration_label)}</div>"
        f"{chips}</div>"
    )
    st.markdown(ribbon, unsafe_allow_html=True)


def _section_header(content: AppContent, key: str) -> None:
    """Render one JSON-backed section heading and escaped lead paragraph."""

    section = content.section(key)
    st.subheader(section.heading)
    st.markdown(section_intro_html(section), unsafe_allow_html=True)


def _formatted_metric(value: object, unit: str, decimals: int) -> str:
    if unit in {"fraction", "rate", "proportion"}:
        return format_fraction(value, decimals=max(1, min(decimals, 2)))
    if unit in {"count", "samples", "sites", "AA"}:
        return format_count(value)
    return format_number(value, decimals=min(decimals, 3))


def _metric_help(registry: pd.DataFrame, name: str) -> str:
    row = registry.loc[registry.name.eq(name)]
    if row.empty:
        return name
    return str(row.iloc[0].help_tip)


def _metric_card(
    container: object,
    row: pd.Series,
    registry: pd.DataFrame,
    name: str,
    label: str,
) -> None:
    definition = registry.loc[registry.name.eq(name)]
    unit = "index" if definition.empty else str(definition.iloc[0].unit)
    decimals = 3 if definition.empty else int(definition.iloc[0].decimals)
    container.metric(
        label,
        _formatted_metric(row.get(name), unit, decimals),
        help=_metric_help(registry, name),
    )


def _metric_guide(catalog: pd.DataFrame, level: str, content: AppContent) -> None:
    """Render the complete JSON-backed glossary for one analytical entity level."""

    guidance = metric_help_for_level(catalog, level)
    section = content.section("metric_guide")
    with st.expander(section.heading):
        st.markdown(section.intro)
        st.markdown(metric_guide_html(guidance), unsafe_allow_html=True)


def _block_plots(
    selection: Selection,
    summary: pd.Series,
    histogram: pd.DataFrame,
    exact: pd.DataFrame,
) -> None:
    from aa_metrics_app.charts.block import (
        matplotlib_callability_histogram,
        matplotlib_rank_distribution,
        plotly_callability_histogram,
        plotly_rank_distribution,
    )

    top_n = 0
    if not exact.empty:
        candidates = sorted({min(len(exact), value) for value in (10, 20, 50, len(exact)) if value})
        top_n = st.selectbox(
            "Haplotypes to show",
            candidates,
            index=min(1, len(candidates) - 1),
            key="block_rank_top_n",
            help="The complete exact rank distribution remains available in the table below.",
        )
    left, right = st.columns(2)
    if selection.renderer is Renderer.PLOTLY:
        _plotly_chart(
            left,
            plotly_callability_histogram(histogram, mean=summary.get("mean_sample_call_rate")),
        )
        if exact.empty:
            right.info("No exact block haplotypes are available for this selection.")
        else:
            _plotly_chart(
                right,
                plotly_rank_distribution(exact, top_n=top_n),
            )
    else:
        callability = matplotlib_callability_histogram(
            histogram, mean=summary.get("mean_sample_call_rate")
        )
        left.pyplot(callability, width="stretch")
        if exact.empty:
            right.info("No exact block haplotypes are available for this selection.")
        else:
            ranks = matplotlib_rank_distribution(exact, top_n=top_n)
            right.pyplot(ranks, width="stretch")


def _sequence_panel(
    block: BlockOption,
    summary: pd.Series,
    exact: pd.DataFrame,
    all_call: pd.DataFrame,
    content: AppContent,
) -> None:
    section = content.section("representative_block_sequences")
    st.markdown(f"#### {section.heading}")
    st.markdown(section.intro)
    rows = align_sequences(
        block.reference_sequence,
        (
            ("Consensus", summary.get("allelic_consensus_sequence")),
            ("Consensus (all)", summary.get("all_call_consensus_sequence")),
            (TOP_HAP_LABEL, summary.get("top_exact_sequence")),
            (TOP_HAP_ALL_LABEL, summary.get("modal_all_call_sequence")),
        ),
    )
    st.markdown(alignment_html(rows, start=block.aa_start), unsafe_allow_html=True)
    st.caption(
        "Orange underline: differs from 3D7 · Purple: mixed call · Grey: missing · Coral: "
        "unresolved. A bracketed mixed call occupies one AA coordinate."
    )
    exact_tab, all_tab = st.tabs(("Top complete haplotypes", "Top haplotypes (all)"))
    columns = (
        "sequence_rank",
        "sequence",
        "observation_count",
        "observation_denominator",
        "frequency",
        "cumulative_frequency",
        "distance_to_reference",
        "is_reference_sequence",
        "is_allelic_consensus_sequence",
    )
    with exact_tab:
        if exact.empty:
            st.info("No complete haplotype observations are available.")
        else:
            exact_columns = available_sequence_rank_columns(exact, columns)
            exact_view = exact.loc[:, exact_columns].head(10)
            st.dataframe(
                exact_view,
                hide_index=True,
                width="stretch",
                column_config={"observation_denominator": "Complete haplotype observations"},
            )
    with all_tab:
        if all_call.empty:
            st.info("No all-sample haplotype patterns are available.")
        else:
            all_call_columns = available_sequence_rank_columns(all_call, columns)
            all_view = all_call.loc[:, all_call_columns].head(10)
            st.dataframe(
                all_view,
                hide_index=True,
                width="stretch",
                column_config={"observation_denominator": "All selected samples"},
            )


def _mutation_explorer(
    path: str,
    selection: Selection,
    aa_start: int,
    aa_end: int,
    *,
    context_key: str,
    content: AppContent,
) -> None:
    """Render the shared geographic mutation comparison for one inclusive interval."""

    section = content.section("mutation_explorer")
    st.markdown(f"#### {section.heading}")
    st.markdown(section.intro)
    include_missing = st.toggle(
        "Include missing-call fractions for every group",
        value=False,
        key=f"mutation_missing_{context_key}",
        help=(
            "Add one missing-fraction column per stored group. This is the fraction of all "
            "samples with a literal '-' at the mutation locus."
        ),
    )
    wide = load_mutation_comparison(
        path,
        aa_start,
        aa_end,
        selection.stratum_type.value,
        selection.stratum_value,
        include_missing,
    )
    tidy = load_mutation_facts(path, aa_start, aa_end)
    if wide.empty:
        st.info(
            f"No canonical non-reference amino-acid mutation is catalogued for AA "
            f"{aa_start}–{aa_end}."
        )
        return

    selected = tidy.loc[
        tidy.stratum_type.eq(selection.stratum_type.value)
        & tidy.stratum_value.eq(selection.stratum_value)
    ]
    observed = selected.loc[selected.is_observed_in_stratum.eq(1)]
    cards = st.columns(4)
    cards[0].metric("Catalogued mutations", format_count(len(wide)))
    cards[1].metric("Observed in selected group", format_count(len(observed)))
    cards[2].metric(
        "Median observed frequency",
        format_fraction(observed.mutation_allele_frequency.median()) if not observed.empty else "-",
    )
    cards[3].metric("Selected interval", f"{aa_start}-{aa_end}")
    if (selected.allele_denominator_fractional.lt(20)).any():
        st.warning(
            "At least one selected-group mutation locus has fewer than 20 callable allele "
            "units. Compare frequency with the stored callable denominator and missingness.",
            icon="⚠️",
        )

    frequency_columns = [column for column in wide if column.startswith("frequency__")]
    missing_columns = [column for column in wide if column.startswith("missing__")]
    selected_label = next(iter(selected.display_label), selection.stratum_value)
    selected_prefix = f"{selection.stratum_type.value}/{selected_label}"
    column_labels: dict[str, str] = {
        "mutation": "Mutation",
        "aa_position": "AA position",
        "reference_residue": "3D7 residue",
        "alternate_residue": "Alternate residue",
        "block_number": "Block",
    }
    for column in frequency_columns:
        group = column.removeprefix("frequency__")
        column_labels[column] = (
            f"Selected frequency · {group}" if group == selected_prefix else f"Frequency · {group}"
        )
    for column in missing_columns:
        group = column.removeprefix("missing__")
        column_labels[column] = (
            f"Selected missing · {group}" if group == selected_prefix else f"Missing · {group}"
        )
    st.dataframe(
        wide,
        hide_index=True,
        width="stretch",
        height=min(620, 110 + 35 * len(wide)),
        column_config=column_labels,
    )
    context = {
        "stratum_type": selection.stratum_type.value,
        "stratum_value": selection.stratum_value,
        "block_number": selection.block_number,
        "window_size": selection.window_size,
        "aa_start": aa_start,
        "aa_end": aa_end,
        "include_missing_in_wide_table": include_missing,
    }
    bundle = build_mutation_bundle(
        wide,
        tidy,
        context=context,
        database_identity=load_identity(path),
    )
    downloads = st.columns(2)
    downloads[0].download_button(
        "Download mutation comparison (.zip)",
        bundle.content,
        file_name=bundle.filename,
        mime=bundle.mime_type,
        key=f"mutation_zip_{context_key}",
        width="stretch",
        help="Includes the wide comparison, tidy denominator-rich facts, query, and provenance.",
    )
    downloads[1].download_button(
        "Download wide table (.csv)",
        safe_csv(wide),
        file_name=bundle.filename.removesuffix(".zip") + ".csv",
        mime="text/csv",
        key=f"mutation_csv_{context_key}",
        width="stretch",
    )
    guidance = load_mutation_help(str(MUTATION_HELP_PATH))
    with st.expander("How mutation frequencies and missingness are calculated"):
        for key in ("canonical_id", "frequency", "missing", "warning"):
            st.markdown(guidance[key])


def _residue_logo(
    path: str,
    selection: Selection,
    aa_start: int,
    aa_end: int,
    *,
    context_key: str,
    content: AppContent,
    wrap: int = 50,
) -> None:
    """Render and export a sample-normalized residue-composition logo."""

    from aa_metrics_app.charts.logo import (
        matplotlib_logo,
        plotly_logo,
        prepare_logo_data,
    )

    section = content.section("residue_logo")
    st.markdown(f"#### {section.heading}")
    st.markdown(section.intro)
    spectrum = load_allele_spectrum(
        path,
        selection.stratum_type.value,
        selection.stratum_value,
        aa_start,
        aa_end,
    )
    logo = prepare_logo_data(spectrum, wrap=wrap)
    title = f"Residue composition · AA {aa_start}–{aa_end}"
    if selection.renderer is Renderer.PLOTLY:
        _plotly_chart(
            st,
            plotly_logo(logo, title=title),
            key=f"residue_logo_{context_key}",
        )
    else:
        st.pyplot(matplotlib_logo(logo, title=title), width="stretch")
    st.download_button(
        "Download plotted residue values (.csv)",
        safe_csv(logo.values),
        file_name=(
            f"msp1_{selection.stratum_type.value}-{selection.stratum_value}_"
            f"aa-{aa_start}-{aa_end}_residue-logo.csv"
        ),
        mime="text/csv",
        key=f"residue_logo_download_{context_key}",
        help="Exact plotted heights, empty mass, coordinates, and consensus flags.",
    )
    guidance = load_logo_help(str(LOGO_HELP_PATH))
    with st.expander("How to read the residue-composition logo"):
        for key in ("mixed", "uncalled", "consensus", "interpretation"):
            st.markdown(guidance[key])


def _resolve_structure_asset(database_path: str, relative_path: str) -> Path:
    """Resolve a database-pinned relative asset without network access."""

    artifact_root = _repository(database_path, validate=False).artifact_root
    candidates = (
        artifact_root / relative_path,
        PROJECT_ROOT / relative_path,
    )
    return next((candidate for candidate in candidates if candidate.is_file()), candidates[0])


def _structure_context(
    path: str,
    selection: Selection,
    aa_start: int,
    aa_end: int,
    reference_sequence: str,
    top_complete_sequence: str | None,
    *,
    context_key: str,
    content: AppContent,
) -> None:
    """Render auditable local reference-structure context for one interval."""

    from aa_metrics_app.structures import (
        map_mismatches,
        plotly_structure_context,
        top_complete_mismatches,
    )

    section = content.section("structure_context")
    st.markdown(f"#### {section.heading}")
    st.markdown(section.intro)
    candidates = load_structure_candidates(path, aa_start, aa_end)
    if candidates.empty:
        st.info("No reviewed local structure is available for this interval.")
        return
    adequate = candidates.loc[candidates.adequate_interval_coverage].copy()
    with st.expander("Reviewed structure candidates and interval coverage"):
        st.dataframe(
            candidates.loc[
                :,
                [
                    "accession",
                    "chain_id",
                    "source_type",
                    "experimental_method",
                    "resolution_angstrom",
                    "interval_coverage",
                    "mapping_identity",
                    "adequate_interval_coverage",
                ],
            ],
            hide_index=True,
            width="stretch",
        )
    if adequate.empty:
        st.info(
            "No reviewed structure maps at least 80% of this interval. Mutation and sequence "
            "tables remain available without an approximate 3D placement."
        )
        return
    labels = {
        int(index): (
            f"{row.accession} · chain {row.chain_id} · {row.source_type} · "
            f"{float(row.interval_coverage):.1%} coverage"
        )
        for index, row in adequate.iterrows()
    }
    selected_index = st.selectbox(
        "Structure model",
        list(labels),
        format_func=labels.get,
        key=f"structure_model_{context_key}",
        help=("Adequate experimental models are listed first, then the pinned AlphaFold fallback."),
    )
    selected = adequate.loc[selected_index]
    residue_map = load_structure_map(
        path,
        int(selected.structure_chain_id),
        aa_start,
        aa_end,
    )
    mismatches = top_complete_mismatches(
        reference_sequence,
        top_complete_sequence,
        aa_start=aa_start,
    )
    mapped, unmapped = map_mismatches(mismatches, residue_map)
    cards = st.columns(4)
    cards[0].metric("Model", str(selected.accession))
    cards[1].metric("Source", str(selected.source_type).title())
    cards[2].metric("Interval coverage", format_fraction(selected.interval_coverage))
    cards[3].metric("Mapped top-sequence differences", format_count(len(mapped)))
    if top_complete_sequence is None or pd.isna(top_complete_sequence):
        st.info("No top complete sequence is available, so no substitutions are highlighted.")
    elif not unmapped.empty:
        st.warning(
            f"{len(unmapped)} top-complete difference(s) are not modeled in this chain and are "
            "listed below rather than placed approximately.",
            icon="⚠️",
        )
    asset = _resolve_structure_asset(path, str(selected.asset_relative_path))
    try:
        backbone = load_structure_backbone(
            str(asset),
            str(selected.asset_sha256),
            str(selected.chain_id),
            residue_map,
        )
        if selection.renderer is Renderer.PLOTLY:
            _plotly_chart(
                st,
                plotly_structure_context(
                    backbone,
                    aa_start=aa_start,
                    aa_end=aa_end,
                    mapped_mismatches=mapped,
                    title=(f"{selected.accession} chain {selected.chain_id} · reference structure"),
                ),
                key=f"structure_context_{context_key}",
            )
        else:
            st.info(
                "Interactive plots are disabled. Use the mapping tables below or enable Plotly "
                "in the sidebar to view the rotatable 3D backbone."
            )
    except (FileNotFoundError, ValueError) as error:
        st.warning(
            "The local 3D coordinate view could not be rendered safely. The reviewed mapping "
            f"tables remain available. ({type(error).__name__})",
            icon="⚠️",
        )
    tabs = st.tabs(("Mapped interval", "Top-sequence differences", "Unmapped differences"))
    with tabs[0]:
        st.dataframe(residue_map, hide_index=True, width="stretch")
    with tabs[1]:
        st.dataframe(mapped, hide_index=True, width="stretch")
    with tabs[2]:
        st.dataframe(unmapped, hide_index=True, width="stretch")
    guidance = load_structure_help(str(STRUCTURE_HELP_PATH))
    with st.expander("How the structure model and highlights are chosen"):
        for key in ("selection", "interpretation", "confidence", "fallback"):
            st.markdown(guidance[key])


def _block_overview(
    path: str,
    selection: Selection,
    stratum: StratumOption,
    block: BlockOption,
    content: AppContent,
) -> None:
    _section_header(content, "block_overview")
    summary = load_block_summary(
        path,
        selection.stratum_type.value,
        selection.stratum_value,
        selection.block_number,
    ).iloc[0]
    histogram = load_block_histogram(
        path,
        selection.stratum_type.value,
        selection.stratum_value,
        selection.block_number,
    )
    exact = load_block_sequences(
        path,
        selection.stratum_type.value,
        selection.stratum_value,
        selection.block_number,
        "exact",
    )
    all_call = load_block_sequences(
        path,
        selection.stratum_type.value,
        selection.stratum_value,
        selection.block_number,
        "all_call",
    )
    registry = load_metric_registry(path, "block")

    metadata = st.columns(4)
    metadata[0].metric("Block", f"{block.block_number:02d}")
    metadata[1].metric("Coordinates", f"{block.aa_start}–{block.aa_end}")
    metadata[2].metric("Length", f"{block.length} AA")
    metadata[3].metric("Selected samples", format_count(stratum.n_samples))

    headline_metrics = (
        ("mean_sample_call_rate", "Mean call rate"),
        ("mean_sample_missing_fraction", "Mean missing"),
        ("mean_sample_mixed_fraction_callable", "Mean mixed / callable"),
        ("n_unique_exact_haplotypes", "Unique complete haplotypes"),
        ("effective_haplotypes_simpson", "Effective haplotypes"),
        ("haplotype_diversity_hd", "Haplotype diversity"),
    )
    for offset in (0, 3):
        headline = st.columns(3)
        for container, (name, label) in zip(
            headline, headline_metrics[offset : offset + 3], strict=True
        ):
            _metric_card(container, summary, registry, name, label)

    evidence = st.columns(3)
    _metric_card(evidence[0], summary, registry, "n_samples", "All sample observations")
    _metric_card(evidence[1], summary, registry, "n_fully_callable_sequences", "Fully callable")
    _metric_card(
        evidence[2],
        summary,
        registry,
        "n_exact_haplotype_observations",
        "Complete haplotype observations",
    )
    frequencies = st.columns(5)
    for container, name, label in zip(
        frequencies,
        (
            "reference_haplotype_frequency_exact",
            "allelic_consensus_exact_frequency",
            "all_call_consensus_sequence_frequency",
            "top_exact_sequence_frequency",
            "modal_all_call_sequence_frequency",
        ),
        (
            "3D7 haplotype frequency",
            "Consensus haplotype frequency",
            "All-call consensus frequency",
            "Top complete haplotype frequency",
            "Top haplotype frequency (all)",
        ),
        strict=True,
    ):
        _metric_card(container, summary, registry, name, label)
    st.caption(
        "3D7, canonical consensus, and top complete frequencies use complete observations. "
        "All-call consensus and top haplotype (all) use every selected sample. Mixed calls "
        "contribute allele evidence to the all-call consensus; missing can be its emitted state."
    )
    if int(summary.get("n_exact_haplotype_observations", 0)) < 20:
        st.warning(
            "The number of complete haplotypes is small; frequency and diversity estimates are "
            "descriptive and may be unstable.",
            icon="⚠️",
        )

    _block_plots(selection, summary, histogram, exact)
    _sequence_panel(block, summary, exact, all_call, content)
    _residue_logo(
        path,
        selection,
        block.aa_start,
        block.aa_end,
        context_key="block",
        content=content,
    )
    _structure_context(
        path,
        selection,
        block.aa_start,
        block.aa_end,
        block.reference_sequence,
        summary.get("top_exact_sequence"),
        context_key="block",
        content=content,
    )
    _mutation_explorer(
        path,
        selection,
        block.aa_start,
        block.aa_end,
        context_key="block",
        content=content,
    )

    with st.expander("All block metric values"):
        details = registry.copy()
        details.insert(2, "value", details.name.map(summary.to_dict()))
        st.dataframe(
            details.loc[
                :,
                [
                    "group",
                    "display_name",
                    "value",
                    "unit",
                ],
            ],
            hide_index=True,
            width="stretch",
        )
        st.markdown(
            "**Complete vs all samples.** Complete distributions exclude missing, unresolved, "
            "mixed, and non-canonical calls. Top haplotype (all) uses every sample and is not an "
            "allelic consensus. The sitewise allelic consensus may not itself be observed. "
            "Singleton residues use sample presence, including presence inside mixed calls."
        )


LOCUS_DEFAULT_TRACKS = (
    "call_rate",
    "consensus_frequency",
    "reference_residue_frequency",
    "is_polymorphic_any",
    "is_polymorphic_ge_0_01",
    "n_observed_residues",
    "second_residue_frequency",
    "n_residues_frequency_ge_0_05",
    "expected_heterozygosity",
    "shannon_entropy",
)

LOCUS_BINARY_TRACKS = {
    "contains_reference_residue",
    "consensus_tied",
    "is_polymorphic_any",
    "is_polymorphic_ge_0_01",
    "is_polymorphic_ge_0_05",
    "has_singleton_residue",
}

def _locus_y_range(name: str, unit: str) -> tuple[float, float] | None:
    if unit == "fraction" or name in {
        "expected_heterozygosity",
        "allele_frequency_gini",
        "pielou_evenness",
        "normalized_entropy_observed",
        "normalized_entropy_20",
    }:
        return (0.0, 1.0)
    if name in LOCUS_BINARY_TRACKS:
        return (-0.1, 1.1)
    return None


def _render_locus_track(
    selection: Selection,
    frame: pd.DataFrame,
    definition: pd.Series,
) -> None:
    metric = str(definition["name"])
    plot_frame, hover = prepare_locus_hover_data(frame)
    spec = TrackSpec(
        title=str(definition["display_name"]),
        x_label="MSP1 amino-acid coordinate",
        y_label=str(definition["unit"]),
        x="aa_position",
        y=metric,
        style=selection.plot_style,
        color=PRIMARY_TOMATO,
        y_range=_locus_y_range(metric, str(definition["unit"])),
        height=115 if metric in LOCUS_BINARY_TRACKS else 165,
        binary=metric in LOCUS_BINARY_TRACKS,
        hover_columns=hover,
        description=str(definition["description"]),
        unit=str(definition["unit"]),
    )
    if selection.renderer is Renderer.PLOTLY:
        from aa_metrics_app.charts.plotly import render_track

        rendered = render_track(plot_frame, spec, show_points=selection.show_points)
        _plotly_chart(
            st,
            rendered.figure,
            key=f"locus_plotly_{metric}",
        )
    else:
        from aa_metrics_app.charts.matplotlib import render_track

        rendered = render_track(plot_frame, spec, show_points=selection.show_points)
        st.pyplot(rendered.figure, width="stretch")


def _locus_tracks(path: str, selection: Selection, block: BlockOption, content: AppContent) -> None:
    _section_header(content, "locus_tracks")
    frame = load_locus_track(
        path,
        selection.stratum_type.value,
        selection.stratum_value,
        selection.block_number,
    )
    registry = load_metric_registry(path, "locus")
    registry = registry.loc[registry.name.isin(frame.columns)].copy()
    label_by_name = {
        row["name"]: f"{str(row['group']).title()} · {row['display_name']}"
        for row in registry.to_dict("records")
    }
    defaults = [name for name in LOCUS_DEFAULT_TRACKS if name in label_by_name]
    show_all = st.toggle(
        "Show all locus tracks",
        value=False,
        key="show_all_locus_tracks",
        help=(
            "By default, a selection of sensible metrics are displayed. Enable this to render "
            "a track for every pre-computed locus metric."
        ),
    )
    if show_all:
        selected = registry.name.tolist()
        st.caption(f"Showing all {len(selected)} stored locus metrics.")
    else:
        selected = st.multiselect(
            "Tracks to show",
            registry.name.tolist(),
            default=defaults,
            format_func=label_by_name.get,
            key="locus_tracks_to_show",
            help="Choose compact tracks by scientific group. Coordinates are shown on every plot.",
        )
    if not selected:
        st.info("Choose at least one locus metric to display.")
    for metric in selected:
        definition = registry.loc[registry.name.eq(metric)].iloc[0]
        if frame[metric].notna().any():
            _render_locus_track(selection, frame, definition)
        else:
            st.info(f"{definition.display_name} is unavailable for this selection.")

    detail_columns = [
        "aa_position",
        "block_offset",
        "reference_residue",
        "consensus_residue",
        "consensus_residues",
        "second_residue",
        "n_samples",
        "allele_denominator_fractional",
        *selected,
    ]
    detail_columns = list(dict.fromkeys(column for column in detail_columns if column in frame))
    with st.expander("Locus values and download"):
        st.dataframe(frame.loc[:, detail_columns], hide_index=True, width="stretch")
        st.download_button(
            "Download selected locus values (.csv)",
            frame.loc[:, detail_columns].to_csv(index=False).encode("utf-8"),
            file_name=(
                f"msp1_{selection.stratum_type.value}-{selection.stratum_value}_"
                f"block-{block.block_number:02d}_loci.csv"
            ),
            mime="text/csv",
            key="download_locus_values",
        )
        st.markdown(
            "A mixed call such as `[A/C]` contributes 0.5 to A and 0.5 to C. "
            "`contains_reference_residue` can be true even when the reference is not the "
            "consensus. Tied consensus residues are retained in `consensus_residues`; entropy, "
            "Gini, and Pielou may be undefined when sample support or residue richness is "
            "insufficient."
        )


WINDOW_DEFAULT_TRACKS = (
    "mean_sample_call_rate",
    "min_locus_call_rate",
    "mean_sample_non_callable_fraction",
    "top_exact_peptide_frequency",
    "fraction_sample_peptides_exact",
    "peptide_diversity_hd",
    "effective_peptides_simpson",
    "fraction_sites_polymorphic_ge_0_01",
)


def _window_y_range(name: str, unit: str) -> tuple[float, float] | None:
    if unit == "fraction" or name in {
        "peptide_diversity_hd",
        "peptide_frequency_gini",
        "pielou_evenness",
    }:
        return (0.0, 1.0)
    return None


def _render_window_track(
    selection: Selection,
    frame: pd.DataFrame,
    definition: pd.Series,
) -> None:
    metric = str(definition["name"])
    hover = tuple(
        column
        for column in (
            "aa_end",
            "block_offset_start",
            "block_offset_end",
            "n_samples",
            "n_fully_callable_peptide_observations",
            "n_exact_peptide_observations",
            "reference_sequence",
            "top_exact_peptide_sequence",
            "allelic_consensus_sequence",
            "all_call_consensus_sequence",
            "modal_all_call_sequence",
        )
        if column in frame.columns
    )
    spec = TrackSpec(
        title=str(definition["display_name"]),
        x_label=f"Window start coordinate ({selection.window_size} AA)",
        y_label=str(definition["unit"]),
        x="aa_start",
        y=metric,
        style=selection.plot_style,
        color=PRIMARY_TOMATO,
        y_range=_window_y_range(metric, str(definition["unit"])),
        height=220,
        hover_columns=hover,
        description=str(definition["description"]),
        unit=str(definition["unit"]),
    )
    if selection.renderer is Renderer.PLOTLY:
        from aa_metrics_app.charts.plotly import render_track

        rendered = render_track(frame, spec, show_points=selection.show_points)
        _plotly_chart(
            st,
            rendered.figure,
            key=f"window_plotly_{metric}",
        )
    else:
        from aa_metrics_app.charts.matplotlib import render_track

        rendered = render_track(frame, spec, show_points=selection.show_points)
        st.pyplot(rendered.figure, width="stretch")


def _window_inspector(frame: pd.DataFrame, block: BlockOption, content: AppContent) -> pd.Series:
    section = content.section("peptide_sequence_inspector")
    st.markdown(f"#### {section.heading}")
    st.markdown(section.intro)
    display_key = st.selectbox(
        "Window to inspect",
        frame.display_key.tolist(),
        key="window_inspector_key",
        help=(
            "Choose one stored window to compare reference, consensus, top complete, and "
            "all-sample sequences."
        ),
    )
    row = frame.loc[frame.display_key.eq(display_key)].iloc[0]
    rows = align_sequences(
        str(row.reference_sequence),
        (
            ("Consensus", row.allelic_consensus_sequence),
            ("Consensus (all)", row.all_call_consensus_sequence),
            (TOP_COMPLETE_PEP_LABEL, row.top_exact_peptide_sequence),
            (TOP_PEPTIDE_ALL_LABEL, row.modal_all_call_sequence),
        ),
        annotations={
            "3D7": format_fraction(row.reference_peptide_frequency_exact),
            "Consensus": format_fraction(row.allelic_consensus_exact_frequency),
            "Consensus (all)": format_fraction(row.all_call_consensus_sequence_frequency),
            TOP_COMPLETE_PEP_LABEL: format_fraction(row.top_exact_peptide_frequency),
            TOP_PEPTIDE_ALL_LABEL: format_fraction(row.modal_all_call_sequence_frequency),
        },
    )
    st.markdown(alignment_html(rows, start=int(row.aa_start), wrap=50), unsafe_allow_html=True)
    summary = st.columns(4)
    summary[0].metric("Coordinates", f"{int(row.aa_start)}–{int(row.aa_end)}")
    summary[1].metric("All samples", format_count(row.n_samples))
    summary[2].metric("Fully callable", format_count(row.n_fully_callable_peptide_observations))
    summary[3].metric(
        "Complete peptide observations", format_count(row.n_exact_peptide_observations)
    )
    if int(row.n_exact_peptide_observations) < 20:
        st.warning(
            "This window has fewer than 20 complete peptide observations; complete-sequence "
            "frequencies and diversity may be unstable.",
            icon="⚠️",
        )
    return row


def _window_tracks(
    path: str, selection: Selection, block: BlockOption, content: AppContent
) -> None:
    _section_header(content, "window_tracks")
    frame = load_window_track(
        path,
        selection.stratum_type.value,
        selection.stratum_value,
        selection.block_number,
        selection.window_size,
    )
    registry = load_metric_registry(path, "window")
    registry = registry.loc[registry.name.isin(frame.columns)].copy()
    label_by_name = {
        row["name"]: f"{str(row['group']).title()} · {row['display_name']}"
        for row in registry.to_dict("records")
    }
    defaults = [name for name in WINDOW_DEFAULT_TRACKS if name in label_by_name]
    show_all = st.toggle(
        "Show all window tracks",
        value=False,
        key="show_all_window_tracks",
        help=(
            "Enable every pre-computed numerical metric, instead of the default selection of "
            "sensible metrics."
        ),
    )
    if show_all:
        selected = registry.name.tolist()
        st.caption(f"Showing all {len(selected)} stored window metrics.")
    else:
        selected = st.multiselect(
            "Tracks to show",
            registry.name.tolist(),
            default=defaults,
            format_func=label_by_name.get,
            key="window_tracks_to_show",
            help="All points use aa_start; hover includes aa_end, sequences, and supporting calls.",
        )
    if not selected:
        st.info("Choose at least one window metric to display.")
    for metric in selected:
        definition = registry.loc[registry.name.eq(metric)].iloc[0]
        if frame[metric].notna().any():
            _render_window_track(selection, frame, definition)
        else:
            st.info(f"{definition.display_name} is unavailable for this selection.")

    inspected = _window_inspector(frame, block, content)
    _residue_logo(
        path,
        selection,
        int(inspected.aa_start),
        int(inspected.aa_end),
        context_key="window",
        content=content,
        wrap=selection.window_size,
    )
    _structure_context(
        path,
        selection,
        int(inspected.aa_start),
        int(inspected.aa_end),
        str(inspected.reference_sequence),
        inspected.top_exact_peptide_sequence,
        context_key="window",
        content=content,
    )
    _mutation_explorer(
        path,
        selection,
        int(inspected.aa_start),
        int(inspected.aa_end),
        context_key="window",
        content=content,
    )
    pinned = [
        "display_key",
        "block_number",
        "block_name",
        "window_size",
        "aa_start",
        "aa_end",
        "reference_sequence",
        "allelic_consensus_sequence",
        "all_call_consensus_sequence",
        "top_exact_peptide_sequence",
        "modal_all_call_sequence",
        "n_samples",
        "n_fully_callable_peptide_observations",
        "n_exact_peptide_observations",
        *selected,
    ]
    pinned = list(dict.fromkeys(column for column in pinned if column in frame))
    with st.expander("Window values and download"):
        window_view = frame.loc[:, pinned]
        st.dataframe(
            window_view,
            hide_index=True,
            width="stretch",
            column_config=_column_config(window_view.columns, registry),
        )
        st.download_button(
            "Download selected window values (.csv)",
            frame.loc[:, pinned].to_csv(index=False).encode("utf-8"),
            file_name=(
                f"msp1_{selection.stratum_type.value}-{selection.stratum_value}_"
                f"block-{block.block_number:02d}_{selection.window_size}aa_windows.csv"
            ),
            mime="text/csv",
            key="download_window_values",
        )
        st.markdown(
            "Adjacent windows overlap and therefore are not independent. Complete peptide "
            "metrics exclude missing, unresolved, mixed, and non-canonical calls. All-sample "
            "patterns keep "
            "those normalized states. Raw pairwise distance is AA mismatches; normalized distance "
            "divides by peptide length. Missingness is a quality concern, not direct evidence of "
            "an indel."
        )


def _region_context_chart(
    selection: Selection,
    rows: pd.DataFrame,
    block: BlockOption,
    query: RegionQuery,
) -> None:
    if selection.renderer is Renderer.PLOTLY:
        from aa_metrics_app.charts.region import plotly_region_context

        _plotly_chart(
            st,
            plotly_region_context(
                rows,
                block_start=block.aa_start,
                block_end=block.aa_end,
                query_start=query.start,
                query_end=query.end,
            ),
            key="region_context_plotly",
        )
    else:
        from aa_metrics_app.charts.region import matplotlib_region_context

        figure = matplotlib_region_context(
            rows,
            block_start=block.aa_start,
            block_end=block.aa_end,
            query_start=query.start,
            query_end=query.end,
        )
        st.pyplot(figure, width="stretch")


def _region_sequence_inspector(rows: pd.DataFrame, content: AppContent) -> pd.Series:
    section = content.section("matching_peptide_inspector")
    st.markdown(f"#### {section.heading}")
    st.markdown(section.intro)
    display_key = st.selectbox(
        "Matching peptide",
        rows.display_key.tolist(),
        key="region_inspector_key",
        help="Inspect the four sequence meanings for one matching peptide.",
    )
    row = rows.loc[rows.display_key.eq(display_key)].iloc[0]
    aligned = align_sequences(
        str(row.reference_sequence),
        (
            ("Consensus", row.allelic_consensus_sequence),
            ("Consensus (all)", row.all_call_consensus_sequence),
            (TOP_COMPLETE_PEP_LABEL, row.top_exact_peptide_sequence),
            (TOP_PEPTIDE_ALL_LABEL, row.modal_all_call_sequence),
        ),
        annotations={
            "3D7": format_fraction(row.reference_peptide_frequency_exact),
            "Consensus": format_fraction(row.allelic_consensus_exact_frequency),
            "Consensus (all)": format_fraction(row.all_call_consensus_sequence_frequency),
            TOP_COMPLETE_PEP_LABEL: format_fraction(row.top_exact_peptide_frequency),
            TOP_PEPTIDE_ALL_LABEL: format_fraction(row.modal_all_call_sequence_frequency),
        },
    )
    st.markdown(alignment_html(aligned, start=int(row.aa_start)), unsafe_allow_html=True)
    return row


def _region_explorer(
    path: str, selection: Selection, block: BlockOption, content: AppContent
) -> None:
    _section_header(content, "region_explorer")
    frame = load_window_track(
        path,
        selection.stratum_type.value,
        selection.stratum_value,
        selection.block_number,
        selection.window_size,
    )
    with st.form("region_query_form", clear_on_submit=False):
        left, right = st.columns(2)
        start = int(
            left.number_input(
                "Range start",
                min_value=block.aa_start,
                max_value=block.aa_end,
                value=block.aa_start,
                step=1,
                key=f"region_start_{block.block_number}",
                help="One-based inclusive MSP1 amino-acid coordinate.",
            )
        )
        end = int(
            right.number_input(
                "Range end",
                min_value=block.aa_start,
                max_value=block.aa_end,
                value=block.aa_start,
                step=1,
                key=f"region_end_{block.block_number}",
                help="Inclusive range end; it must be greater than or equal to the start.",
            )
        )
        match_label = st.radio(
            "Peptide relationship",
            ("Overlap any part", "Contain the entire range"),
            horizontal=True,
            help=(
                "Overlap includes any shared coordinate. Containment requires the peptide to "
                "span both submitted boundaries."
            ),
        )
        submit_col, clear_col = st.columns(2)
        submitted = submit_col.form_submit_button(
            "Show overlapping peptides", type="primary", width="stretch"
        )
        cleared = clear_col.form_submit_button("Clear region result", width="stretch")
    if cleared:
        for key in ("region_query", "region_result", "region_warnings"):
            st.session_state.pop(key, None)
        st.rerun()
    if submitted:
        try:
            query = RegionQuery(
                start=start,
                end=end,
                mode=(
                    RegionMode.OVERLAP if match_label == "Overlap any part" else RegionMode.CONTAIN
                ),
            )
            result = select_region(
                frame,
                query,
                block_start=block.aa_start,
                block_end=block.aa_end,
                window_size=selection.window_size,
            )
            st.session_state.region_query = query
            st.session_state.region_result = result.rows
            st.session_state.region_warnings = result.warnings
        except ValueError as error:
            st.error(str(error), icon="🚫")

    query = st.session_state.get("region_query")
    rows = st.session_state.get("region_result")
    warnings = st.session_state.get("region_warnings", ())
    if query is None or rows is None:
        st.markdown(
            "<div class='status-note'>Set an inclusive range and press <strong>Show "
            "overlapping peptides</strong>. Region results are independent from threshold "
            "filters.</div>",
            unsafe_allow_html=True,
        )
        return
    for warning in warnings:
        st.warning(warning, icon="⚠️")
    if rows.empty:
        st.info("The submitted interval has no matching peptides for this relationship.")
        return

    summary = st.columns(5)
    summary[0].metric("Submitted interval", f"{query.start}–{query.end}")
    summary[1].metric("Matching peptides", format_count(len(rows)))
    summary[2].metric("Coordinate span", f"{int(rows.aa_start.min())}–{int(rows.aa_end.max())}")
    summary[3].metric(
        "Median call rate", format_fraction(rows.mean_sample_call_rate.median(), decimals=1)
    )
    summary[4].metric(
        "Median top frequency",
        format_fraction(rows.top_exact_peptide_frequency.median(), decimals=1),
    )
    _region_context_chart(selection, rows, block, query)
    pinned = [
        "display_key",
        "block_number",
        "window_size",
        "aa_start",
        "aa_end",
        "overlap_length",
        "query_relation",
        "reference_sequence",
        "allelic_consensus_sequence",
        "all_call_consensus_sequence",
        "top_exact_peptide_sequence",
        "modal_all_call_sequence",
        "n_samples",
        "n_exact_peptide_observations",
        "mean_sample_call_rate",
        "top_exact_peptide_frequency",
        "peptide_diversity_hd",
    ]
    remaining = [column for column in rows.columns if column not in pinned]
    ordered = rows.loc[:, [column for column in (*pinned, *remaining) if column in rows]]
    region_registry = load_metric_registry(path, "window")
    st.dataframe(
        ordered,
        hide_index=True,
        width="stretch",
        column_config=_column_config(ordered.columns, region_registry),
    )
    query_payload = {
        "start": query.start,
        "end": query.end,
        "mode": query.mode.value,
        "stratum_type": selection.stratum_type.value,
        "stratum_value": selection.stratum_value,
        "block_number": selection.block_number,
        "window_size": selection.window_size,
        "matching_peptides": len(ordered),
        "display_keys": ordered.display_key.tolist(),
    }
    region_bundle = build_region_bundle(
        ordered,
        region_registry,
        query=query_payload,
        database_identity=load_identity(path),
        filename_context=query_payload,
    )
    downloads = st.columns(2)
    downloads[0].download_button(
        "Download region bundle (.zip)",
        region_bundle.content,
        file_name=region_bundle.filename,
        mime=region_bundle.mime_type,
        key="download_region_zip",
        width="stretch",
        help=(
            "Canonical bundle with data, query, metric definitions, database identity, and README."
        ),
    )
    downloads[1].download_button(
        "Download overlapping peptides (.csv)",
        safe_csv(ordered),
        file_name=(
            f"msp1_block-{block.block_number:02d}_{selection.window_size}aa_"
            f"region-{query.start}-{query.end}.csv"
        ),
        mime="text/csv",
        key="download_region_csv",
        width="stretch",
        help="Convenience table only; use the ZIP for reproducibility metadata.",
    )
    inspected = _region_sequence_inspector(rows, content)
    _residue_logo(
        path,
        selection,
        int(inspected.aa_start),
        int(inspected.aa_end),
        context_key="region",
        content=content,
        wrap=selection.window_size,
    )
    _structure_context(
        path,
        selection,
        int(inspected.aa_start),
        int(inspected.aa_end),
        str(inspected.reference_sequence),
        inspected.top_exact_peptide_sequence,
        context_key="region",
        content=content,
    )
    _mutation_explorer(
        path,
        selection,
        int(inspected.aa_start),
        int(inspected.aa_end),
        context_key="region",
        content=content,
    )
    with st.expander("How region matching works"):
        st.markdown(
            "Coordinates are one-based and inclusive. **Overlap** uses `peptide_start ≤ "
            "query_end` and `peptide_end ≥ query_start`. **Containment** uses `peptide_start ≤ "
            "query_start` and `peptide_end ≥ query_end`. A boundary-touching peptide therefore "
            "matches overlap. Every result has the peptide length selected in the sidebar."
        )


CONSERVATION_TARGET_LABELS = {
    ConservationTarget.TOP_EXACT: "Top complete peptide",
    ConservationTarget.REFERENCE: "3D7/reference peptide",
    ConservationTarget.ALLELIC_CONSENSUS: "Allelic consensus peptide",
}

FILTER_GROUPS = (
    (
        "Quality and evidence",
        (
            "mean_sample_call_rate",
            "min_locus_call_rate",
            "fraction_sample_peptides_exact",
            "mean_sample_missing_fraction",
            "mean_sample_mixed_fraction_callable",
            "max_locus_mixed_fraction_callable",
            "fraction_samples_missing_ge_1",
        ),
    ),
    (
        "Conservation and diversity",
        (
            "$conservation_frequency",
            "peptide_diversity_hd",
            "effective_peptides_simpson",
            "mean_pairwise_aa_distance",
            "fraction_sites_polymorphic_ge_0_01",
        ),
    ),
)


def _constraint_expression(constraints: tuple[MetricConstraint, ...]) -> str:
    parts = []
    for item in constraints:
        if item.minimum is not None:
            parts.append(f"{item.metric} ≥ {item.minimum:g}")
        if item.maximum is not None:
            parts.append(f"{item.metric} ≤ {item.maximum:g}")
    return " AND ".join(parts) if parts else "No active constraints"


def _constraints_match(
    submitted: tuple[MetricConstraint, ...], expected: tuple[MetricConstraint, ...]
) -> bool:
    """Compare widget-produced floats to preset values at display precision."""

    if len(submitted) != len(expected):
        return False
    submitted_by_metric = {item.metric: item for item in submitted}
    expected_by_metric = {item.metric: item for item in expected}
    if submitted_by_metric.keys() != expected_by_metric.keys():
        return False
    for metric, left in submitted_by_metric.items():
        right = expected_by_metric[metric]
        for left_value, right_value in (
            (left.minimum, right.minimum),
            (left.maximum, right.maximum),
        ):
            if (left_value is None) != (right_value is None):
                return False
            if left_value is not None and not math.isclose(
                left_value, right_value, rel_tol=0.0, abs_tol=5e-5
            ):
                return False
    return True


def _filter_number_input_options(
    definition: pd.Series,
) -> dict[str, float | str]:
    """Return bounded Streamlit number-input options for one metric."""

    options: dict[str, float | str] = {
        "step": float(definition.input_step),
        "format": "%.0f" if definition.data_type == "INTEGER" else "%.4f",
    }
    if pd.notna(definition.valid_minimum):
        options["min_value"] = float(definition.valid_minimum)
    maximum = definition.valid_maximum
    if pd.notna(maximum):
        options["max_value"] = float(maximum)
    return options


def _filter_explorer(
    path: str, selection: Selection, block: BlockOption, content: AppContent
) -> None:
    _section_header(content, "filter_explorer")
    frame = load_window_track(
        path,
        selection.stratum_type.value,
        selection.stratum_value,
        selection.block_number,
        selection.window_size,
    )
    registry = load_metric_registry(path, "window")
    filterable = registry.loc[registry.filterable & registry.name.isin(frame.columns)].copy()
    metric_definitions = tuple(
        item for item in load_metric_definitions(path, "window") if item.filterable
    )
    catalog = load_filter_presets(path)

    preset_id = st.selectbox(
        "Starting preset",
        [item.identifier for item in catalog.presets],
        index=1,
        format_func=lambda value: catalog.get(value).label,
        key="filter_preset_id",
        help="Presets are versioned preliminary thresholds. You can edit or deactivate any rule.",
    )
    preset = catalog.get(preset_id)
    st.caption(f"{preset.description} {catalog.scientific_caution}")
    target = st.radio(
        "Conservation target",
        tuple(ConservationTarget),
        format_func=CONSERVATION_TARGET_LABELS.get,
        horizontal=True,
        key="filter_conservation_target",
        help=(
            "Choose which complete-peptide frequency must meet the conservation threshold: the "
            "most common observed peptide, the exact 3D7 peptide, or the site-wise consensus. "
            "This choice does not change the other filters."
        ),
    )
    target_help = content.section("conservation_target")
    st.caption(
        "The selected target replaces the preset's generic conservation-frequency rule; all "
        "other active constraints remain unchanged."
    )
    with st.expander(target_help.heading):
        st.markdown(target_help.intro)
    resolved = preset.recipe(target)
    resolved_by_metric = {item.metric: item for item in resolved.constraints}
    core_names = {item.metric for item in resolved.constraints}
    advanced_options = [name for name in filterable.name if name not in core_names]
    advanced = st.multiselect(
        "Add advanced metric filters",
        advanced_options,
        format_func=lambda name: filterable.loc[filterable.name.eq(name), "display_name"].iloc[0],
        key="filter_advanced_metrics",
        help=(
            "Add only the extra database metrics you need; this avoids a wall of inactive sliders."
        ),
    )

    submitted_constraints: list[MetricConstraint] = []
    submitted = False
    cleared = False
    with st.form("peptide_filter_form", clear_on_submit=False):
        for title, raw_names in FILTER_GROUPS:
            with st.expander(title, expanded=True):
                for raw_name in raw_names:
                    metric = target.value if raw_name == "$conservation_frequency" else raw_name
                    item = resolved_by_metric.get(metric)
                    if item is None:
                        continue
                    definition = filterable.loc[filterable.name.eq(item.metric)].iloc[0]
                    direction = "minimum" if item.minimum is not None else "maximum"
                    default = item.minimum if item.minimum is not None else item.maximum
                    row = st.columns((0.20, 0.48, 0.32))
                    active = row[0].checkbox(
                        "Active",
                        value=True,
                        key=f"active_{preset_id}_{target.value}_{item.metric}_{direction}",
                        help="Inactive rules are omitted completely from the submitted recipe.",
                    )
                    row[1].markdown(f"**{definition.display_name}**  \n`{item.metric}`")
                    value = float(
                        row[2].number_input(
                            direction.title(),
                            value=float(default),
                            key=f"value_{preset_id}_{target.value}_{item.metric}_{direction}",
                            help=_metric_help(filterable, item.metric),
                            **_filter_number_input_options(definition),
                        )
                    )
                    if active:
                        submitted_constraints.append(
                            MetricConstraint(
                                item.metric,
                                minimum=value if direction == "minimum" else None,
                                maximum=value if direction == "maximum" else None,
                            )
                        )
        if advanced:
            with st.expander("Advanced custom bounds", expanded=True):
                for metric in advanced:
                    definition = filterable.loc[filterable.name.eq(metric)].iloc[0]
                    values = pd.to_numeric(frame[metric], errors="coerce")
                    observed_min = float(values.min()) if values.notna().any() else 0.0
                    observed_max = float(values.max()) if values.notna().any() else 1.0
                    input_options = _filter_number_input_options(definition)
                    if "min_value" in input_options:
                        observed_min = max(observed_min, float(input_options["min_value"]))
                        observed_max = max(observed_max, float(input_options["min_value"]))
                    if "max_value" in input_options:
                        observed_min = min(observed_min, float(input_options["max_value"]))
                        observed_max = min(observed_max, float(input_options["max_value"]))
                    st.markdown(f"**{definition.display_name}** · `{metric}`")
                    columns = st.columns(4)
                    use_min = columns[0].checkbox(
                        "Use minimum",
                        key=f"advanced_min_active_{metric}",
                        help="Require values greater than or equal to the entered minimum.",
                    )
                    minimum = float(
                        columns[1].number_input(
                            "Minimum",
                            value=observed_min,
                            key=f"advanced_min_value_{metric}",
                            help=_metric_help(filterable, metric),
                            **input_options,
                        )
                    )
                    use_max = columns[2].checkbox(
                        "Use maximum",
                        key=f"advanced_max_active_{metric}",
                        help="Require values less than or equal to the entered maximum.",
                    )
                    maximum = float(
                        columns[3].number_input(
                            "Maximum",
                            value=observed_max,
                            key=f"advanced_max_value_{metric}",
                            help=_metric_help(filterable, metric),
                            **input_options,
                        )
                    )
                    if use_min or use_max:
                        submitted_constraints.append(
                            MetricConstraint(
                                metric,
                                minimum=minimum if use_min else None,
                                maximum=maximum if use_max else None,
                            )
                        )
        st.code(_constraint_expression(tuple(submitted_constraints)), language=None)
        apply_col, clear_col = st.columns(2)
        submitted = apply_col.form_submit_button("Apply filters", type="primary", width="stretch")
        cleared = clear_col.form_submit_button("Clear filter result", width="stretch")

    if cleared:
        for key in ("filter_recipe", "filter_result"):
            st.session_state.pop(key, None)
        st.rerun()
    if submitted:
        ordered_constraints = tuple(submitted_constraints)
        validate_metric_constraints(
            list(ordered_constraints),
            list(metric_definitions),
        )
        matches_preset = _constraints_match(ordered_constraints, resolved.constraints)
        recipe = FilterRecipe(
            constraints=ordered_constraints,
            conservation_target=target,
            preset_id=preset.identifier if matches_preset else f"custom_from:{preset.identifier}",
            preset_version=preset.version,
        )
        result = apply_filters(
            frame,
            recipe,
            filterable_metrics=filterable.name.tolist(),
            metric_definitions=metric_definitions,
        )
        st.session_state.filter_recipe = recipe
        st.session_state.filter_result = result

    result = st.session_state.get("filter_result")
    if result is None:
        st.markdown(
            "<div class='status-note'>Review the editable rules and press <strong>Apply "
            "filters</strong>. Changing controls alone does not alter a previous result.</div>",
            unsafe_allow_html=True,
        )
        return
    conservation_metric = result.recipe.conservation_target.value
    ranked = canonical_filter_rows(result.rows, conservation_metric)
    cards = st.columns(4)
    cards[0].metric("Initial windows", format_count(result.initial_count))
    cards[1].metric("Retained windows", format_count(result.retained_count))
    cards[2].metric("Removed windows", format_count(result.initial_count - result.retained_count))
    cards[3].metric(
        "Retained fraction",
        format_fraction(
            result.retained_count / result.initial_count if result.initial_count else 0.0
        ),
    )
    evidence = st.columns(4)
    coordinate_span = (
        f"{int(ranked.aa_start.min())}-{int(ranked.aa_end.max())}" if not ranked.empty else "-"
    )
    evidence[0].metric("Retained coordinate span", coordinate_span)
    evidence[1].metric(
        "Median call rate",
        format_fraction(ranked.mean_sample_call_rate.median()) if not ranked.empty else "-",
    )
    evidence[2].metric(
        "Median conservation",
        format_fraction(ranked[conservation_metric].median()) if not ranked.empty else "-",
    )
    evidence[3].metric(
        "Median complete peptide observations",
        format_count(ranked.n_exact_peptide_observations.median()) if not ranked.empty else "-",
    )
    st.caption(
        f"Submitted recipe: `{result.recipe.preset_id}` · conservation field "
        f"`{result.recipe.conservation_target.value}` · null policy "
        f"`{result.recipe.null_policy}`"
    )
    steps = result.steps.copy()
    if not steps.empty:
        steps.insert(
            2,
            "rule",
            [
                f"{row.metric} ≥ {row.minimum:g}"
                if pd.notna(row.minimum)
                else f"{row.metric} ≤ {row.maximum:g}"
                for row in steps.itertuples()
            ],
        )
        st.markdown("#### How each rule reduced the candidate set")
        st.dataframe(steps, hide_index=True, width="stretch")

    if result.rows.empty:
        st.warning("No peptide passes every active constraint. Try a looser preset or fewer rules.")
    else:
        from aa_metrics_app.charts.results import (
            diagnostic_frame,
            matplotlib_coordinate_strip,
            matplotlib_diagnostic_scatter,
            plotly_coordinate_strip,
            plotly_diagnostic_scatter,
        )

        diagnostics_section = content.section("filter_diagnostics")
        st.markdown(f"#### {diagnostics_section.heading}")
        st.markdown(section_intro_html(diagnostics_section), unsafe_allow_html=True)
        coordinate_diagnostic = diagnostic_frame(
            frame,
            ranked,
            conservation_metric=conservation_metric,
            comparison_metric="mean_sample_call_rate",
        )
        if selection.renderer is Renderer.PLOTLY:
            _plotly_chart(
                st,
                plotly_coordinate_strip(coordinate_diagnostic),
                key="filter_coordinate_strip_plotly",
            )
        else:
            st.pyplot(matplotlib_coordinate_strip(coordinate_diagnostic), width="stretch")

        st.markdown("##### Conservation vs. peptide evidence")
        st.markdown(
            "Choose the supporting metric shown on the vertical axis. The horizontal axis keeps "
            "the conservation target used by the submitted filter recipe."
        )
        comparison_metric = st.selectbox(
            "Diagnostic comparison metric",
            (
                "mean_sample_call_rate",
                "peptide_diversity_hd",
                "mean_sample_missing_fraction",
                "fraction_sample_peptides_exact",
            ),
            format_func=lambda name: registry.loc[registry.name.eq(name), "display_name"].iloc[0],
            help="Compare the submitted conservation target with quality, diversity, or evidence.",
        )
        diagnostic = diagnostic_frame(
            frame,
            ranked,
            conservation_metric=conservation_metric,
            comparison_metric=comparison_metric,
        )
        if selection.renderer is Renderer.PLOTLY:
            _plotly_chart(
                st,
                plotly_diagnostic_scatter(
                    diagnostic,
                    conservation_metric=conservation_metric,
                    comparison_metric=comparison_metric,
                ),
                key="filter_diagnostic_scatter_plotly",
            )
        else:
            st.pyplot(
                matplotlib_diagnostic_scatter(
                    diagnostic,
                    conservation_metric=conservation_metric,
                    comparison_metric=comparison_metric,
                ),
                width="stretch",
            )
        pinned = [
            "display_key",
            "stratum_type",
            "stratum_value",
            "block_number",
            "block_name",
            "window_size",
            "aa_start",
            "aa_end",
            "reference_sequence",
            "allelic_consensus_sequence",
            "all_call_consensus_sequence",
            "top_exact_peptide_sequence",
            "modal_all_call_sequence",
            "n_samples",
            "n_fully_callable_peptide_observations",
            "n_exact_peptide_observations",
            *[item.metric for item in result.recipe.constraints],
            conservation_metric,
            "mean_sample_call_rate",
            "mean_sample_missing_fraction",
            "mean_sample_mixed_fraction_callable",
            "peptide_diversity_hd",
        ]
        remaining = [name for name in ranked.columns if name not in pinned]
        complete = ranked.loc[
            :, [name for name in dict.fromkeys([*pinned, *remaining]) if name in ranked]
        ]
        st.markdown("#### Complete retained-peptide table")
        st.dataframe(
            complete,
            hide_index=True,
            width="stretch",
            height=min(620, 110 + 35 * len(complete)),
            column_config=_column_config(complete.columns, registry),
        )
        inspected = _region_sequence_inspector(ranked, content)
        _residue_logo(
            path,
            selection,
            int(inspected.aa_start),
            int(inspected.aa_end),
            context_key="filter",
            content=content,
            wrap=selection.window_size,
        )
        _structure_context(
            path,
            selection,
            int(inspected.aa_start),
            int(inspected.aa_end),
            str(inspected.reference_sequence),
            inspected.top_exact_peptide_sequence,
            context_key="filter",
            content=content,
        )
        _mutation_explorer(
            path,
            selection,
            int(inspected.aa_start),
            int(inspected.aa_end),
            context_key="filter",
            content=content,
        )

    identity = load_identity(path)
    selection_payload = {
        "stratum_type": selection.stratum_type.value,
        "stratum_value": selection.stratum_value,
        "block_number": selection.block_number,
        "window_size": selection.window_size,
    }
    summary_payload = recipe_manifest(
        result.recipe,
        selection=selection_payload,
        database_identity={
            "schema_version": identity.schema_version,
            "artifact_sha256": identity.artifact_sha256,
            "build_timestamp": identity.build_timestamp,
        },
        application_version="0.1.0",
    )
    summary_payload.update(
        {
            "initial_count": result.initial_count,
            "retained_count": result.retained_count,
            "removed_count": result.initial_count - result.retained_count,
            "display_keys": ranked.display_key.tolist(),
            "steps": json.loads(result.steps.to_json(orient="records")),
        }
    )
    bundle = build_filter_bundle(
        ranked,
        registry,
        summary=summary_payload,
        database_identity=identity,
        filename_context={
            **selection_payload,
            "preset_id": result.recipe.preset_id or "custom",
        },
    )
    st.markdown("#### Download this result")
    download_columns = st.columns(2)
    download_columns[0].download_button(
        "Download retained peptides (.zip)",
        bundle.content,
        file_name=bundle.filename,
        mime=bundle.mime_type,
        width="stretch",
        help=(
            "Canonical bundle with complete data, recipe, definitions, database identity, "
            "and README."
        ),
    )
    download_columns[1].download_button(
        "Download retained peptides (.csv)",
        safe_csv(ranked),
        file_name=bundle.filename.removesuffix(".zip") + ".csv",
        mime="text/csv",
        width="stretch",
        help="Convenience table only; use the ZIP when reproducibility metadata matters.",
    )
    with st.expander("Submitted rules and cumulative accounting"):
        st.code(_constraint_expression(result.recipe.constraints), language=None)
        st.dataframe(result.steps, hide_index=True, width="stretch")
        st.markdown(
            "Comparisons are inclusive. Nulls fail an active rule and are counted at the step "
            "where they are first encountered. Step exclusions depend on rule order; the final "
            "AND result does not. Missingness is evidence against confident conservation, not "
            "direct proof of an insertion or deletion."
        )


def _about(
    path: str,
    identity: DatabaseIdentity,
    blocks: tuple[BlockOption, ...],
    content: AppContent,
) -> None:
    _section_header(content, "about")
    st.markdown(
        """
### Interpreting results and limitations

- Complete peptide frequencies exclude missing, unresolved, and mixed calls. Always
  inspect how many complete peptide observations support them (excessive missing data could be an
  indication of poor mapping to the 3D7 reference because the sequence is too divergent).
- Mixed calls contribute fractionally (e.g., 0.5 for each allele if two alleles are observed) to
  site-level allele frequencies and diversity, but not to
  complete peptide frequencies (i.e., haplotypes with mixed calls are not considered complete).
- Missingness and low callability weaken evidence for conservation. They are not direct evidence
  of an insertion/deletion but could indicate poor data quality or a sequence too dissimilar from
  the 3D7 reference.
- Country and Population are alternative sample groupings (i.e., samples in a country could be
  part of two different populations).
- Overlapping windows are statistically dependent.
- Bear in mind the classical disclaimer: a high conservation score alone does not
  establish epitope presentation, immunogenicity, structural accessibility, or protection.
  
  
### Indels
This prototype doesn't consider indels (to be included in future iterations) but below you can see
a very conservative preliminary summary of indel burden for MSP1 derived from Pf7 data (the
previous MalariaGEN data release).
**Here we highlight intervals in the gene for which aggregated indel frequency was above 1% in
African samples** and we included the associated 3D7 reference sequence for that interval.
**Being very conservative, we advise that you don't include peptides that overlap these intervals** or are within any of the
gene blocks documented to be very divergent and enriched for indels (e.g., block 2).

- **57–58** (block_01): `EI`
- **61–96** (block_01 / block_02): `KGASAQSGASAQSGASAQSGASAQSGASAQSGASAQ`
- **99–104** (block_02): `TSGPSG`
- **125–131** (block_03): `GASPPAD`
- **335–335** (block_04): `D`
- **337–338** (block_04): `IN`
- **342–345** (block_04): `NPPP`
- **357–357** (block_05): `D`
- **740–753** (block_09): `GETEVTEETEETEE`
- **757–758** (block_09): `HT`
- **765–771** (block_09 / block_10): `LPPTQPS`
- **1471–1471** (block_16): `E`
  
### Database overview
"""
    )
    genes = load_gene_metadata(path)
    statistics = load_database_statistics(path)
    peptide_rows = load_peptide_rows_by_stratum(path)
    statistic_values = statistics.set_index("statistic").value

    columns = st.columns(6)
    columns[0].metric("Database schema", identity.schema_version)
    columns[1].metric("Storage", identity.storage_backend.replace("_", " + "))
    columns[2].metric("Dataset size", format_bytes(identity.size_bytes))
    columns[3].metric("Blocks", format_count(statistic_values.get("Blocks")))
    columns[4].metric("Peptide windows", format_count(statistic_values.get("Peptide windows")))

    st.markdown("### Gene metadata")
    st.markdown(
        "This prototype only focuses on the MSP1 gene. The table below shows the gene "
        "identifier and reference sequence."
    )
    gene_view = genes.loc[
        :,
        [
            "gene_identifier",
            "gene_name",
            "reference_name",
            "aa_start",
            "aa_end",
            "aa_length",
            "reference_sequence",
        ],
    ].rename(
        columns={
            "gene_identifier": "Gene identifier",
            "gene_name": "Gene name",
            "reference_name": "Reference",
            "aa_start": "AA start",
            "aa_end": "AA end",
            "aa_length": "AA length",
            "reference_sequence": "3D7/reference AA sequence",
        }
    )
    st.dataframe(gene_view, hide_index=True, width="stretch")

    st.markdown("### Block coordinates and reference sequences")
    block_view = pd.DataFrame.from_records(
        [
            {
                "Block": f"{item.block_number:02d}",
                "Name": item.name,
                "AA start": item.aa_start,
                "AA end": item.aa_end,
                "Length": item.length,
                "3D7/reference AA sequence": item.reference_sequence,
            }
            for item in blocks
        ]
    )
    st.dataframe(block_view, hide_index=True, width="stretch")

    st.markdown("### Database inventory")
    st.markdown(
        "Tables and entities in the DB (for debugging purposes)."
    )
    st.dataframe(
        statistics.rename(
            columns={"statistic": "Statistic", "value": "Rows", "description": "Meaning"}
        ),
        hide_index=True,
        width="stretch",
    )

    st.markdown("### Peptide metric rows by stratum")
    st.markdown(
        "Each row below is one stored independent Continent, Population, or Country cohort. "
        "Peptide metric "
        "rows count the stored stratum-specific summaries available to plots and filters."
    )
    st.dataframe(
        peptide_rows.rename(
            columns={
                "stratum_type": "Stratum type",
                "stratum_value": "Stratum value",
                "display_label": "Stratum",
                "n_samples": "QC-pass samples",
                "peptide_metric_rows": "Peptide metric rows",
                "physical_windows": "Physical windows",
                "peptide_lengths": "Peptide lengths",
                "minimum_length": "Minimum length",
                "maximum_length": "Maximum length",
            }
        ),
        hide_index=True,
        width="stretch",
        height=min(720, 75 + 35 * len(peptide_rows)),
    )
    st.info(
        "The database stores aggregate metrics and sequence summaries, not sample-level "
        "identifiers. For exploring complete haplotypes interactively go to "
        "[https://pf-haploatlas.streamlit.app/](https://pf-haploatlas.streamlit.app/)",
        icon="ℹ️",
    )


def main() -> None:
    """Render the complete application shell."""

    try:
        app_content = load_app_content(APP_CONTENT_PATH)
    except (FileNotFoundError, ValueError) as error:
        st.set_page_config(
            page_title="Amino-acid conservation explorer",
            page_icon="🧬",
            layout="wide",
            initial_sidebar_state="auto",
        )
        LOGGER.error("Cannot load application content: %s", type(error).__name__)
        st.error("The application content configuration is unavailable or invalid.", icon="🚫")
        st.stop()

    st.set_page_config(
        page_title=app_content.title,
        page_icon=app_content.page_icon,
        layout="wide",
        initial_sidebar_state="auto",
    )
    if app_content.logo_path is not None:
        logo_icon = (
            app_content.logo_icon_path.read_bytes()
            if app_content.logo_icon_path is not None
            else None
        )
        st.logo(
            app_content.logo_path.read_bytes(),
            icon_image=logo_icon,
        )
    _load_css()
    try:
        path = _data_source_token()
        identity = load_identity(path)
        strata = load_strata(path)
        blocks = load_blocks(path)
        metric_help = load_metric_help_catalog(str(METRIC_HELP_PATH), path)
        selection, stratum, block, navigation = _sidebar(path, strata, blocks, app_content)
    except (FileNotFoundError, DatabaseCompatibilityError, ValueError) as error:
        LOGGER.error("Cannot initialize conservation explorer: %s", type(error).__name__)
        st.error("The metrics database is unavailable or incompatible.", icon="🚫")
        st.markdown(
            "Set `AA_METRICS_SOURCE_CONFIG` to a valid source JSON file, then restart the app. "
            "The legacy `AA_METRICS_DB`, `AA_METRICS_DB_SIZE`, and `AA_METRICS_DB_SHA256` "
            "variables remain supported during migration. No sample-level data are required."
        )
        with st.expander("Technical summary"):
            st.code(f"{type(error).__name__}: startup validation failed")
            st.caption("See server logs for operational details; local paths are not shown here.")
        st.stop()
    _header(app_content, navigation, identity, stratum, block, selection.window_size)
    if navigation == "Block overview":
        _block_overview(path, selection, stratum, block, app_content)
    elif navigation == "Locus tracks":
        _locus_tracks(path, selection, block, app_content)
    elif navigation == "Window tracks":
        _window_tracks(path, selection, block, app_content)
    elif navigation == "Explore a region":
        _region_explorer(path, selection, block, app_content)
    elif navigation == "Filter & export":
        _filter_explorer(path, selection, block, app_content)
    elif navigation == "About & limitations":
        _about(path, identity, blocks, app_content)
    metric_help_level = {
        "Block overview": "block",
        "Locus tracks": "locus",
        "Window tracks": "window",
        "Explore a region": "window",
        "Filter & export": "window",
    }.get(navigation)
    if metric_help_level is not None:
        _metric_guide(metric_help, metric_help_level, app_content)


if __name__ == "__main__":
    main()
