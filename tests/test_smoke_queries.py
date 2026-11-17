"""Exercise representative read-only application queries and exports."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from aa_metrics_app import (
    MetricConstraint,
    MetricsRepository,
    RegionQuery,
    build_region_bundle,
    load_data_source_config,
    select_region,
)

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = MetricsRepository(load_data_source_config(ROOT / "app/config/data-source.json"))


def test_representative_block_locus_window_and_mutation_queries() -> None:
    block = REPOSITORY.block_summary("continent", "Africa", block_number=17)
    locus = REPOSITORY.locus_track("continent", "Africa", block_number=17)
    windows = REPOSITORY.window_track(
        "continent", "Africa", block_number=17, window_size=15
    )
    mutations = REPOSITORY.mutation_comparison(
        1644, 1716, "continent", "Africa", include_missing=True
    )
    spectrum = REPOSITORY.locus_allele_spectrum_interval(
        "continent", "Africa", 1644, 1716
    )

    assert len(block) == 1
    assert locus.aa_position.is_monotonic_increasing
    assert len(locus) == int(block.aa_length.iloc[0])
    assert not windows.empty and windows.aa_start.is_monotonic_increasing
    assert {"mutation", "frequency__continent/Africa", "missing__continent/Africa"} <= set(
        mutations.columns
    )
    assert not spectrum.empty


def test_region_filter_and_in_memory_export() -> None:
    windows = REPOSITORY.window_track(
        "continent", "Africa", block_number=17, window_size=15
    )
    coordinate = int(windows.aa_start.iloc[len(windows) // 2])
    selected = select_region(
        windows,
        RegionQuery(coordinate, coordinate),
        block_start=int(windows.aa_start.min()),
        block_end=int(windows.aa_end.max()),
        window_size=15,
    )
    filtered = REPOSITORY.filter_windows(
        "continent",
        "Africa",
        block_number=17,
        window_size=15,
        constraints=[MetricConstraint("mean_sample_call_rate", minimum=0.0, maximum=1.0)],
    )
    definitions = pd.DataFrame(
        [
            {
                "name": definition.name,
                "display_name": definition.display_name,
                "description": definition.description,
                "unit": definition.unit,
            }
            for definition in REPOSITORY.metric_definitions("window")
        ]
    )
    bundle = build_region_bundle(
        selected.rows,
        definitions,
        query={"start": coordinate, "end": coordinate},
        database_identity=REPOSITORY.identity(include_sha256=True),
        filename_context={
            "stratum_type": "continent",
            "stratum_value": "Africa",
            "block_number": 17,
            "window_size": 15,
        },
    )

    assert not selected.rows.empty
    assert len(filtered) == len(windows)
    assert bundle.filename.endswith(".zip")
    assert bundle.content.startswith(b"PK")


def test_structure_query_and_mapping_are_available() -> None:
    structures = REPOSITORY.structures_for_interval(1626, 1720)
    assert not structures.empty
    selected = structures.iloc[0]
    mapping = REPOSITORY.structure_residue_map(
        int(selected.structure_chain_id), 1626, 1720
    )
    assert not mapping.empty
    assert mapping.aa_position.between(1626, 1720).all()
