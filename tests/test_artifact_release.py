"""Validate the immutable standalone database and its companion structures."""

from __future__ import annotations

import hashlib
from pathlib import Path

from aa_metrics_app import (
    MetricsRepository,
    load_data_source_config,
    load_parquet_manifest,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "app/config/data-source.json"
EXPECTED_DIGEST = "97e46d356b4b035946e7391289a15674d131c6031b0093c4243b55246a2adf61"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_qualified_parquet_artifact_is_complete_and_bounded() -> None:
    source = load_data_source_config(CONFIG)
    assert source.backend.value == "parquet_duckdb"
    assert source.expected_total_size_bytes == 59_627_648
    assert source.expected_dataset_sha256 == EXPECTED_DIGEST

    manifest = load_parquet_manifest(source.artifact)
    shards = [shard for table in manifest.tables for shard in table.files]
    expected_paths = {shard.path for shard in shards}
    actual_paths = {
        path.relative_to(source.artifact.parent).as_posix()
        for path in source.artifact.parent.rglob("*.parquet")
    }

    assert len(manifest.tables) == 23
    assert len(shards) == 73
    assert manifest.total_size_bytes == 59_627_648
    assert manifest.dataset_sha256 == EXPECTED_DIGEST
    assert max(shard.size_bytes for shard in shards) == 4_743_526
    assert all(shard.size_bytes <= 20_000_000 for shard in shards)
    assert actual_paths == expected_paths

    for shard in shards:
        path = source.artifact.parent / shard.path
        assert path.stat().st_size == shard.size_bytes
        assert _sha256(path) == shard.sha256


def test_release_counts_and_default_stratum() -> None:
    repository = MetricsRepository(load_data_source_config(CONFIG))
    statistics = repository.database_statistics().set_index("statistic")["value"]
    strata = repository.strata()

    assert statistics["Genes"] == 1
    assert statistics["Blocks"] == 17
    assert statistics["Loci"] == 1_720
    assert statistics["Peptide windows"] == 34_140
    assert statistics["Strata"] == 26
    assert statistics["Peptide metric rows"] == 887_640
    assert len(strata) == 26
    assert [(item.stratum_type.value, item.value) for item in strata if item.is_default] == [
        ("continent", "Africa")
    ]


def test_all_database_referenced_structures_are_packaged_and_verified() -> None:
    source = load_data_source_config(CONFIG)
    repository = MetricsRepository(source)
    candidates = repository.structures_for_interval(1, 1720)
    assets = candidates.loc[:, ["asset_relative_path", "asset_sha256"]].drop_duplicates()

    assert len(assets) == 6
    for row in assets.itertuples(index=False):
        path = source.resolved_artifact_root / row.asset_relative_path
        assert path.is_file()
        assert _sha256(path) == row.asset_sha256
