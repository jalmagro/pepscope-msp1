"""Typed, deployment-owned data-source configuration for the metrics app."""

from __future__ import annotations

import json
import os
import re
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

SOURCE_CONFIG_ENV = "AA_METRICS_SOURCE_CONFIG"
LEGACY_DATABASE_ENV = "AA_METRICS_DB"
LEGACY_DATABASE_SIZE_ENV = "AA_METRICS_DB_SIZE"
LEGACY_DATABASE_SHA_ENV = "AA_METRICS_DB_SHA256"
SUPPORTED_SOURCE_SCHEMA_VERSION = "1.0.0"


class DataBackend(str, Enum):
    """Storage engines that can satisfy the application repository contract."""

    SQLITE = "sqlite"
    PARQUET_DUCKDB = "parquet_duckdb"


def _parse_memory_bytes(value: str, *, label: str) -> int:
    match = re.fullmatch(r"([1-9][0-9]*)(MB|GB)", value.strip().upper())
    if match is None:
        raise ValueError(f"{label} must use an integer MB or GB value")
    multiplier = 1_000_000 if match.group(2) == "MB" else 1_000_000_000
    return int(match.group(1)) * multiplier


@dataclass(frozen=True, slots=True)
class DuckDBRuntimeConfig:
    """Conservative DuckDB limits for a low-memory Streamlit deployment."""

    memory_limit: str = "256MB"
    threads: int = 1
    temp_directory: str = "/tmp/aa_metrics_duckdb"
    max_temp_directory_size: str = "1GB"

    def __post_init__(self) -> None:
        memory_bytes = _parse_memory_bytes(self.memory_limit, label="memory_limit")
        if not 64_000_000 <= memory_bytes <= 512_000_000:
            raise ValueError("memory_limit must be between 64MB and 512MB")
        if self.threads not in {1, 2}:
            raise ValueError("threads must be 1 or 2")
        temp_path = Path(self.temp_directory)
        if not temp_path.is_absolute() or not temp_path.is_relative_to(Path("/tmp")):
            raise ValueError("temp_directory must be an absolute path below /tmp")
        max_temp_bytes = _parse_memory_bytes(
            self.max_temp_directory_size, label="max_temp_directory_size"
        )
        if not 64_000_000 <= max_temp_bytes <= 2_000_000_000:
            raise ValueError(
                "max_temp_directory_size must be between 64MB and 2GB"
            )


@dataclass(frozen=True, slots=True)
class DataSourceConfig:
    """Resolved artifact location and optional release-verification metadata."""

    backend: DataBackend
    artifact: Path
    config_path: Path | None = None
    artifact_root: Path | None = None
    expected_size_bytes: int | None = None
    expected_sha256: str | None = None
    expected_total_size_bytes: int | None = None
    expected_dataset_sha256: str | None = None
    duckdb: DuckDBRuntimeConfig = DuckDBRuntimeConfig()

    def __post_init__(self) -> None:
        if self.expected_size_bytes is not None and self.expected_size_bytes < 1:
            raise ValueError("expected_size_bytes must be positive")
        if (
            self.expected_total_size_bytes is not None
            and self.expected_total_size_bytes < 1
        ):
            raise ValueError("expected_total_size_bytes must be positive")
        for label, value in (
            ("expected_sha256", self.expected_sha256),
            ("expected_dataset_sha256", self.expected_dataset_sha256),
        ):
            if value is not None and re.fullmatch(r"[0-9a-fA-F]{64}", value) is None:
                raise ValueError(f"{label} must contain 64 hexadecimal characters")
        if self.backend is DataBackend.SQLITE:
            if self.expected_total_size_bytes is not None or self.expected_dataset_sha256:
                raise ValueError("SQLite sources cannot use dataset-level expectations")
        elif self.expected_size_bytes is not None or self.expected_sha256:
            raise ValueError("Parquet sources cannot use single-file expectations")

    @property
    def resolved_artifact_root(self) -> Path:
        """Return the configured asset root, or the artifact's containing directory."""

        return (self.artifact_root or self.artifact.parent).resolve()


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _resolve_relative(value: object, base: Path, *, label: str) -> Path:
    text = str(value).strip()
    if not text:
        raise ValueError(f"{label} cannot be empty")
    path = Path(text)
    return (path if path.is_absolute() else base / path).resolve()


def load_data_source_config(path: str | Path) -> DataSourceConfig:
    """Load a source JSON file and resolve all relative paths from that file."""

    source = Path(path).resolve()
    payload = json.loads(
        source.read_text(encoding="utf-8"), object_pairs_hook=_unique_json_object
    )
    if not isinstance(payload, dict):
        raise ValueError("Data source configuration must be a JSON object")
    if payload.get("schema_version") != SUPPORTED_SOURCE_SCHEMA_VERSION:
        raise ValueError("Unsupported data-source schema version")
    backend = DataBackend(str(payload["backend"]))
    artifact = _resolve_relative(payload["artifact"], source.parent, label="artifact")
    artifact_root_value = payload.get("artifact_root")
    artifact_root = (
        _resolve_relative(artifact_root_value, source.parent, label="artifact_root")
        if artifact_root_value is not None
        else None
    )
    duckdb_payload = payload.get("duckdb", {})
    if not isinstance(duckdb_payload, dict):
        raise ValueError("duckdb must be a JSON object")
    duckdb = DuckDBRuntimeConfig(
        memory_limit=str(duckdb_payload.get("memory_limit", "256MB")),
        threads=int(duckdb_payload.get("threads", 1)),
        temp_directory=str(
            duckdb_payload.get("temp_directory", "/tmp/aa_metrics_duckdb")
        ),
        max_temp_directory_size=str(
            duckdb_payload.get("max_temp_directory_size", "1GB")
        ),
    )
    return DataSourceConfig(
        backend=backend,
        artifact=artifact,
        config_path=source,
        artifact_root=artifact_root,
        expected_size_bytes=(
            int(payload["expected_size_bytes"])
            if payload.get("expected_size_bytes") is not None
            else None
        ),
        expected_sha256=payload.get("expected_sha256"),
        expected_total_size_bytes=(
            int(payload["expected_total_size_bytes"])
            if payload.get("expected_total_size_bytes") is not None
            else None
        ),
        expected_dataset_sha256=payload.get("expected_dataset_sha256"),
        duckdb=duckdb,
    )


def data_source_from_environment(
    default_sqlite_path: str | Path,
    *,
    environment: Mapping[str, str] | None = None,
) -> DataSourceConfig:
    """Resolve preferred source-config, legacy SQLite variables, then the default path."""

    values = os.environ if environment is None else environment
    source_path = values.get(SOURCE_CONFIG_ENV)
    if source_path:
        return load_data_source_config(source_path)
    database_path = Path(values.get(LEGACY_DATABASE_ENV, default_sqlite_path)).resolve()
    size_value = values.get(LEGACY_DATABASE_SIZE_ENV)
    return DataSourceConfig(
        backend=DataBackend.SQLITE,
        artifact=database_path,
        expected_size_bytes=int(size_value) if size_value else None,
        expected_sha256=values.get(LEGACY_DATABASE_SHA_ENV),
    )

