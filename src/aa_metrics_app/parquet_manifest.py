"""Versioned manifest contract for sharded Parquet metric artifacts."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

ARTIFACT_FORMAT_VERSION = "1.0.0"
PARQUET_STORAGE_BACKEND = "parquet_duckdb"
DEFAULT_MAX_FILE_BYTES = 20_000_000


def _require_sha256(value: str, *, label: str) -> str:
    normalized = value.strip().lower()
    if re.fullmatch(r"[0-9a-f]{64}", normalized) is None:
        raise ValueError(f"{label} must contain 64 hexadecimal characters")
    return normalized


def _safe_relative_path(value: str) -> str:
    if "\\" in value:
        raise ValueError("Manifest shard paths must use POSIX separators")
    path = PurePosixPath(value)
    if not value or path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError(f"Unsafe manifest shard path: {value!r}")
    if path.as_posix() != value:
        raise ValueError(f"Manifest shard path is not canonical: {value!r}")
    return value


@dataclass(frozen=True, slots=True)
class ParquetShardManifest:
    """Integrity, size, and ordering metadata for one Parquet shard."""

    path: str
    row_count: int
    size_bytes: int
    sha256: str
    key_min: tuple[Any, ...] | None = None
    key_max: tuple[Any, ...] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", _safe_relative_path(self.path))
        object.__setattr__(self, "sha256", _require_sha256(self.sha256, label="sha256"))
        if self.row_count < 0:
            raise ValueError("Shard row_count cannot be negative")
        if self.size_bytes < 1:
            raise ValueError("Shard size_bytes must be positive")
        if (self.key_min is None) != (self.key_max is None):
            raise ValueError("Shard key_min and key_max must be provided together")
        if self.key_min is not None and len(self.key_min) != len(self.key_max or ()):
            raise ValueError("Shard key bounds must have equal arity")

    def as_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "path": self.path,
            "row_count": self.row_count,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
        }
        if self.key_min is not None:
            result["key_min"] = list(self.key_min)
            result["key_max"] = list(self.key_max or ())
        return result


@dataclass(frozen=True, slots=True)
class ParquetTableManifest:
    """Logical schema and shards for one persisted relational table."""

    name: str
    row_count: int
    column_names: tuple[str, ...]
    column_types: tuple[str, ...]
    primary_key: tuple[str, ...]
    partition_columns: tuple[str, ...]
    schema_fingerprint: str
    files: tuple[ParquetShardManifest, ...]

    def __post_init__(self) -> None:
        if not self.name or not self.name.replace("_", "").isalnum():
            raise ValueError(f"Unsafe table name: {self.name!r}")
        if self.row_count < 0:
            raise ValueError("Table row_count cannot be negative")
        if not self.column_names or len(set(self.column_names)) != len(self.column_names):
            raise ValueError("Table column_names must be non-empty and unique")
        if len(self.column_types) != len(self.column_names) or any(
            value not in {"INTEGER", "REAL", "TEXT", "BLOB"}
            for value in self.column_types
        ):
            raise ValueError(
                "Table column_types must align with columns and use supported logical types"
            )
        if not set(self.primary_key).issubset(self.column_names):
            raise ValueError("Table primary_key must reference declared columns")
        if not set(self.partition_columns).issubset(self.column_names):
            raise ValueError("Table partition_columns must reference declared columns")
        if not self.schema_fingerprint.strip():
            raise ValueError("Table schema_fingerprint cannot be empty")
        if len({item.path for item in self.files}) != len(self.files):
            raise ValueError(f"Table {self.name!r} contains duplicate shard paths")
        if sum(item.row_count for item in self.files) != self.row_count:
            raise ValueError(f"Table {self.name!r} shard row counts do not match row_count")
        if self.row_count and not self.files:
            raise ValueError(f"Non-empty table {self.name!r} has no shards")

    def as_dict(self) -> dict[str, Any]:
        return {
            "row_count": self.row_count,
            "column_names": list(self.column_names),
            "column_types": list(self.column_types),
            "primary_key": list(self.primary_key),
            "partition_columns": list(self.partition_columns),
            "schema_fingerprint": self.schema_fingerprint,
            "files": [item.as_dict() for item in self.files],
        }


def dataset_digest(tables: tuple[ParquetTableManifest, ...]) -> str:
    """Return a deterministic digest of sorted shard path/size/hash records."""

    records = sorted(
        (shard.path, shard.size_bytes, shard.sha256)
        for table in tables
        for shard in table.files
    )
    payload = json.dumps(records, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True, slots=True)
class ParquetArtifactManifest:
    """Authoritative description of one immutable Parquet dataset."""

    scientific_schema_version: str
    dataset_id: str
    created_at_utc: str
    duckdb_min_version: str
    compression: str
    compression_level: int
    max_file_bytes: int
    tables: tuple[ParquetTableManifest, ...]
    total_size_bytes: int | None = None
    dataset_sha256: str | None = None
    artifact_format_version: str = ARTIFACT_FORMAT_VERSION
    storage_backend: str = PARQUET_STORAGE_BACKEND
    build_profile_id: str | None = None
    default_stratum_type: str | None = None
    default_stratum_value: str | None = None
    gene_identifier: str | None = None
    structure_asset_root: str | None = None

    def __post_init__(self) -> None:
        if self.artifact_format_version != ARTIFACT_FORMAT_VERSION:
            raise ValueError("Unsupported Parquet artifact format version")
        if self.storage_backend != PARQUET_STORAGE_BACKEND:
            raise ValueError("Manifest storage_backend must be parquet_duckdb")
        for label, value in (
            ("scientific_schema_version", self.scientific_schema_version),
            ("dataset_id", self.dataset_id),
            ("created_at_utc", self.created_at_utc),
            ("duckdb_min_version", self.duckdb_min_version),
        ):
            if not value.strip():
                raise ValueError(f"{label} cannot be empty")
        if self.compression != "zstd":
            raise ValueError("The first Parquet artifact version requires zstd compression")
        if not 1 <= self.compression_level <= 22:
            raise ValueError("compression_level must be between 1 and 22")
        if self.max_file_bytes < 64 * 1024:
            raise ValueError("max_file_bytes is too small")
        names = [item.name for item in self.tables]
        if not names or len(names) != len(set(names)):
            raise ValueError("Manifest tables must be non-empty and unique")
        paths = [shard.path for table in self.tables for shard in table.files]
        if len(paths) != len(set(paths)):
            raise ValueError("Manifest shard paths must be unique across all tables")
        computed_size = sum(
            shard.size_bytes for table in self.tables for shard in table.files
        )
        if self.total_size_bytes is None:
            object.__setattr__(self, "total_size_bytes", computed_size)
        elif self.total_size_bytes != computed_size:
            raise ValueError("Manifest total_size_bytes does not match shard sizes")
        computed_digest = dataset_digest(self.tables)
        if self.dataset_sha256 is None:
            object.__setattr__(self, "dataset_sha256", computed_digest)
        else:
            normalized_digest = _require_sha256(
                self.dataset_sha256, label="dataset_sha256"
            )
            if normalized_digest != computed_digest:
                raise ValueError("Manifest dataset_sha256 does not match shard records")
            object.__setattr__(self, "dataset_sha256", normalized_digest)
        for shard in (item for table in self.tables for item in table.files):
            if shard.size_bytes > self.max_file_bytes:
                raise ValueError(f"Parquet shard exceeds max_file_bytes: {shard.path}")
        if self.structure_asset_root is not None:
            _safe_relative_path(self.structure_asset_root)

    def as_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "artifact_format_version": self.artifact_format_version,
            "scientific_schema_version": self.scientific_schema_version,
            "storage_backend": self.storage_backend,
            "dataset_id": self.dataset_id,
            "created_at_utc": self.created_at_utc,
            "duckdb_min_version": self.duckdb_min_version,
            "compression": self.compression,
            "compression_level": self.compression_level,
            "max_file_bytes": self.max_file_bytes,
            "total_size_bytes": self.total_size_bytes,
            "dataset_sha256": self.dataset_sha256,
            "tables": {item.name: item.as_dict() for item in self.tables},
        }
        for key in (
            "build_profile_id",
            "default_stratum_type",
            "default_stratum_value",
            "gene_identifier",
            "structure_asset_root",
        ):
            value = getattr(self, key)
            if value is not None:
                result[key] = value
        return result

    def validate_paths(self, artifact_root: str | Path, *, require_files: bool = True) -> None:
        """Confirm every listed shard remains below the artifact root."""

        root = Path(artifact_root).resolve()
        for table in self.tables:
            for shard in table.files:
                target = (root / shard.path).resolve()
                if not target.is_relative_to(root):
                    raise ValueError(f"Shard escapes artifact root: {shard.path}")
                if require_files and not target.is_file():
                    raise FileNotFoundError(f"Manifest shard is missing: {shard.path}")


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def manifest_from_dict(payload: dict[str, Any]) -> ParquetArtifactManifest:
    """Parse and validate a manifest mapping."""

    tables_payload = payload.get("tables")
    if not isinstance(tables_payload, dict):
        raise ValueError("Manifest tables must be a JSON object")
    tables: list[ParquetTableManifest] = []
    for name, table_payload in sorted(tables_payload.items()):
        if not isinstance(table_payload, dict):
            raise ValueError(f"Manifest table {name!r} must be a JSON object")
        files_payload = table_payload.get("files", [])
        if not isinstance(files_payload, list):
            raise ValueError(f"Manifest table {name!r} files must be a JSON list")
        files = tuple(
            ParquetShardManifest(
                path=str(item["path"]),
                row_count=int(item["row_count"]),
                size_bytes=int(item["size_bytes"]),
                sha256=str(item["sha256"]),
                key_min=(
                    tuple(item["key_min"]) if item.get("key_min") is not None else None
                ),
                key_max=(
                    tuple(item["key_max"]) if item.get("key_max") is not None else None
                ),
            )
            for item in files_payload
        )
        tables.append(
            ParquetTableManifest(
                name=str(name),
                row_count=int(table_payload["row_count"]),
                column_names=tuple(str(value) for value in table_payload["column_names"]),
                column_types=tuple(str(value) for value in table_payload["column_types"]),
                primary_key=tuple(str(value) for value in table_payload["primary_key"]),
                partition_columns=tuple(
                    str(value) for value in table_payload["partition_columns"]
                ),
                schema_fingerprint=str(table_payload["schema_fingerprint"]),
                files=files,
            )
        )
    return ParquetArtifactManifest(
        artifact_format_version=str(payload["artifact_format_version"]),
        scientific_schema_version=str(payload["scientific_schema_version"]),
        storage_backend=str(payload["storage_backend"]),
        dataset_id=str(payload["dataset_id"]),
        created_at_utc=str(payload["created_at_utc"]),
        duckdb_min_version=str(payload["duckdb_min_version"]),
        compression=str(payload["compression"]),
        compression_level=int(payload["compression_level"]),
        max_file_bytes=int(payload["max_file_bytes"]),
        total_size_bytes=int(payload["total_size_bytes"]),
        dataset_sha256=str(payload["dataset_sha256"]),
        build_profile_id=payload.get("build_profile_id"),
        default_stratum_type=payload.get("default_stratum_type"),
        default_stratum_value=payload.get("default_stratum_value"),
        gene_identifier=payload.get("gene_identifier"),
        structure_asset_root=payload.get("structure_asset_root"),
        tables=tuple(tables),
    )


def load_parquet_manifest(path: str | Path) -> ParquetArtifactManifest:
    """Load an artifact manifest while rejecting duplicate JSON keys."""

    payload = json.loads(
        Path(path).read_text(encoding="utf-8"), object_pairs_hook=_unique_json_object
    )
    if not isinstance(payload, dict):
        raise ValueError("Parquet artifact manifest must be a JSON object")
    return manifest_from_dict(payload)


def write_parquet_manifest(
    manifest: ParquetArtifactManifest, path: str | Path
) -> None:
    """Write canonical, human-readable JSON. Call only after all shards validate."""

    Path(path).write_text(
        json.dumps(manifest.as_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
