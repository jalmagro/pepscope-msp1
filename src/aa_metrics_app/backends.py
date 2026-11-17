"""Backend-neutral query sessions and the read-only SQLite implementation."""

from __future__ import annotations

import hashlib
import re
import sqlite3
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

import pandas as pd

from .compatibility import (
    REQUIRED_COLUMNS,
    SUPPORTED_SCHEMA_VERSION,
    DatabaseCompatibilityError,
    database_identity,
    open_readonly,
    verify_artifact,
)
from .data_sources import DataBackend, DataSourceConfig
from .models import DatabaseIdentity
from .parquet_manifest import load_parquet_manifest
from .views import SCHEMA_V2_VIEW_SQL


@runtime_checkable
class QuerySession(Protocol):
    """Small result-conversion boundary shared by supported SQL engines."""

    def frame(
        self, sql: str, parameters: Sequence[object] = ()
    ) -> pd.DataFrame: ...

    def rows(
        self, sql: str, parameters: Sequence[object] = ()
    ) -> tuple[dict[str, Any], ...]: ...

    def row(
        self, sql: str, parameters: Sequence[object] = ()
    ) -> dict[str, Any] | None: ...

    def scalar(self, sql: str, parameters: Sequence[object] = ()) -> Any: ...


@runtime_checkable
class QueryBackend(Protocol):
    """Artifact validation, identity, inspection, and query-session contract."""

    @property
    def artifact_root(self) -> Path: ...

    @property
    def artifact_path(self) -> Path: ...

    @property
    def cache_identity(self) -> str: ...

    @property
    def storage_backend(self) -> str: ...

    def validate(self) -> None: ...

    def identity(self, *, include_sha256: bool = False) -> DatabaseIdentity: ...

    def table_columns(self, name: str) -> tuple[str, ...]: ...

    def float64_cast(self, expression: str) -> str: ...

    def session(self) -> Any: ...


class SQLiteQuerySession:
    """Engine-neutral result conversion over one owned SQLite connection."""

    def __init__(self, connection: sqlite3.Connection) -> None:
        self._connection = connection

    @staticmethod
    def _columns(cursor: sqlite3.Cursor) -> list[str]:
        return [str(item[0]) for item in cursor.description or ()]

    def frame(
        self, sql: str, parameters: Sequence[object] = ()
    ) -> pd.DataFrame:
        cursor = self._connection.execute(sql, tuple(parameters))
        columns = self._columns(cursor)
        return pd.DataFrame.from_records(
            [tuple(row) for row in cursor.fetchall()], columns=columns
        )

    def rows(
        self, sql: str, parameters: Sequence[object] = ()
    ) -> tuple[dict[str, Any], ...]:
        cursor = self._connection.execute(sql, tuple(parameters))
        columns = self._columns(cursor)
        return tuple(dict(zip(columns, tuple(row), strict=True)) for row in cursor.fetchall())

    def row(
        self, sql: str, parameters: Sequence[object] = ()
    ) -> dict[str, Any] | None:
        cursor = self._connection.execute(sql, tuple(parameters))
        columns = self._columns(cursor)
        result = cursor.fetchone()
        return None if result is None else dict(zip(columns, tuple(result), strict=True))

    def scalar(self, sql: str, parameters: Sequence[object] = ()) -> Any:
        cursor = self._connection.execute(sql, tuple(parameters))
        result = cursor.fetchone()
        if result is None:
            return None
        return result[0]


class SQLiteQueryBackend:
    """Immutable, short-lived SQLite query backend."""

    def __init__(self, source: DataSourceConfig | str | Path) -> None:
        if isinstance(source, DataSourceConfig):
            if source.backend is not DataBackend.SQLITE:
                raise ValueError("SQLiteQueryBackend requires a SQLite data source")
            self.source = source
        else:
            self.source = DataSourceConfig(
                backend=DataBackend.SQLITE,
                artifact=Path(source).resolve(),
            )
        self.path = self.source.artifact.resolve()

    @property
    def artifact_root(self) -> Path:
        return self.source.resolved_artifact_root

    @property
    def artifact_path(self) -> Path:
        return self.path

    @property
    def storage_backend(self) -> str:
        return DataBackend.SQLITE.value

    @property
    def cache_identity(self) -> str:
        digest = self.source.expected_sha256
        if digest:
            return f"sqlite:{digest.lower()}"
        if not self.path.is_file():
            return f"sqlite:{self.path.as_posix()}:missing"
        stat = self.path.stat()
        return f"sqlite:{self.path.as_posix()}:{stat.st_size}:{stat.st_mtime_ns}"

    def validate(self) -> None:
        verify_artifact(
            self.path,
            expected_size=self.source.expected_size_bytes,
            expected_sha256=self.source.expected_sha256,
        )

    def identity(self, *, include_sha256: bool = False) -> DatabaseIdentity:
        return database_identity(
            self.path,
            include_sha256=include_sha256,
            artifact_root=self.artifact_root,
            cache_identity=self.cache_identity,
        )

    def table_columns(self, name: str) -> tuple[str, ...]:
        if not name or not name.replace("_", "").isalnum():
            raise ValueError("Unsafe table name")
        escaped = name.replace('"', '""')
        with open_readonly(self.path) as connection:
            rows = connection.execute(f'PRAGMA table_info("{escaped}")').fetchall()
        return tuple(str(row[1]) for row in rows)

    def float64_cast(self, expression: str) -> str:
        return f"CAST({expression} AS REAL)"

    @contextmanager
    def session(self) -> Iterator[SQLiteQuerySession]:
        connection = open_readonly(self.path)
        try:
            yield SQLiteQuerySession(connection)
        finally:
            connection.close()


class DuckDBParquetQuerySession:
    """Engine-neutral result conversion over one owned DuckDB connection."""

    def __init__(self, connection: Any) -> None:
        self._connection = connection

    @staticmethod
    def _columns(cursor: Any) -> list[str]:
        return [str(item[0]) for item in cursor.description or ()]

    def frame(
        self, sql: str, parameters: Sequence[object] = ()
    ) -> pd.DataFrame:
        cursor = self._connection.execute(sql, list(parameters))
        columns = self._columns(cursor)
        return pd.DataFrame.from_records(cursor.fetchall(), columns=columns)

    def rows(
        self, sql: str, parameters: Sequence[object] = ()
    ) -> tuple[dict[str, Any], ...]:
        cursor = self._connection.execute(sql, list(parameters))
        columns = self._columns(cursor)
        return tuple(dict(zip(columns, row, strict=True)) for row in cursor.fetchall())

    def row(
        self, sql: str, parameters: Sequence[object] = ()
    ) -> dict[str, Any] | None:
        cursor = self._connection.execute(sql, list(parameters))
        columns = self._columns(cursor)
        result = cursor.fetchone()
        return None if result is None else dict(zip(columns, result, strict=True))

    def scalar(self, sql: str, parameters: Sequence[object] = ()) -> Any:
        result = self._connection.execute(sql, list(parameters)).fetchone()
        return None if result is None else result[0]


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _version_tuple(value: str) -> tuple[int, ...]:
    numbers = re.match(r"^(\d+(?:\.\d+)*)", value)
    if numbers is None:
        raise ValueError(f"Cannot parse version: {value!r}")
    return tuple(int(item) for item in numbers.group(1).split("."))


def _sql_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _sql_identifier(value: str) -> str:
    if not value or not value.replace("_", "").isalnum():
        raise ValueError(f"Unsafe SQL identifier: {value!r}")
    return '"' + value.replace('"', '""') + '"'


_DUCKDB_TYPES = {
    "INTEGER": "BIGINT",
    "REAL": "DOUBLE",
    "TEXT": "VARCHAR",
    "BLOB": "BLOB",
}


class DuckDBParquetQueryBackend:
    """Lazy, in-memory DuckDB catalog over exact manifest-listed Parquet shards."""

    def __init__(self, source: DataSourceConfig) -> None:
        if source.backend is not DataBackend.PARQUET_DUCKDB:
            raise ValueError("DuckDBParquetQueryBackend requires a Parquet data source")
        self.source = source
        self.manifest_path = source.artifact.resolve()
        self.manifest = load_parquet_manifest(self.manifest_path)

    @staticmethod
    def _duckdb():
        try:
            import duckdb
        except ImportError as error:
            raise RuntimeError(
                "DuckDB is required for a parquet_duckdb data source"
            ) from error
        return duckdb

    @property
    def artifact_root(self) -> Path:
        return self.source.resolved_artifact_root

    @property
    def artifact_path(self) -> Path:
        return self.manifest_path

    @property
    def storage_backend(self) -> str:
        return DataBackend.PARQUET_DUCKDB.value

    @property
    def cache_identity(self) -> str:
        return f"parquet_duckdb:{self.manifest.dataset_sha256}"

    def _validate_release_expectations(self) -> None:
        if (
            self.source.expected_total_size_bytes is not None
            and self.source.expected_total_size_bytes != self.manifest.total_size_bytes
        ):
            raise DatabaseCompatibilityError(
                "Parquet dataset size does not match the configured release"
            )
        if self.source.expected_dataset_sha256 is not None and (
            self.source.expected_dataset_sha256.lower()
            != self.manifest.dataset_sha256
        ):
            raise DatabaseCompatibilityError(
                "Parquet dataset SHA-256 does not match the configured release"
            )

    def _validate_shards(self) -> None:
        self.manifest.validate_paths(self.manifest_path.parent)
        for table in self.manifest.tables:
            for shard in table.files:
                path = self.manifest_path.parent / shard.path
                if path.stat().st_size != shard.size_bytes:
                    raise DatabaseCompatibilityError(
                        f"Parquet shard size mismatch: {shard.path}"
                    )
                if _file_sha256(path) != shard.sha256:
                    raise DatabaseCompatibilityError(
                        f"Parquet shard checksum mismatch: {shard.path}"
                    )

    def _apply_runtime_settings(self, connection: Any) -> None:
        runtime = self.source.duckdb
        Path(runtime.temp_directory).mkdir(parents=True, exist_ok=True)
        connection.execute(f"SET memory_limit = {_sql_string(runtime.memory_limit)}")
        connection.execute(f"SET threads = {runtime.threads}")
        connection.execute(
            f"SET temp_directory = {_sql_string(runtime.temp_directory)}"
        )
        connection.execute(
            "SET max_temp_directory_size = "
            + _sql_string(runtime.max_temp_directory_size)
        )

    def _register_table_views(self, connection: Any) -> None:
        for table in self.manifest.tables:
            name = _sql_identifier(table.name)
            if table.files:
                paths = ", ".join(
                    _sql_string(str((self.manifest_path.parent / shard.path).resolve()))
                    for shard in table.files
                )
                connection.execute(
                    f"CREATE VIEW {name} AS SELECT * FROM read_parquet(["
                    f"{paths}], union_by_name=false, hive_partitioning="
                    f"{'true' if table.partition_columns else 'false'})"
                )
            else:
                columns = ", ".join(
                    f"CAST(NULL AS {_DUCKDB_TYPES[data_type]}) AS "
                    f"{_sql_identifier(column)}"
                    for column, data_type in zip(
                        table.column_names, table.column_types, strict=True
                    )
                )
                connection.execute(
                    f"CREATE VIEW {name} AS SELECT {columns} WHERE FALSE"
                )
        for statement in SCHEMA_V2_VIEW_SQL:
            connection.execute(statement)

    def _connect(self) -> Any:
        duckdb = self._duckdb()
        if _version_tuple(duckdb.__version__) < _version_tuple(
            self.manifest.duckdb_min_version
        ):
            raise DatabaseCompatibilityError(
                f"DuckDB {duckdb.__version__} is older than required "
                f"{self.manifest.duckdb_min_version}"
            )
        connection = duckdb.connect(":memory:")
        try:
            self._apply_runtime_settings(connection)
            self._register_table_views(connection)
        except Exception:
            connection.close()
            raise
        return connection

    def validate(self) -> None:
        if not self.manifest_path.is_file():
            raise FileNotFoundError(
                f"Parquet artifact manifest does not exist: {self.manifest_path}"
            )
        if self.manifest.scientific_schema_version != SUPPORTED_SCHEMA_VERSION:
            raise DatabaseCompatibilityError(
                "Unsupported Parquet scientific schema version "
                f"{self.manifest.scientific_schema_version!r}"
            )
        self._validate_release_expectations()
        self._validate_shards()
        connection = self._connect()
        try:
            physical_tables = int(
                connection.execute(
                    "SELECT COUNT(*) FROM duckdb_tables() WHERE internal = FALSE"
                ).fetchone()[0]
            )
            if physical_tables:
                raise DatabaseCompatibilityError(
                    "Parquet backend must register lazy views, not physical tables"
                )
            for table in self.manifest.tables:
                quoted = _sql_identifier(table.name)
                description = connection.execute(f"DESCRIBE SELECT * FROM {quoted}").fetchall()
                names = tuple(str(row[0]) for row in description)
                types = tuple(str(row[1]) for row in description)
                expected_types = tuple(_DUCKDB_TYPES[item] for item in table.column_types)
                if names != table.column_names or types != expected_types:
                    raise DatabaseCompatibilityError(
                        f"Parquet table schema drift: {table.name}"
                    )
                observed = int(
                    connection.execute(f"SELECT COUNT(*) FROM {quoted}").fetchone()[0]
                )
                if observed != table.row_count:
                    raise DatabaseCompatibilityError(
                        f"Parquet table row-count drift: {table.name}"
                    )
            for relation, required in REQUIRED_COLUMNS.items():
                columns = {
                    str(row[0])
                    for row in connection.execute(
                        f"DESCRIBE SELECT * FROM {_sql_identifier(relation)}"
                    ).fetchall()
                }
                missing = required.difference(columns)
                if missing:
                    raise DatabaseCompatibilityError(
                        f"{relation} is missing required columns: {sorted(missing)}"
                    )
            run = connection.execute(
                "SELECT schema_version, status FROM analysis_run "
                "ORDER BY run_id DESC LIMIT 1"
            ).fetchone()
            if run is None or run[0] != SUPPORTED_SCHEMA_VERSION or run[1] != "validated":
                raise DatabaseCompatibilityError(
                    "Parquet database has no compatible validated analysis run"
                )
            defaults = int(
                connection.execute(
                    "SELECT COUNT(*) FROM stratum WHERE is_default = 1"
                ).fetchone()[0]
            )
            if defaults != 1:
                raise DatabaseCompatibilityError(
                    "Database must contain exactly one default stratum"
                )
        finally:
            connection.close()

    def identity(self, *, include_sha256: bool = False) -> DatabaseIdentity:
        with self.session() as session:
            run = session.row(
                "SELECT * FROM analysis_run ORDER BY run_id DESC LIMIT 1"
            )
            gene = session.row(
                "SELECT gene_identifier, gene_name FROM gene ORDER BY gene_id LIMIT 1"
            )
            if run is None or gene is None:
                raise DatabaseCompatibilityError(
                    "Parquet artifact is missing identity metadata"
                )
            sources = session.rows(
                "SELECT logical_role, relative_path, size_bytes, sha256, row_count "
                "FROM source_file WHERE run_id = ? ORDER BY logical_role, relative_path",
                (int(run["run_id"]),),
            )
        return DatabaseIdentity(
            path=str(self.manifest_path),
            size_bytes=int(self.manifest.total_size_bytes or 0),
            schema_version=str(run["schema_version"]),
            sqlite_user_version=None,
            build_timestamp=str(run["created_at_utc"]),
            build_status=str(run["status"]),
            package_version=str(run["package_version"]),
            gene_identifier=str(gene["gene_identifier"]),
            gene_name=None if gene["gene_name"] is None else str(gene["gene_name"]),
            artifact_sha256=(
                str(self.manifest.dataset_sha256) if include_sha256 else None
            ),
            source_files=sources,
            storage_backend=self.storage_backend,
            artifact_root=str(self.artifact_root),
            cache_identity=self.cache_identity,
        )

    def table_columns(self, name: str) -> tuple[str, ...]:
        quoted = _sql_identifier(name)
        connection = self._connect()
        try:
            return tuple(
                str(row[0])
                for row in connection.execute(
                    f"DESCRIBE SELECT * FROM {quoted}"
                ).fetchall()
            )
        finally:
            connection.close()

    def float64_cast(self, expression: str) -> str:
        return f"CAST({expression} AS DOUBLE)"

    @contextmanager
    def session(self) -> Iterator[DuckDBParquetQuerySession]:
        connection = self._connect()
        try:
            yield DuckDBParquetQuerySession(connection)
        finally:
            connection.close()


def backend_from_source(source: DataSourceConfig | str | Path) -> QueryBackend:
    """Construct the selected backend without importing optional engines eagerly."""

    if isinstance(source, DataSourceConfig) and source.backend is DataBackend.PARQUET_DUCKDB:
        return DuckDBParquetQueryBackend(source)
    return SQLiteQueryBackend(source)
