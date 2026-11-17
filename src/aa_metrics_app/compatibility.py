"""Database compatibility and identity checks for the standalone application."""

from __future__ import annotations

import hashlib
import re
import sqlite3
from pathlib import Path

from .models import DatabaseIdentity

SUPPORTED_SCHEMA_VERSION = "2.0.0"
SUPPORTED_SQLITE_USER_VERSION = 20_000

REQUIRED_COLUMNS: dict[str, set[str]] = {
    "analysis_run": {"schema_version", "created_at_utc", "status", "package_version"},
    "gene": {"gene_identifier", "gene_name", "aa_start", "aa_end", "aa_length"},
    "block": {"block_number", "aa_start", "aa_end", "aa_length"},
    "stratum": {
        "stratum_type", "stratum_value", "n_samples", "is_default", "sort_order"
    },
    "metric_definition": {"entity_level", "column_name", "description"},
    "v_block_metrics": {"stratum_type", "stratum_value", "block_number"},
    "v_locus_metrics": {"stratum_type", "stratum_value", "block_number", "aa_position"},
    "v_window_metrics": {
        "stratum_type",
        "stratum_value",
        "block_number",
        "window_size",
        "aa_start",
        "aa_end",
    },
    "block_sequence_frequency": {"stratum_id", "block_id", "sequence_scope"},
    "block_callability_histogram": {"stratum_id", "block_id", "bin_index"},
    "window_metric_bounds": {"stratum_id", "window_size", "metric_name"},
    "mutation": {"canonical_mutation_id", "aa_position", "alternate_residue"},
    "mutation_metric": {"stratum_id", "mutation_id", "mutation_allele_frequency"},
    "structure": {"accession", "asset_relative_path", "asset_sha256"},
    "structure_residue_map": {"structure_chain_id", "project_locus_id"},
}


class DatabaseCompatibilityError(RuntimeError):
    """Raised when an artifact cannot safely serve the application."""


def open_readonly(path: str | Path) -> sqlite3.Connection:
    """Open one immutable, read-only SQLite connection.

    Parameters
    ----------
    path
        Existing database artifact.

    Returns
    -------
    sqlite3.Connection
        Row-oriented connection owned by the caller.
    """

    target = Path(path).resolve()
    if not target.is_file():
        raise FileNotFoundError(f"Metrics database does not exist: {target}")
    try:
        connection = sqlite3.connect(
            f"file:{target.as_posix()}?mode=ro&immutable=1",
            uri=True,
            timeout=10,
        )
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA query_only = ON")
        # SELECT 1 does not force SQLite to parse the database header. Reading schema metadata does.
        connection.execute("PRAGMA schema_version").fetchone()
        connection.execute("SELECT 1").fetchone()
        return connection
    except sqlite3.DatabaseError as error:
        raise DatabaseCompatibilityError(f"Cannot open metrics database: {error}") from error


def _table_columns(connection: sqlite3.Connection, table: str) -> set[str]:
    return {str(row[1]) for row in connection.execute(f'PRAGMA table_info("{table}")')}


def validate_database(path: str | Path) -> None:
    """Validate the schema and scientific contract required by the app."""

    with open_readonly(path) as connection:
        user_version = int(connection.execute("PRAGMA user_version").fetchone()[0])
        if user_version != SUPPORTED_SQLITE_USER_VERSION:
            raise DatabaseCompatibilityError(
                f"Unsupported SQLite user_version {user_version}; "
                f"expected {SUPPORTED_SQLITE_USER_VERSION}"
            )
        run = connection.execute(
            "SELECT schema_version, status FROM analysis_run ORDER BY run_id DESC LIMIT 1"
        ).fetchone()
        if run is None:
            raise DatabaseCompatibilityError("Database has no analysis run metadata")
        if run["schema_version"] != SUPPORTED_SCHEMA_VERSION:
            raise DatabaseCompatibilityError(
                f"Unsupported schema version {run['schema_version']!r}; "
                f"expected {SUPPORTED_SCHEMA_VERSION!r}"
            )
        if run["status"] != "validated":
            raise DatabaseCompatibilityError("Database has not passed validation")
        for table, required in REQUIRED_COLUMNS.items():
            columns = _table_columns(connection, table)
            if not columns:
                raise DatabaseCompatibilityError(f"Required table or view is missing: {table}")
            missing = required.difference(columns)
            if missing:
                raise DatabaseCompatibilityError(
                    f"{table} is missing required columns: {sorted(missing)}"
                )
        default_rows = connection.execute(
            "SELECT COUNT(*) FROM stratum WHERE is_default = 1"
        ).fetchone()[0]
        if default_rows != 1:
            raise DatabaseCompatibilityError("Database must contain exactly one default stratum")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_artifact(
    path: str | Path,
    *,
    expected_size: int | None = None,
    expected_sha256: str | None = None,
) -> str | None:
    """Validate an immutable database artifact against optional release metadata.

    Size is checked before the more expensive checksum. The function never modifies or replaces
    the supplied artifact and returns the computed digest only when one was requested.
    """

    target = Path(path).resolve()
    validate_database(target)
    if expected_size is not None:
        if expected_size < 1:
            raise ValueError("Expected database size must be a positive integer")
        actual_size = target.stat().st_size
        if actual_size != expected_size:
            raise DatabaseCompatibilityError(
                f"Database size mismatch: expected {expected_size} bytes, observed {actual_size}"
            )
    if expected_sha256 is None:
        return None
    expected = expected_sha256.strip().lower()
    if not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ValueError("Expected database SHA-256 must contain exactly 64 hexadecimal characters")
    actual = _sha256(target)
    if actual != expected:
        raise DatabaseCompatibilityError("Database SHA-256 does not match the configured release")
    return actual


def database_identity(
    path: str | Path,
    *,
    include_sha256: bool = False,
    artifact_root: str | Path | None = None,
    cache_identity: str | None = None,
) -> DatabaseIdentity:
    """Return artifact identity and build provenance after compatibility validation.

    SHA-256 calculation is optional because the production database is large. Deployment should
    compare a configured release checksum once at startup rather than during every page rerun.
    """

    target = Path(path).resolve()
    validate_database(target)
    with open_readonly(target) as connection:
        run = connection.execute(
            "SELECT * FROM analysis_run ORDER BY run_id DESC LIMIT 1"
        ).fetchone()
        gene = connection.execute(
            "SELECT gene_identifier, gene_name FROM gene ORDER BY gene_id LIMIT 1"
        ).fetchone()
        sources = tuple(
            dict(row)
            for row in connection.execute(
                "SELECT logical_role, relative_path, size_bytes, sha256, row_count "
                "FROM source_file WHERE run_id = ? ORDER BY logical_role, relative_path",
                (int(run["run_id"]),),
            )
        )
    return DatabaseIdentity(
        path=str(target),
        size_bytes=target.stat().st_size,
        schema_version=str(run["schema_version"]),
        sqlite_user_version=SUPPORTED_SQLITE_USER_VERSION,
        build_timestamp=str(run["created_at_utc"]),
        build_status=str(run["status"]),
        package_version=str(run["package_version"]),
        gene_identifier=str(gene["gene_identifier"]),
        gene_name=None if gene["gene_name"] is None else str(gene["gene_name"]),
        artifact_sha256=_sha256(target) if include_sha256 else None,
        source_files=sources,
        storage_backend="sqlite",
        artifact_root=str(Path(artifact_root).resolve()) if artifact_root else str(target.parent),
        cache_identity=cache_identity,
    )
