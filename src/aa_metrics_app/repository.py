"""Short-lived, read-only queries for the conservation application."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from .backends import QueryBackend, QuerySession, backend_from_source
from .data_sources import DataSourceConfig
from .metric_catalog import build_metric_catalog
from .models import (
    BlockOption,
    DatabaseIdentity,
    MetricConstraint,
    MetricDefinition,
    StratumOption,
    StratumType,
    WindowSizeOption,
    validate_metric_constraints,
)

STRATUM_TYPES = {item.value for item in StratumType}
COORDINATE_ORDER_FIELDS = {
    "window_id": "window_id",
    "block_number": "block_number",
    "aa_start": "aa_start",
    "aa_end": "aa_end",
    "window_size": "window_size",
}


class MetricsRepository:
    """Query one immutable aggregate database with one connection per operation.

    Parameters
    ----------
    path
        Path, typed data-source descriptor, or backend satisfying the query contract.
    validate
        Validate the artifact at construction. Applications should leave this enabled and cache
        the successfully constructed repository configuration, not a live engine connection.
    """

    def __init__(
        self,
        path: DataSourceConfig | QueryBackend | str | Path,
        *,
        validate: bool = True,
    ) -> None:
        self.backend = (
            path if isinstance(path, QueryBackend) else backend_from_source(path)
        )
        self.path = self.backend.artifact_path
        if validate:
            self.backend.validate()

    @property
    def artifact_root(self) -> Path:
        """Root used to resolve companion structure assets."""

        return self.backend.artifact_root

    @property
    def cache_identity(self) -> str:
        """Stable backend-aware key for application data caches."""

        return self.backend.cache_identity

    @property
    def storage_backend(self) -> str:
        """Public storage implementation identifier."""

        return self.backend.storage_backend

    @staticmethod
    def _validate_stratum_type(stratum_type: str) -> None:
        if stratum_type not in STRATUM_TYPES:
            raise ValueError(
                "stratum_type must be global, continent, country, or population"
            )

    def _stratum_id(
        self,
        session: QuerySession,
        stratum_type: str,
        stratum_value: str,
    ) -> int:
        self._validate_stratum_type(stratum_type)
        row = session.row(
            "SELECT stratum_id FROM stratum WHERE stratum_type = ? AND stratum_value = ?",
            (stratum_type, stratum_value),
        )
        if row is None:
            raise ValueError(f"Unknown stratum {stratum_type!r}/{stratum_value!r}")
        return int(row["stratum_id"])

    def identity(self, *, include_sha256: bool = False) -> DatabaseIdentity:
        """Return build and artifact identity."""

        return self.backend.identity(include_sha256=include_sha256)

    def genes(self) -> pd.DataFrame:
        """Return gene metadata and decoded reference sequences."""

        with self.backend.session() as session:
            return session.frame(
                "SELECT g.*, s.sequence_text AS reference_sequence FROM gene AS g "
                "JOIN sequence AS s ON s.sequence_id = g.reference_sequence_id "
                "ORDER BY g.gene_id",
            )

    def database_statistics(self) -> pd.DataFrame:
        """Return ordered high-level row counts for the aggregate artifact."""

        definitions = (
            ("Genes", "gene", "Gene metadata records"),
            ("Blocks", "block", "Non-overlapping gene blocks"),
            ("Loci", "locus", "Reference amino-acid coordinates"),
            ("Peptide windows", "window", "Physical coordinate/length windows"),
            ("Strata", "stratum", "Stored independent geographic cohorts"),
            ("Block metric rows", "block_metric", "Stratum-specific block summaries"),
            ("Locus metric rows", "locus_metric", "Stratum-specific locus summaries"),
            ("Peptide metric rows", "window_metric", "Stratum-specific peptide summaries"),
            ("Stored sequences", "sequence", "Deduplicated sequence strings"),
            ("Mutations", "mutation", "Canonical non-reference residues"),
            ("Mutation metric rows", "mutation_metric", "Dense stratum mutation facts"),
            ("Structures", "structure", "Reviewed local structure models"),
            ("Structure residue maps", "structure_residue_map", "Project-coordinate mappings"),
        )
        with self.backend.session() as session:
            records = [
                {
                    "statistic": label,
                    "value": int(session.scalar(f'SELECT COUNT(*) FROM "{table}"')),
                    "description": description,
                }
                for label, table, description in definitions
            ]
            stratum_counts = {
                str(row["stratum_type"]): int(row["n_strata"])
                for row in session.rows(
                    "SELECT stratum_type, COUNT(*) AS n_strata "
                    "FROM stratum GROUP BY stratum_type"
                )
            }
        insert_at = 5
        for label, key in (
            ("Global strata", "global"),
            ("Continent strata", "continent"),
            ("Country strata", "country"),
            ("Population strata", "population"),
        ):
            records.insert(
                insert_at,
                {
                    "statistic": label,
                    "value": stratum_counts.get(key, 0),
                    "description": f"Stored {key} cohort definitions",
                },
            )
            insert_at += 1
        return pd.DataFrame.from_records(records)

    def peptide_rows_by_stratum(self) -> pd.DataFrame:
        """Return peptide-window metric coverage for every independent stratum."""

        with self.backend.session() as session:
            return session.frame(
                "SELECT s.stratum_type, s.stratum_value, s.display_label, s.n_samples, "
                "COUNT(wm.window_id) AS peptide_metric_rows, "
                "COUNT(DISTINCT wm.window_id) AS physical_windows, "
                "COUNT(DISTINCT w.window_size) AS peptide_lengths, "
                "MIN(w.window_size) AS minimum_length, MAX(w.window_size) AS maximum_length "
                "FROM stratum AS s LEFT JOIN window_metric AS wm "
                "ON wm.stratum_id = s.stratum_id LEFT JOIN \"window\" AS w "
                "ON w.window_id = wm.window_id "
                "GROUP BY s.stratum_id, s.stratum_type, s.stratum_value, "
                "s.display_label, s.n_samples, s.sort_order "
                "ORDER BY s.sort_order",
            )

    def strata(self, stratum_type: str | None = None) -> tuple[StratumOption, ...]:
        """Return ordered independent stratum choices."""

        parameters: tuple[object, ...] = ()
        where = ""
        if stratum_type is not None:
            self._validate_stratum_type(stratum_type)
            where = "WHERE stratum_type = ?"
            parameters = (stratum_type,)
        with self.backend.session() as session:
            rows = session.rows(
                f"SELECT * FROM stratum {where} ORDER BY sort_order",
                parameters,
            )
        return tuple(
            StratumOption(
                stratum_id=int(row["stratum_id"]),
                stratum_type=StratumType(str(row["stratum_type"])),
                value=str(row["stratum_value"]),
                label=str(row["display_label"]),
                n_samples=int(row["n_samples"]),
                small_sample_warning=bool(row["small_sample_warning"]),
                is_default=bool(row["is_default"]),
                sort_order=int(row["sort_order"]),
            )
            for row in rows
        )

    def blocks(self) -> tuple[BlockOption, ...]:
        """Return ordered block choices with decoded reference sequences."""

        with self.backend.session() as session:
            rows = session.rows(
                "SELECT b.*, s.sequence_text AS reference_sequence FROM block AS b "
                "JOIN sequence AS s ON s.sequence_id = b.reference_sequence_id "
                "ORDER BY b.gene_id, b.block_number"
            )
        return tuple(
            BlockOption(
                block_id=int(row["block_id"]),
                block_number=int(row["block_number"]),
                name=str(row["block_name"]),
                aa_start=int(row["aa_start"]),
                aa_end=int(row["aa_end"]),
                length=int(row["aa_length"]),
                reference_sequence=str(row["reference_sequence"]),
            )
            for row in rows
        )

    def window_sizes(self, block_number: int | None = None) -> tuple[WindowSizeOption, ...]:
        """Return available peptide lengths and row counts, optionally for one block."""

        parameters: tuple[object, ...] = ()
        where = ""
        if block_number is not None:
            if block_number < 1:
                raise ValueError("block_number must be positive")
            where = "WHERE b.block_number = ?"
            parameters = (block_number,)
        with self.backend.session() as session:
            rows = session.rows(
                "SELECT w.window_size, COUNT(*) AS n_windows FROM \"window\" AS w "
                "JOIN block AS b ON b.block_id = w.block_id "
                f"{where} GROUP BY w.window_size ORDER BY w.window_size",
                parameters,
            )
        return tuple(
            WindowSizeOption(
                window_size=int(row["window_size"]),
                n_windows=int(row["n_windows"]),
            )
            for row in rows
        )

    def metric_definitions(self, level: str | None = None) -> tuple[MetricDefinition, ...]:
        """Return authoritative metric definitions augmented with UI groups."""

        parameters: tuple[object, ...] = ()
        where = ""
        if level is not None:
            if level not in {"block", "locus", "window"}:
                raise ValueError("level must be block, locus, or window")
            where = "WHERE entity_level = ?"
            parameters = (level,)
        with self.backend.session() as session:
            frame = session.frame(
                f"SELECT * FROM metric_definition {where} "
                "ORDER BY entity_level, column_name",
                parameters,
            )
        return build_metric_catalog(frame)

    def mutations(
        self,
        aa_start: int,
        aa_end: int,
        *,
        stratum_type: str | None = None,
        stratum_value: str | None = None,
    ) -> pd.DataFrame:
        """Return dense mutation facts for an inclusive amino-acid interval.

        Supplying a stratum type and value limits the tidy result to that stored group. Omitting
        both returns every stored group. Zero-frequency facts are retained and remain distinct
        from unavailable frequencies caused by a zero callable-allele denominator.
        """

        if aa_start < 1 or aa_end < aa_start:
            raise ValueError("Mutation interval coordinates are invalid")
        if (stratum_type is None) != (stratum_value is None):
            raise ValueError("stratum_type and stratum_value must be provided together")
        where = ["mu.aa_position BETWEEN ? AND ?"]
        parameters: list[object] = [aa_start, aa_end]
        with self.backend.session() as session:
            if stratum_type is not None and stratum_value is not None:
                self._stratum_id(session, stratum_type, stratum_value)
                where.extend(("s.stratum_type = ?", "s.stratum_value = ?"))
                parameters.extend((stratum_type, stratum_value))
            return session.frame(
                "SELECT mu.canonical_mutation_id AS mutation, mu.aa_position, "
                "mu.reference_residue, mu.alternate_residue, b.block_number, "
                "s.stratum_type, s.stratum_value, s.display_label, s.sort_order, "
                "mm.n_samples, mm.allele_denominator_fractional, "
                "mm.alternate_fractional_count, mm.mutation_allele_frequency, "
                "mm.sample_presence_count, mm.sample_presence_fraction, mm.n_missing, "
                "mm.missing_fraction, mm.n_unresolved, mm.unresolved_fraction, "
                "mm.n_non_callable, mm.non_callable_fraction, mm.is_observed_in_stratum "
                "FROM mutation_metric AS mm "
                "JOIN mutation AS mu ON mu.mutation_id = mm.mutation_id "
                "JOIN block AS b ON b.block_id = mu.block_id "
                "JOIN stratum AS s ON s.stratum_id = mm.stratum_id WHERE "
                + " AND ".join(where)
                + " ORDER BY mu.aa_position, mu.alternate_residue, s.sort_order",
                parameters,
            )

    def mutation_comparison(
        self,
        aa_start: int,
        aa_end: int,
        stratum_type: str,
        stratum_value: str,
        *,
        include_missing: bool = False,
    ) -> pd.DataFrame:
        """Return one row per interval mutation with the selected group first."""

        with self.backend.session() as session:
            selected_id = self._stratum_id(session, stratum_type, stratum_value)
            selected = session.row(
                "SELECT stratum_type, display_label FROM stratum WHERE stratum_id = ?",
                (selected_id,),
            )
        if selected is None:
            raise ValueError(f"Unknown stratum {stratum_type!r}/{stratum_value!r}")
        tidy = self.mutations(aa_start, aa_end)
        if tidy.empty:
            return pd.DataFrame(
                columns=(
                    "mutation",
                    "aa_position",
                    "reference_residue",
                    "alternate_residue",
                    "block_number",
                )
            )
        tidy = tidy.copy()
        tidy["group"] = tidy.stratum_type + "/" + tidy.display_label
        selected_group = f"{selected['stratum_type']}/{selected['display_label']}"
        stored_order = tidy.sort_values("sort_order", kind="stable").group.drop_duplicates()
        group_order = [selected_group, *stored_order.loc[stored_order.ne(selected_group)].tolist()]
        index = [
            "mutation",
            "aa_position",
            "reference_residue",
            "alternate_residue",
            "block_number",
        ]
        frequencies = tidy.pivot(
            index=index, columns="group", values="mutation_allele_frequency"
        ).reindex(columns=group_order)
        result = frequencies.add_prefix("frequency__")
        if include_missing:
            missing = tidy.pivot(
                index=index, columns="group", values="missing_fraction"
            ).reindex(columns=group_order)
            result = result.join(missing.add_prefix("missing__"))
        result.columns.name = None
        result = result.reset_index()
        selected_frequency = f"frequency__{selected_group}"
        selected_missing = f"missing__{selected_group}"
        leading = ["mutation", selected_frequency]
        if selected_missing in result:
            leading.append(selected_missing)
        coordinates = [
            "aa_position",
            "reference_residue",
            "alternate_residue",
            "block_number",
        ]
        remaining = [column for column in result if column not in (*leading, *coordinates)]
        return result.loc[:, [*leading, *coordinates, *remaining]]

    def locus_allele_spectrum_interval(
        self,
        stratum_type: str,
        stratum_value: str,
        aa_start: int,
        aa_end: int,
    ) -> pd.DataFrame:
        """Return sample-normalized canonical residue weights for an inclusive interval."""

        if aa_start < 1 or aa_end < aa_start:
            raise ValueError("Allele-spectrum interval coordinates are invalid")
        missing_cast = self.backend.float64_cast("lm.n_missing")
        unresolved_cast = self.backend.float64_cast("lm.n_unresolved")
        with self.backend.session() as session:
            stratum_id = self._stratum_id(session, stratum_type, stratum_value)
            return session.frame(
                "SELECT l.aa_position, l.reference_residue, lm.n_samples, lm.call_rate, "
                f"{missing_cast} / lm.n_samples AS missing_fraction, "
                f"{unresolved_cast} / lm.n_samples AS unresolved_fraction, "
                "a.residue, a.fractional_count, "
                "a.fractional_count / lm.n_samples AS sample_normalized_height, "
                "a.fractional_frequency AS allele_frequency_callable "
                "FROM locus AS l JOIN locus_metric AS lm ON lm.locus_id = l.locus_id "
                "LEFT JOIN locus_allele_frequency AS a ON a.locus_id = l.locus_id "
                "AND a.stratum_id = lm.stratum_id WHERE lm.stratum_id = ? "
                "AND l.aa_position BETWEEN ? AND ? "
                "ORDER BY l.aa_position, a.residue",
                (stratum_id, aa_start, aa_end),
            )

    def structures_for_interval(self, aa_start: int, aa_end: int) -> pd.DataFrame:
        """Return reviewed structure-chain coverage for an inclusive project interval."""

        if aa_start < 1 or aa_end < aa_start:
            raise ValueError("Structure interval coordinates are invalid")
        length = aa_end - aa_start + 1
        count_expression = "COALESCE(mapped.mapped_residue_count_interval, 0)"
        count_cast = self.backend.float64_cast(count_expression)
        with self.backend.session() as session:
            frame = session.frame(
                "SELECT s.*, sc.structure_chain_id, sc.chain_id, sc.entity_id, "
                "sc.mapping_identity, sc.mapping_coverage, sc.mapping_status, "
                f"{count_expression} AS mapped_residue_count_interval, "
                f"{count_cast} / ? AS interval_coverage "
                "FROM structure AS s JOIN structure_chain AS sc "
                "ON sc.structure_id = s.structure_id "
                "LEFT JOIN (SELECT rm.structure_chain_id, "
                "COUNT(*) AS mapped_residue_count_interval "
                "FROM structure_residue_map AS rm JOIN locus AS l "
                "ON l.locus_id = rm.project_locus_id "
                "WHERE l.aa_position BETWEEN ? AND ? "
                "GROUP BY rm.structure_chain_id) AS mapped "
                "ON mapped.structure_chain_id = sc.structure_chain_id",
                (length, aa_start, aa_end),
            )
        if frame.empty:
            return frame
        frame["adequate_interval_coverage"] = frame.interval_coverage.ge(0.80)
        frame["selection_priority"] = 2
        frame.loc[
            frame.adequate_interval_coverage & frame.source_type.eq("experimental"),
            "selection_priority",
        ] = 0
        frame.loc[
            frame.adequate_interval_coverage & frame.source_type.eq("alphafold"),
            "selection_priority",
        ] = 1
        return frame.sort_values(
            [
                "selection_priority",
                "interval_coverage",
                "mapping_identity",
                "resolution_angstrom",
                "accession",
                "chain_id",
            ],
            ascending=[True, False, False, True, True, True],
            na_position="last",
            kind="stable",
            ignore_index=True,
        )

    def structure_residue_map(
        self,
        structure_chain_id: int,
        aa_start: int,
        aa_end: int,
    ) -> pd.DataFrame:
        """Return the reviewed project-to-structure map for one chain and interval."""

        if structure_chain_id < 1 or aa_start < 1 or aa_end < aa_start:
            raise ValueError("Structure map query values are invalid")
        with self.backend.session() as session:
            return session.frame(
                "SELECT s.accession, s.source_type, sc.chain_id, sc.entity_id, "
                "l.aa_position, rm.structure_residue_number, rm.insertion_code, "
                "rm.structure_residue, rm.project_reference_residue, rm.mapping_status, "
                "rm.confidence_value FROM structure_residue_map AS rm "
                "JOIN structure_chain AS sc "
                "ON sc.structure_chain_id = rm.structure_chain_id "
                "JOIN structure AS s ON s.structure_id = sc.structure_id "
                "JOIN locus AS l ON l.locus_id = rm.project_locus_id "
                "WHERE rm.structure_chain_id = ? AND l.aa_position BETWEEN ? AND ? "
                "ORDER BY l.aa_position",
                (structure_chain_id, aa_start, aa_end),
            )

    def block_summary(
        self,
        stratum_type: str,
        stratum_value: str,
        block_number: int | None = None,
    ) -> pd.DataFrame:
        """Return all block metrics for one stratum, optionally one block."""

        self._validate_stratum_type(stratum_type)
        if block_number is not None:
            block_where = " AND block_number = ?"
        else:
            block_where = ""
        with self.backend.session() as session:
            stratum_id = self._stratum_id(session, stratum_type, stratum_value)
            parameters: list[object] = [stratum_id]
            if block_number is not None:
                parameters.append(block_number)
            frame = session.frame(
                "SELECT * FROM v_block_metrics WHERE stratum_id = ?"
                + block_where
                + " ORDER BY block_number",
                parameters,
            )
        if frame.empty:
            raise ValueError("Unknown stratum or block selection")
        return frame

    def locus_track(
        self,
        stratum_type: str,
        stratum_value: str,
        block_number: int,
    ) -> pd.DataFrame:
        """Return one ordered locus row per coordinate in a selected block."""

        self._validate_stratum_type(stratum_type)
        with self.backend.session() as session:
            stratum_id = self._stratum_id(session, stratum_type, stratum_value)
            frame = session.frame(
                "SELECT * FROM v_locus_metrics WHERE stratum_id = ? "
                "AND block_number = ? ORDER BY aa_position",
                (stratum_id, block_number),
            )
        if frame.empty:
            raise ValueError("Unknown stratum or block selection")
        return frame

    def window_track(
        self,
        stratum_type: str,
        stratum_value: str,
        block_number: int,
        window_size: int,
    ) -> pd.DataFrame:
        """Return every stored metric for one stratum, block, and peptide length."""

        self._validate_stratum_type(stratum_type)
        if window_size < 1:
            raise ValueError("window_size must be positive")
        with self.backend.session() as session:
            stratum_id = self._stratum_id(session, stratum_type, stratum_value)
            frame = session.frame(
                "SELECT * FROM v_window_metrics WHERE stratum_id = ? "
                "AND block_number = ? AND window_size = ? "
                "ORDER BY aa_start, aa_end, window_id",
                (stratum_id, block_number, window_size),
            )
        if frame.empty:
            raise ValueError("Unknown stratum, block, or window-size selection")
        frame.insert(
            0,
            "display_key",
            frame.apply(
                lambda row: (
                    f"b{int(row.block_number):02d}:"
                    f"{int(row.aa_start)}-{int(row.aa_end)}:{int(row.window_size)}aa"
                ),
                axis=1,
            ),
        )
        return frame

    def window_metric_bounds(
        self,
        stratum_type: str,
        stratum_value: str,
        window_size: int,
    ) -> pd.DataFrame:
        """Return stored data-aware bounds and help for filterable window metrics."""

        with self.backend.session() as session:
            stratum_id = self._stratum_id(session, stratum_type, stratum_value)
            return session.frame(
                "SELECT b.*, d.display_name, d.description, d.unit, d.higher_is, "
                "d.display_decimals, d.denominator_definition, d.null_definition "
                "FROM window_metric_bounds AS b JOIN metric_definition AS d "
                "ON d.entity_level = 'window' AND d.column_name = b.metric_name "
                "WHERE b.stratum_id = ? AND b.window_size = ? ORDER BY b.metric_name",
                (stratum_id, window_size),
            )

    def block_sequence_frequencies(
        self,
        stratum_type: str,
        stratum_value: str,
        block_number: int,
        *,
        scope: str = "exact",
    ) -> pd.DataFrame:
        """Return ranked exact or normalized all-call block sequences."""

        if scope not in {"exact", "all_call"}:
            raise ValueError("scope must be exact or all_call")
        with self.backend.session() as session:
            stratum_id = self._stratum_id(session, stratum_type, stratum_value)
            return session.frame(
                "SELECT b.block_number, b.block_name, f.*, q.sequence_text AS sequence "
                "FROM block_sequence_frequency AS f "
                "JOIN block AS b ON b.block_id = f.block_id "
                "JOIN sequence AS q ON q.sequence_id = f.sequence_id "
                "WHERE f.stratum_id = ? AND b.block_number = ? AND f.sequence_scope = ? "
                "ORDER BY f.sequence_rank, f.sequence_id",
                (stratum_id, block_number, scope),
            )

    def block_callability_histogram(
        self,
        stratum_type: str,
        stratum_value: str,
        block_number: int,
    ) -> pd.DataFrame:
        """Return fixed-bin sample callability for one selected block."""

        with self.backend.session() as session:
            stratum_id = self._stratum_id(session, stratum_type, stratum_value)
            return session.frame(
                "SELECT b.block_number, b.block_name, h.* "
                "FROM block_callability_histogram AS h "
                "JOIN block AS b ON b.block_id = h.block_id "
                "WHERE h.stratum_id = ? AND b.block_number = ? ORDER BY h.bin_index",
                (stratum_id, block_number),
            )

    def filter_windows(
        self,
        stratum_type: str,
        stratum_value: str,
        block_number: int,
        window_size: int,
        constraints: Sequence[MetricConstraint],
    ) -> pd.DataFrame:
        """Apply validated inclusive constraints in the selected query backend."""

        definitions = {
            item.name: item for item in self.metric_definitions("window") if item.filterable
        }
        validate_metric_constraints(
            list(constraints),
            list(definitions.values()),
        )
        where = ["stratum_id = ?", "block_number = ?", "window_size = ?"]
        for constraint in constraints:
            if constraint.metric not in definitions:
                raise ValueError(f"Unknown filterable window metric: {constraint.metric}")
            if constraint.minimum is not None:
                where.append(f'"{constraint.metric}" >= ?')
            if constraint.maximum is not None:
                where.append(f'"{constraint.metric}" <= ?')
        with self.backend.session() as session:
            stratum_id = self._stratum_id(session, stratum_type, stratum_value)
            parameters: list[object] = [stratum_id, block_number, window_size]
            for constraint in constraints:
                if constraint.minimum is not None:
                    parameters.append(constraint.minimum)
                if constraint.maximum is not None:
                    parameters.append(constraint.maximum)
            return session.frame(
                "SELECT * FROM v_window_metrics WHERE "
                + " AND ".join(where)
                + " ORDER BY aa_start, aa_end, window_id",
                parameters,
            )
