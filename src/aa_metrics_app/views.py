"""Trusted schema-2 analytical views shared by supported SQL backends."""

from __future__ import annotations

SCHEMA_V2_VIEW_SQL: tuple[str, ...] = (
    """
    CREATE VIEW v_block_metrics AS
    SELECT
        s.stratum_type,
        s.stratum_value,
        s.display_label,
        s.small_sample_warning,
        b.block_number, b.block_name, b.aa_start, b.aa_end, b.aa_length,
        seq1.sequence_text AS top_exact_sequence,
        seq2.sequence_text AS allelic_consensus_sequence,
        seq3.sequence_text AS all_call_consensus_sequence,
        seq4.sequence_text AS modal_all_call_sequence,
        m.*
    FROM block_metric AS m
    JOIN stratum AS s ON s.stratum_id = m.stratum_id
    JOIN block AS b ON b.block_id = m.block_id
    LEFT JOIN sequence AS seq1 ON seq1.sequence_id = m.top_exact_sequence_id
    LEFT JOIN sequence AS seq2 ON seq2.sequence_id = m.allelic_consensus_sequence_id
    LEFT JOIN sequence AS seq3 ON seq3.sequence_id = m.all_call_consensus_sequence_id
    LEFT JOIN sequence AS seq4 ON seq4.sequence_id = m.modal_all_call_sequence_id
    """.strip(),
    """
    CREATE VIEW v_locus_metrics AS
    SELECT
        s.stratum_type,
        s.stratum_value,
        s.display_label,
        s.small_sample_warning,
        b.block_number, b.block_name, l.aa_position, l.block_offset,
        l.reference_residue,
        m.*
    FROM locus_metric AS m
    JOIN stratum AS s ON s.stratum_id = m.stratum_id
    JOIN locus AS l ON l.locus_id = m.locus_id
    JOIN block AS b ON b.block_id = l.block_id
    """.strip(),
    """
    CREATE VIEW v_window_metrics AS
    SELECT
        s.stratum_type,
        s.stratum_value,
        s.display_label,
        s.small_sample_warning,
        b.block_number, b.block_name, w.window_size, w.aa_start, w.aa_end,
        w.block_offset_start, w.block_offset_end,
        ref.sequence_text AS reference_sequence,
        seq1.sequence_text AS allelic_consensus_sequence,
        seq2.sequence_text AS all_call_consensus_sequence,
        seq3.sequence_text AS top_exact_peptide_sequence,
        seq4.sequence_text AS modal_all_call_sequence,
        m.*
    FROM window_metric AS m
    JOIN stratum AS s ON s.stratum_id = m.stratum_id
    JOIN "window" AS w ON w.window_id = m.window_id
    JOIN block AS b ON b.block_id = w.block_id
    JOIN sequence AS ref ON ref.sequence_id = w.reference_sequence_id
    LEFT JOIN sequence AS seq1 ON seq1.sequence_id = m.allelic_consensus_sequence_id
    LEFT JOIN sequence AS seq2 ON seq2.sequence_id = m.all_call_consensus_sequence_id
    LEFT JOIN sequence AS seq3 ON seq3.sequence_id = m.top_exact_peptide_sequence_id
    LEFT JOIN sequence AS seq4 ON seq4.sequence_id = m.modal_all_call_sequence_id
    """.strip(),
)
