"""
Module: blast_parser.py

Description
-----------
    Parses, filters, and exports BLAST result sequences as per-species FASTA files.
    Applies a strict, auditable filter chain:

        table_filter  →  orf_aware_filter  →  hmm_aware_filter  →  FASTA export

    Each gate follows the same contract:
        - Column missing        → step was disabled → pass through unchanged.
        - Column present, ≥1 True  → keep only True rows.
        - Column present, 0 True   → return empty DF (species legitimately produced 0 hits).

    An audit parquet is written to fastas/{species}.parquet before any filtering.
    It contains every input row plus four boolean columns:
        Quality_Pass  – passed table_filter thresholds
        ORF_Pass      – True/False/NaN (NaN = ORF step not run)
        HMM_Pass      – True/False/NaN (NaN = HMMER step not run)
        Selected      – made it through all gates to FASTA export

Requirements
------------
    - pandas
    - biopython
    - defaults (user-defined thresholds and paths)
    - colored_logging (custom logging utility)
"""

# -------------------------------------------------------------------------
# Imports
# -------------------------------------------------------------------------

import argparse
import logging
from pathlib import Path
from typing import List

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import defaults
from colored_logging import colored_logging


# -------------------------------------------------------------------------
# Functions
# -------------------------------------------------------------------------


def table_filter(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filters the table based on identity, e-value, alignment length, and bit score thresholds.

        Parameters
        ----------
            :param df: The input DataFrame containing BLAST results to be filtered.

        Returns
        -------
            :returns: A DataFrame containing only the rows that pass the filtering thresholds.

        Raises
        ------
            :raises KeyError: If required columns are missing in the input DataFrame.
    """
    return df[
        (df['Pct Identity'] >= defaults.PERC_IDENTITY_THRESHOLD) &
        (df['E-value'] <= defaults.E_VALUE_THRESHOLD) &
        (df['Alignment Length'] >= defaults.SEQ_LENGTH_THRESHOLD) &
        (df['Bit Score'] >= defaults.BITSCORE_THRESHOLD)
    ].copy()


def orf_aware_filter(df: pd.DataFrame) -> pd.DataFrame:
    """
    Strict ORF gate: keeps only rows where ORF_Filtered is True.

    Contract:
        - Column missing        → ORF step was disabled → pass through unchanged.
        - Column present, ≥1 True  → keep only True rows.
        - Column present, 0 True   → return empty DF (real result: ORF ran, nothing passed).

        Parameters
        ----------
            :param df: The input DataFrame containing BLAST results.

        Returns
        -------
            :returns: Filtered DataFrame (may be empty).
    """
    if 'ORF_Filtered' not in df.columns:
        logging.info("No ORF_Filtered column — ORF step was not run; passing through")
        return df

    filtered_df = df[df['ORF_Filtered'] == True].copy()
    logging.info(f"ORF gate: {len(filtered_df)}/{len(df)} rows passed")
    return filtered_df


def hmm_aware_filter(df: pd.DataFrame) -> pd.DataFrame:
    """
    Strict HMM gate: keeps only rows where HMM_Filtered is True.

    Same contract as orf_aware_filter but on the HMM_Filtered column.

        Parameters
        ----------
            :param df: The input DataFrame containing BLAST results.

        Returns
        -------
            :returns: Filtered DataFrame (may be empty).
    """
    if 'HMM_Filtered' not in df.columns:
        logging.info("No HMM_Filtered column — HMMER step was not run; passing through")
        return df

    filtered_df = df[df['HMM_Filtered'] == True].copy()
    logging.info(f"HMM gate: {len(filtered_df)}/{len(df)} rows passed")
    return filtered_df


# -------------------------------------------------------------------------
# Main Execution
# -------------------------------------------------------------------------

if __name__ == '__main__':
    colored_logging(log_file_name='blast_parser.txt')

    parser = argparse.ArgumentParser(description='Parse and filter BLAST results for a single species.')
    parser.add_argument('--species', type=str, required=True,
                        help='Species name to process.')
    parser.add_argument('--input-parquet', type=str, required=True,
                        help='Path to the per-species enriched Parquet file.')
    parser.add_argument('--output-fasta', type=str, required=True,
                        help='Path to the output FASTA file.')
    parser.add_argument('--output-audit', type=str, required=True,
                        help='Path to the per-species audit Parquet file.')

    args = parser.parse_args()

    # ── Read full input ─────────────────────────────────────────────────────
    raw_df: pd.DataFrame = pd.read_parquet(args.input_parquet)

    # ── Compute audit flags on the full input (before any row is dropped) ───
    quality_mask = (
        (raw_df['Pct Identity']      >= defaults.PERC_IDENTITY_THRESHOLD) &
        (raw_df['E-value']           <= defaults.E_VALUE_THRESHOLD) &
        (raw_df['Alignment Length']  >= defaults.SEQ_LENGTH_THRESHOLD) &
        (raw_df['Bit Score']         >= defaults.BITSCORE_THRESHOLD)
    )
    raw_df['Quality_Pass'] = quality_mask

    # ORF_Pass / HMM_Pass: True/False if the column exists, NaN if the step was not run
    if 'ORF_Filtered' in raw_df.columns:
        raw_df['ORF_Pass'] = raw_df['ORF_Filtered'].astype(bool)
    else:
        raw_df['ORF_Pass'] = pd.array([pd.NA] * len(raw_df), dtype=pd.BooleanDtype())

    if 'HMM_Filtered' in raw_df.columns:
        raw_df['HMM_Pass'] = raw_df['HMM_Filtered'].astype(bool)
    else:
        raw_df['HMM_Pass'] = pd.array([pd.NA] * len(raw_df), dtype=pd.BooleanDtype())

    # ── Run the filter chain ────────────────────────────────────────────────
    blast_df = table_filter(raw_df)
    blast_df = orf_aware_filter(blast_df)
    blast_df = hmm_aware_filter(blast_df)

    # ── Mark Selected on the full audit frame ──────────────────────────────
    selected_idx = blast_df.index
    raw_df['Selected'] = raw_df.index.isin(selected_idx)

    # ── Write audit parquet ─────────────────────────────────────────────────
    output_audit_path: Path = Path(args.output_audit)
    output_audit_path.parent.mkdir(parents=True, exist_ok=True)
    raw_df.to_parquet(output_audit_path, index=False)
    logging.info(f"Audit parquet written: {len(raw_df)} rows, {raw_df['Selected'].sum()} selected → {output_audit_path}")

    # ── Write FASTA ─────────────────────────────────────────────────────────
    output_fasta_path: Path = Path(args.output_fasta)
    output_fasta_path.parent.mkdir(parents=True, exist_ok=True)

    if 'Subject Sequence' in blast_df.columns and not blast_df.empty:
        blast_df['Subject Sequence'] = (
            blast_df['Subject Sequence']
            .str.replace('-', '', regex=False)
            .str.replace('*', '', regex=False)
        )

        seq_records: List[SeqRecord] = []
        for _, row in blast_df.iterrows():
            header: str = f"{row['Subject ID']}:{row['S. Start']}-{row['S. End']}|tag:{row['Tag']}"
            seq_records.append(SeqRecord(Seq(str(row['Subject Sequence'])), id=header, description=''))

        SeqIO.write(seq_records, output_fasta_path, 'fasta')
        logging.info(f"Exported {len(seq_records)} sequences to {output_fasta_path}")
    else:
        # Touch an empty FASTA so Snakemake output is satisfied
        output_fasta_path.touch()
        logging.warning(f"No sequences passed all gates for {args.species} — empty FASTA written")
