"""
Module: blast_parser.py

Description
-----------
    Parses, filters, and optionally exports BLAST result sequences as per-species FASTA files.
    Applies filtering based on configurable identity, e-value, alignment length, and bitscore thresholds.

    ORF-aware behavior (unified pipeline):
    - If blast.parquet contains an 'ORF_Filtered' column with True values,
      only those sequences are exported to FASTA (ORF-filtered sequences).
    - If the column doesn't exist or has no True values, all sequences are used
      (backward compatible behavior).

Requirements
------------
    - pandas
    - biopython
    - tqdm
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
    Filters the table to ORF-filtered sequences if ORF analysis was run.

    This enables the unified pipeline where ORF analysis enriches blast.parquet
    with an 'ORF_Filtered' column. When this column exists and has True values,
    only those sequences are used for downstream processing.

        Parameters
        ----------
            :param df: The input DataFrame containing BLAST results.

        Returns
        -------
            :returns: A DataFrame filtered to ORF-passing sequences if ORF analysis
                      was run, otherwise the original DataFrame unchanged.
    """
    if 'ORF_Filtered' not in df.columns:
        logging.info("No ORF_Filtered column found - using all sequences (no ORF analysis run)")
        return df

    # Check if there are any ORF-filtered sequences
    orf_filtered_count = df['ORF_Filtered'].sum() if df['ORF_Filtered'].dtype == bool else (df['ORF_Filtered'] == True).sum()

    if orf_filtered_count > 0:
        filtered_df = df[df['ORF_Filtered'] == True].copy()
        logging.info(f"ORF analysis detected - using {len(filtered_df)} ORF-filtered sequences "
                     f"(from {len(df)} total rows)")
        return filtered_df
    else:
        logging.warning("ORF_Filtered column exists but no TRUE values found - using all sequences")
        return df


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

    args = parser.parse_args()

    blast_df: pd.DataFrame = pd.read_parquet(args.input_parquet)

    # Apply quality thresholds
    blast_df = table_filter(blast_df)

    # Apply ORF filtering if ORF analysis was run (unified pipeline)
    blast_df = orf_aware_filter(blast_df)

    # Write FASTA for this species
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

        if seq_records:
            SeqIO.write(seq_records, output_fasta_path, 'fasta')
            logging.info(f"Exported {len(seq_records)} sequences to {output_fasta_path}")
    else:
        logging.warning(f"No sequences to export for {args.species}")
