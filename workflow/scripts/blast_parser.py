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

import os
import argparse
import logging
from typing import List, Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

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


def species_divider(df: pd.DataFrame) -> List[pd.DataFrame]:
    """
    Divides the filtered DataFrame into a list of species-specific DataFrames.

        Parameters
        ----------
            :param df: A filtered DataFrame that includes a 'Species' column.

        Returns
        -------
            :returns: A list of DataFrames, each corresponding to a unique species in the input data.

        Raises
        ------
            :raises KeyError: If the 'Species' column is not present in the DataFrame.
    """
    species_list = df['Species'].unique()
    species_dfs: List[pd.DataFrame] = []

    for species in species_list:
        species_df = df[df['Species'] == species].copy()
        species_dfs.append(species_df)

    return species_dfs


def fasta_generator(df: pd.DataFrame) -> Optional[None]:
    """
    Generates and writes a FASTA file for the provided species-specific DataFrame.

        Parameters
        ----------
            :param df: A DataFrame containing BLAST sequences for a single species.
                       Must contain 'Subject Sequence', 'Subject ID', 'S. Start', 'S. End', and 'Species' columns.

        Returns
        -------
            :returns: None. Writes a FASTA file to disk for the species represented in the DataFrame.

        Raises
        ------
            :raises KeyError: If expected columns are missing from the DataFrame.
            :raises IOError: If the FASTA file cannot be written to disk.
    """
    if 'Subject Sequence' not in df.columns:
        return

    species: Optional[str] = None if df.empty else df['Species'].unique()[0]

    df['Subject Sequence'] = (
        df['Subject Sequence']
        .str.replace('-', '', regex=False)
        .str.replace('*', '', regex=False)
    )

    seq_records: List[SeqRecord] = []
    for _, row in df.iterrows():
        header: str = f"{row['Subject ID']}:{row['S. Start']}-{row['S. End']}|tag:{row['Tag']}"
        sequence: str = str(row['Subject Sequence'])
        seq_records.append(SeqRecord(Seq(sequence), id=header, description=''))

    if seq_records:
        output_fasta_path: str = os.path.join(
            defaults.PATH_DICT['FASTA_OUTPUT_DIR'],
            f'{species}.fa'
        )
        SeqIO.write(seq_records, output_fasta_path, 'fasta')


# -------------------------------------------------------------------------
# Main Execution
# -------------------------------------------------------------------------

if __name__ == '__main__':
    colored_logging(log_file_name='blast_parser.txt')

    parser = argparse.ArgumentParser(description='Parse and filter BLAST results.')
    parser.add_argument(
        '--input_parquet_file',
        type=str,
        required=True,
        help='Path to the input Parquet file containing BLAST results.'
    )
    parser.add_argument(
        '--export_fasta',
        action='store_true',
        help='Export filtered sequences to FASTA files.'
    )

    args = parser.parse_args()

    input_parquet_path: str = os.path.join(
        defaults.PATH_DICT['TABLE_OUTPUT_DIR'],
        args.input_parquet_file
    )
    blast_df: pd.DataFrame = pd.read_parquet(input_parquet_path)

    # Apply quality thresholds
    blast_df = table_filter(blast_df)

    # Apply ORF filtering if ORF analysis was run (unified pipeline)
    blast_df = orf_aware_filter(blast_df)

    if args.export_fasta:
        species_dfs: List[pd.DataFrame] = species_divider(blast_df)

        for df in species_dfs:
            fasta_generator(df)

        logging.info(f"Exported FASTA files for {len(species_dfs)} species")