import logging
import os
import subprocess
from typing import List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

import defaults
from colored_logging import colored_logging


def blaster(
    command: str,
    input_database_path: str,
    query_file_path: str,
    subject: str,
    evalue: float,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Runs a BLAST search for a given object against a given database.

    Parameters
    ----------
    command : str
        The command to run BLAST.
    input_database_path : str
        The path to the database.
    query_file_path : str
        The path to the query file.
    subject : str
        The particular genome against whose database it's being BLASTed.
    evalue : float
        The e-value threshold for BLAST.

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        A tuple containing the tabular BLAST output and the path to the ASN.1 file.
    """
    input_path = os.path.join(input_database_path, subject, subject)
    try:
        # Ensure the ASN output directory exists
        os.makedirs(defaults.ASN_OUTPUT_DIR, exist_ok=True)
        # Create ASN.1 file path
        asn_file_name = os.path.join(defaults.ASN_OUTPUT_DIR, f"{subject}.asn")

        # Construct the BLAST command
        blast_command = [
            command,
            '-db',
            input_path,
            '-query',
            query_file_path,
            '-evalue',
            str(evalue),
            '-outfmt',
            '11',
            '-out',
            asn_file_name
        ]

        # Run the BLAST command
        result = subprocess.run(blast_command, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f'Error running {command} for {subject}: {result.stderr}')
            return None, None

        # Run blast_formatter to get tabular output
        blast_formatter_command = [
            'blast_formatter',
            '-archive',
            asn_file_name,
            '-outfmt',
            '6',
        ]

        formatter_result = subprocess.run(blast_formatter_command, capture_output=True, text=True)
        if formatter_result.returncode != 0:
            logging.error(f'Error running blast_formatter for {subject}: {formatter_result.stderr}')
            return None, None

        blast_output = formatter_result.stdout

        return blast_output, asn_file_name

    except Exception as e:
        logging.error(f'An exception occurred while running BLAST for {subject}: {str(e)}')
        return None, None


def parse_blast_output(blast_output: str) -> pd.DataFrame:
    """Parses the BLAST output into a DataFrame."""
    columns = ['Query ID', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches', 'Gap Openings',
               'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value', 'Bit Score']
    if not blast_output:
        return pd.DataFrame(columns=columns)

    rows = [line.split('\t') for line in blast_output.strip().split('\n')]
    df = pd.DataFrame(rows, columns=columns)
    return df


def parse_sequence_output(sequence_output: str) -> pd.DataFrame:
    """Parses the sequence output from blast_formatter into a DataFrame."""
    columns = ['Subject ID', 'Hit Description', 'S. Start', 'S. End', 'Subject Sequence']
    if not sequence_output:
        return pd.DataFrame(columns=columns)

    rows = [line.split('\t') for line in sequence_output.strip().split('\n') if line]
    df = pd.DataFrame(rows, columns=columns)
    df['Subject Sequence'] = df['Subject Sequence'].str.replace('-', '').str.replace('*', '')
    df['S. Start'] = df['S. Start'].astype(int)
    df['S. End'] = df['S. End'].astype(int)
    return df


def process_species(species: str) -> Tuple[Optional[pd.DataFrame], List[SeqRecord]]:
    """Processes BLAST results for a single species."""
    # logging.info(f'Processing species: {species}')
    blast_output, asn_file_name = blaster(
        command='tblastn',
        input_database_path=defaults.SPECIES_DB,
        query_file_path=defaults.QUERY_SEQ,
        subject=species,
        evalue=defaults.E_VALUE_THRESHOLD
    )

    if not blast_output:
        logging.warning(f'No BLAST output for species {species}.')
        return None, []

    # Parse BLAST output
    blast_df = parse_blast_output(blast_output)
    if blast_df.empty:
        logging.warning(f'No data returned from BLAST for species {species}.')
        return None, []

    # Convert data types
    blast_df = blast_df.astype({
        'Pct Identity': float,
        'E-value': float,
        'Alignment Length': int,
        'S. Start': int,
        'S. End': int
    })

    # Apply filters and create a copy to avoid SettingWithCopyWarning
    filtered_df = blast_df[
        (blast_df['Pct Identity'] >= defaults.PERC_IDENTITY_THRESHOLD) &
        (blast_df['E-value'] <= defaults.E_VALUE_THRESHOLD) &
        (blast_df['Alignment Length'] >= defaults.SEQ_LENGTH_THRESHOLD)
    ].copy()

    if filtered_df.empty:
        logging.warning(f'No hits after filtering for species {species}.')
        return None, []

    filtered_df['Species'] = species

    # Run blast_formatter to get sequences
    sequence_formatter_command = [
        'blast_formatter',
        '-archive',
        asn_file_name,
        '-outfmt',
        '6 sseqid salltitles sstart send sseq',
    ]

    sequence_formatter_result = subprocess.run(sequence_formatter_command, capture_output=True, text=True)
    if sequence_formatter_result.returncode != 0:
        logging.error(f'Error running blast_formatter for sequences of {species}: {sequence_formatter_result.stderr}')
        return filtered_df, []

    sequence_df = parse_sequence_output(sequence_formatter_result.stdout)

    # Merge DataFrames
    merged_df = pd.merge(
        filtered_df,
        sequence_df,
        on=['Subject ID', 'S. Start', 'S. End'],
        how='left'
    )

    # Collect sequences
    seq_records = []
    for _, row in merged_df.iterrows():
        header = f"{row['Hit Description']} {row['Subject ID']}:{row['S. Start']}-{row['S. End']}"
        sequence = row['Subject Sequence']
        seq_records.append(SeqRecord(Seq(sequence), id=header, description=''))

    return filtered_df, seq_records


def main():
    # Set up logging
    colored_logging(log_file_name='blast.txt')

    # Initialize list to collect DataFrames
    dfs = []

    # Ensure output directories exist
    os.makedirs(defaults.FASTA_OUTPUT_DIR, exist_ok=True)
    os.makedirs(defaults.TABLE_OUTPUT_DIR, exist_ok=True)

    with tqdm(total=len(defaults.SPECIES)) as pbar:
        for species in defaults.SPECIES:
            pbar.set_description(f'Processing {species.replace("_", " ")}...')
            filtered_df, seq_records = process_species(species)
            if filtered_df is not None and not filtered_df.empty:
                dfs.append(filtered_df)
                # Save sequences to FASTA file for this species
                if seq_records:
                    output_fasta_path = os.path.join(defaults.FASTA_OUTPUT_DIR, f'{species}_sequences.fasta')
                    SeqIO.write(seq_records, output_fasta_path, 'fasta')
                    logging.info(f'Saved sequences to {output_fasta_path}')
            pbar.update(1)

    # Concatenate all DataFrames
    if dfs:
        df = pd.concat(dfs, axis=0, ignore_index=True)
        # Save the DataFrame to PARQUET and CSV file
        output_csv_path = os.path.join(defaults.TABLE_OUTPUT_DIR, 'blast.csv')
        output_parquet_path = os.path.join(defaults.TABLE_OUTPUT_DIR, 'blast.parquet')

        df.to_csv(output_csv_path, index=False)
        logging.info(f'Saved DataFrame as CSV to {output_csv_path}')

        df.to_parquet(output_parquet_path, index=False)
        logging.info(f'Saved DataFrame as Parquet to {output_parquet_path}')
    else:
        logging.warning('No data frames to concatenate. No output files were created.')


if __name__ == '__main__':
    main()
