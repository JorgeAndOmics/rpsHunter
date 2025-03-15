import logging
import os
import subprocess
from typing import List, Optional, Tuple
import concurrent.futures

import pandas as pd
from tqdm import tqdm

import defaults
from colored_logging import colored_logging


def rpsblaster(
    command: str,
    input_database_path: str,
    query_file_path: str,
    species: str,
    evalue: float,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Runs an RPS-BLAST search for given sequences against a given database.

    Parameters
    ----------
    command : str
        The command to run RPS-BLAST.
    input_database_path : str
        The path to the RPS database.
    query_file_path : str
        The path to the query file (species FASTA file).
    species : str
        The species name, used to name the ASN.1 file.
    evalue : float
        The e-value threshold for RPS-BLAST.

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        A tuple containing the tabular RPS-BLAST output and the path to the ASN.1 file.
    """
    try:
        # Create ASN.1 file path
        asn_file_name = os.path.join(defaults.ASN_RPSBLAST_DIR, f"{species}.asn")

        # Construct the RPS-BLAST command
        rpsblast_command = [
            command,
            '-db',
            input_database_path,
            '-query',
            query_file_path,
            '-evalue',
            str(evalue),
            '-outfmt',
            '11',
            '-out',
            asn_file_name
        ]

        # Run the RPS-BLAST command
        result = subprocess.run(rpsblast_command, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f'Error running {command} for {species}: {result.stderr}')
            return None, None

        # Run blast_formatter to get tabular output
        blast_formatter_command = [
            'blast_formatter',
            '-archive',
            asn_file_name,
            '-outfmt',
            '6 std stitle',
        ]

        formatter_result = subprocess.run(blast_formatter_command, capture_output=True, text=True)
        if formatter_result.returncode != 0:
            logging.error(f'Error running blast_formatter for {species}: {formatter_result.stderr}')
            return None, None

        blast_output = formatter_result.stdout

        return blast_output, asn_file_name

    except Exception as e:
        logging.error(f'An exception occurred while running RPS-BLAST for {species}: {str(e)}')
        return None, None


def parse_blast_output(blast_output: str) -> pd.DataFrame:
    """Parses the BLAST output into a DataFrame."""
    columns = ['Query ID', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches', 'Gap Openings',
               'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value', 'Bit Score', 'Subject Title']
    if not blast_output:
        return pd.DataFrame(columns=columns)

    rows = [line.split('\t') for line in blast_output.strip().split('\n')]
    df = pd.DataFrame(rows, columns=columns)
    return df


def process_rps_species(species: str) -> Optional[pd.DataFrame]:
    """Processes RPS-BLAST results for a single species."""
    fasta_file_path = os.path.join(defaults.FASTA_OUTPUT_DIR, f'{species}.fa')

    if not os.path.exists(fasta_file_path):
        logging.warning(f'FASTA file for species {species} does not exist at {fasta_file_path}.')
        return None

    # Run rpsblaster
    blast_output, asn_file_name = rpsblaster(
        command='rpsblast',
        input_database_path=defaults.RPS_DB,
        query_file_path=fasta_file_path,
        species=species,
        evalue=defaults.E_VALUE_THRESHOLD
    )

    if not blast_output:
        logging.warning(f'No RPS-BLAST output for species {species}.')
        return None

    # Parse BLAST output
    blast_df = parse_blast_output(blast_output)
    if blast_df.empty:
        logging.warning(f'No data returned from RPS-BLAST for species {species}.')
        return None

    # Convert data types using pd.to_numeric with errors='coerce'
    numeric_columns = ['Pct Identity', 'E-value', 'Alignment Length', 'S. Start', 'S. End', 'Bit Score']
    for col in numeric_columns:
        blast_df[col] = pd.to_numeric(blast_df[col], errors='coerce')

    # Drop rows with NaN values in the numeric columns
    blast_df.dropna(subset=numeric_columns, inplace=True)

    # Apply filters and create a copy to avoid SettingWithCopyWarning
    filtered_df = blast_df[
        (blast_df['Pct Identity'] >= defaults.PERC_IDENTITY_THRESHOLD) &
        (blast_df['E-value'] <= defaults.E_VALUE_THRESHOLD) &
        (blast_df['Alignment Length'] >= defaults.SEQ_LENGTH_THRESHOLD) &
        (blast_df['Bit Score'] >= defaults.BITSCORE_THRESHOLD)
    ].copy()

    if filtered_df.empty:
        logging.warning(f'No hits after filtering for species {species}.')
        return None

    filtered_df['Species'] = species

    return filtered_df


def main():
    # Set up logging
    colored_logging(log_file_name='rpsblast.txt')

    # Initialize list to collect DataFrames
    dfs = []

    # Use ThreadPoolExecutor to process species concurrently
    with tqdm(total=len(defaults.SPECIES), desc='Running RPSBLAST against species...') as pbar:
        with concurrent.futures.ThreadPoolExecutor(max_workers=defaults.MAX_THREADPOOL_WORKERS) as executor:
            future_to_species = {
                executor.submit(process_rps_species, species): species
                for species in defaults.SPECIES
            }
            for future in concurrent.futures.as_completed(future_to_species):
                species = future_to_species[future]
                try:
                    filtered_df = future.result()
                    if filtered_df is not None and not filtered_df.empty:
                        dfs.append(filtered_df)
                except Exception as exc:
                    logging.error(f'Error processing species {species}: {exc}')
                finally:
                    pbar.update(1)

    # Concatenate all DataFrames
    if dfs:
        df = pd.concat(dfs, axis=0, ignore_index=True)
        # Save the DataFrame to PARQUET and CSV file
        output_csv_path = os.path.join(defaults.TABLE_OUTPUT_DIR, 'rpsblast.csv')
        output_parquet_path = os.path.join(defaults.TABLE_OUTPUT_DIR, 'rpsblast.parquet')

        df.to_csv(output_csv_path, index=False)
        logging.info(f'Saved DataFrame as CSV to {output_csv_path}')

        df.to_parquet(output_parquet_path, index=False)
        logging.info(f'Saved DataFrame as Parquet to {output_parquet_path}')
    else:
        logging.warning('No data frames to concatenate. No output files were created.')


if __name__ == '__main__':
    main()
