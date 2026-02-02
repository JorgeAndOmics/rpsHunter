"""
    Executes RPS-BLAST searches for multiple species, processes results, and saves them to disk.

    This script runs `rpsblast` on species-specific FASTA files using a shared RPS database.
    It formats and filters results, then writes the merged output to CSV and Parquet formats.
"""

import argparse
import concurrent.futures
import logging
import os
import re
import subprocess
from typing import List, Optional, Tuple, Dict, Any

from Bio import SeqIO

import pandas as pd
from tqdm import tqdm

import defaults
from colored_logging import colored_logging


# -----------------------------------------------------------------------------
# RPS-BLAST Execution Function
# -----------------------------------------------------------------------------

def rpsblaster(
    command: str,
    input_database_path: str,
    query_file_path: str,
    species: str,
    evalue: float,
    asn_output_dir: str = None
) -> Tuple[Optional[str], Optional[str]]:
    """
    Runs an RPS-BLAST search for given sequences against a given database.

        Parameters
        ----------
            :param command: The RPS-BLAST command to execute (e.g., 'rpsblast').
            :param input_database_path: Path to the RPS-BLAST database.
            :param query_file_path: Path to the species-specific FASTA query file.
            :param species: Name of the species being processed.
            :param evalue: E-value threshold for filtering BLAST results.
            :param asn_output_dir: ASN output directory. Defaults to ASN_RPSBLAST_DIR.

        Returns
        -------
            :returns: A tuple containing (BLAST output as tabular string, ASN.1 file path).
                      If BLAST fails, returns (None, None).

        Raises
        ------
            :raises subprocess.CalledProcessError: If BLAST command fails.
            :raises Exception: For unexpected errors during execution.
    """
    try:
        if asn_output_dir is None:
            asn_output_dir = defaults.PATH_DICT['ASN_RPSBLAST_DIR']
        asn_file_name: str = os.path.join(asn_output_dir, f'{species}.asn')

        rpsblast_command: List[str] = [
            command,
            '-db', input_database_path,
            '-query', query_file_path,
            '-evalue', str(evalue),
            '-outfmt', '11',
            '-out', asn_file_name
        ]

        result = subprocess.run(rpsblast_command, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f'Error running {command} for {species}: {result.stderr}')
            return None, None

        blast_formatter_command: List[str] = [
            'blast_formatter',
            '-archive', asn_file_name,
            '-outfmt', '6 std stitle'
        ]

        formatter_result: subprocess.CompletedProcess = subprocess.run(
            blast_formatter_command, capture_output=True, text=True
        )

        if formatter_result.returncode != 0:
            logging.error(f'Error running blast_formatter for {species}: {formatter_result.stderr}')
            return None, None

        blast_output: str = formatter_result.stdout
        return blast_output, asn_file_name

    except Exception as e:
        logging.error(f'An exception occurred while running RPS-BLAST for {species}: {str(e)}')
        return None, None


# -----------------------------------------------------------------------------
# RPS-BLAST Output Parser
# -----------------------------------------------------------------------------

def parse_blast_output(blast_output: str) -> pd.DataFrame:
    """
    Parses the RPS-BLAST tabular output into a pandas DataFrame.

        Parameters
        ----------
            :param blast_output: String containing tab-delimited BLAST results.

        Returns
        -------
            :returns: DataFrame containing parsed BLAST result rows.

        Raises
        ------
            :raises ValueError: If the output format is unexpected.
    """
    columns: List[str] = [
        'Query ID', 'Tag', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches', 'Gap Openings',
        'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value', 'Bit Score', 'Subject Title'
    ]

    if not blast_output:
        return pd.DataFrame(columns=columns)

    rows: List[List[str]] = [line.split('\t') for line in blast_output.strip().split('\n')]
    df: pd.DataFrame = pd.DataFrame(rows, columns=columns)
    return df


# -----------------------------------------------------------------------------
# Per-Species RPS-BLAST Processing
# -----------------------------------------------------------------------------

def process_rps_species(species: str, fasta_input_dir: str = None, asn_output_dir: str = None) -> Optional[pd.DataFrame]:
    """
    Processes a single species FASTA file with RPS-BLAST and returns filtered results.

        Parameters
        ----------
            :param species: Species name whose sequences will be analyzed via RPS-BLAST.
            :param fasta_input_dir: Directory containing input FASTA files. Defaults to FASTA_OUTPUT_DIR.
            :param asn_output_dir: ASN output directory. Defaults to ASN_RPSBLAST_DIR.

        Returns
        -------
            :returns: DataFrame of filtered BLAST results or None if no results found.

        Raises
        ------
            :raises FileNotFoundError: If the FASTA file does not exist.
            :raises Exception: If any unexpected error occurs.
    """
    if fasta_input_dir is None:
        fasta_input_dir = defaults.PATH_DICT['FASTA_OUTPUT_DIR']
    if asn_output_dir is None:
        asn_output_dir = defaults.PATH_DICT['ASN_RPSBLAST_DIR']
    fasta_file_path: str = os.path.join(fasta_input_dir, f'{species}.fa')

    if not os.path.exists(fasta_file_path):
        logging.warning(f'FASTA file for species {species} does not exist at {fasta_file_path}.')
        return None

    blast_output, asn_file_name = rpsblaster(
        command='rpsblast',
        input_database_path=defaults.PATH_DICT['RPS_DB'],
        query_file_path=fasta_file_path,
        species=species,
        evalue=defaults.E_VALUE_THRESHOLD,
        asn_output_dir=asn_output_dir
    )

    if not blast_output:
        logging.warning(f'No RPS-BLAST output for species {species}.')
        return None

    blast_df: pd.DataFrame = parse_blast_output(blast_output)

    if blast_df.empty:
        logging.warning(f'No data returned from RPS-BLAST for species {species}.')
        return None

    # Convert numeric columns, coercing errors
    numeric_columns: List[str] = ['Pct Identity', 'E-value', 'Alignment Length', 'S. Start', 'S. End', 'Bit Score']
    for col in numeric_columns:
        blast_df[col] = pd.to_numeric(blast_df[col], errors='coerce')

    blast_df = blast_df.dropna(subset=numeric_columns)
    blast_df['Species'] = species

    blast_df['Tag'] = blast_df['Query ID'].str.extract(rf'\|tag:(\w{{{defaults.RANDOM_ID_LENGTH}}})', expand=False)
    blast_df['Query ID'] = blast_df['Query ID'].str.replace(rf'\|tag:\w{{{defaults.RANDOM_ID_LENGTH}}}', '', regex=True)

    return blast_df


# -----------------------------------------------------------------------------
# Main Workflow
# -----------------------------------------------------------------------------

def main(
    species_list: List[str],
    fasta_input_dir: str = None,
    output_csv: str = None,
    output_parquet: str = None,
    asn_output_dir: str = None
) -> None:
    """
    Runs RPS-BLAST across all species in the provided list and saves results.

        Parameters
        ----------
            :param species_list: List of species names to run RPS-BLAST against.
            :param fasta_input_dir: Directory containing input FASTA files. Defaults to FASTA_OUTPUT_DIR.
            :param output_csv: Output CSV path. Defaults to TABLE_OUTPUT_DIR/rpsblast.csv.
            :param output_parquet: Output Parquet path. Defaults to TABLE_OUTPUT_DIR/rpsblast.parquet.
            :param asn_output_dir: ASN output directory. Defaults to ASN_RPSBLAST_DIR.

        Returns
        -------
            :returns: None. Writes CSV and Parquet outputs to disk.

        Raises
        ------
            :raises IOError: If result files cannot be written to disk.
    """
    colored_logging(log_file_name='rpsblast.txt')

    # Set default output paths
    if output_csv is None:
        output_csv = os.path.join(defaults.PATH_DICT['TABLE_OUTPUT_DIR'], 'rpsblast.csv')
    if output_parquet is None:
        output_parquet = os.path.join(defaults.PATH_DICT['TABLE_OUTPUT_DIR'], 'rpsblast.parquet')
    if asn_output_dir is None:
        asn_output_dir = defaults.PATH_DICT['ASN_RPSBLAST_DIR']

    # Ensure ASN output directory exists
    os.makedirs(asn_output_dir, exist_ok=True)

    dfs: List[pd.DataFrame] = []

    with tqdm(total=len(species_list), desc='Running RPSBLAST against species...') as pbar:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=defaults.MAX_THREADPOOL_WORKERS
        ) as executor:
            future_to_species: Dict[concurrent.futures.Future, str] = {
                executor.submit(process_rps_species, species, fasta_input_dir, asn_output_dir): species
                for species in species_list
            }

            for future in concurrent.futures.as_completed(future_to_species):
                species: str = future_to_species[future]
                try:
                    blast_df: Optional[pd.DataFrame] = future.result()
                    if blast_df is not None and not blast_df.empty:
                        dfs.append(blast_df)
                except Exception as exc:
                    logging.error(f'Error processing species {species}: {exc}')
                finally:
                    pbar.update(1)

    if dfs:
        df: pd.DataFrame = pd.concat(dfs, axis=0, ignore_index=True)

        # Remove duplicates
        df = df.drop_duplicates()

        # Save outputs
        df.to_csv(output_csv, index=False)
        logging.info(f'Saved DataFrame as CSV to {output_csv}')

        df.to_parquet(output_parquet, index=False)
        logging.info(f'Saved DataFrame as Parquet to {output_parquet}')
    else:
        logging.warning('No data frames to concatenate. No output files were created.')


# -----------------------------------------------------------------------------
# Entry Point
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description='Run RPS-BLAST for a list of species.'
    )

    parser.add_argument(
        '--species_list',
        type=str,
        nargs='+',
        required=False,
        default=None,
        help='List of species to process. Defaults to species from config.'
    )

    parser.add_argument(
        '--input_dir',
        type=str,
        default=None,
        help='Directory containing input FASTA files. Defaults to FASTA_OUTPUT_DIR.'
    )

    parser.add_argument(
        '--output_csv',
        type=str,
        default=None,
        help='Output CSV path. Defaults to TABLE_OUTPUT_DIR/rpsblast.csv'
    )

    parser.add_argument(
        '--output_parquet',
        type=str,
        default=None,
        help='Output Parquet path. Defaults to TABLE_OUTPUT_DIR/rpsblast.parquet'
    )

    parser.add_argument(
        '--asn_output_dir',
        type=str,
        default=None,
        help='ASN output directory. Defaults to ASN_RPSBLAST_DIR'
    )

    args: argparse.Namespace = parser.parse_args()
    species_list: List[str] = args.species_list if args.species_list else list(defaults.SPECIES)

    main(
        species_list=species_list,
        fasta_input_dir=args.input_dir,
        output_csv=args.output_csv,
        output_parquet=args.output_parquet,
        asn_output_dir=args.asn_output_dir
    )
