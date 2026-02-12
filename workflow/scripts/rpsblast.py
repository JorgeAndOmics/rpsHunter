"""
    Executes RPS-BLAST searches for multiple species, processes results, and saves them to disk.

    This script runs `rpsblast` on species-specific FASTA files using a shared RPS database.
    It formats and filters results, then writes the merged output to CSV and Parquet formats.
"""

import argparse
import logging
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

import defaults
from colored_logging import colored_logging


# -----------------------------------------------------------------------------
# RPS-BLAST Execution Function
# -----------------------------------------------------------------------------

def rpsblaster(
    command: str,
    input_database_path: Path,
    query_file_path: Path,
    species: str,
    evalue: float,
    asn_output_dir: Path = None
) -> Tuple[Optional[str], Optional[Path]]:
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
        asn_file_path: Path = asn_output_dir / f'{species}.asn'

        rpsblast_command: List[str] = [
            command,
            '-db', str(input_database_path),
            '-query', str(query_file_path),
            '-evalue', str(evalue),
            '-comp_based_stats', str(defaults.RPSBLAST_COMP_BASED_STATS),
            '-seg', defaults.RPSBLAST_SEG,
            '-outfmt', '11',
            '-out', str(asn_file_path)
        ]

        if defaults.RPSBLAST_WINDOW_SIZE and defaults.RPSBLAST_WINDOW_SIZE > 0:
            rpsblast_command[6:6] = ['-window_size', str(defaults.RPSBLAST_WINDOW_SIZE)]

        result = subprocess.run(rpsblast_command, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f'Error running {command} for {species}: {result.stderr}')
            return None, None

        blast_formatter_command: List[str] = [
            defaults.BLAST_FORMATTER_CMD,
            '-archive', str(asn_file_path),
            '-outfmt', '6 std stitle'
        ]

        formatter_result: subprocess.CompletedProcess = subprocess.run(
            blast_formatter_command, capture_output=True, text=True
        )

        if formatter_result.returncode != 0:
            logging.error(f'Error running blast_formatter for {species}: {formatter_result.stderr}')
            return None, None

        blast_output: str = formatter_result.stdout
        return blast_output, asn_file_path

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
        'Query ID', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches', 'Gap Openings',
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

def process_rps_species(species: str, fasta_input_dir: Path = None, asn_output_dir: Path = None, db_override: Path = None) -> Optional[pd.DataFrame]:
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

    fasta_file_path: Path = fasta_input_dir / f'{species}.fa'

    if not fasta_file_path.exists():
        logging.warning(f'FASTA file for species {species} does not exist at {fasta_file_path}.')
        return None

    rps_db = db_override if db_override else defaults.PATH_DICT['RPS_DB']
    blast_output, asn_file_name = rpsblaster(
        command=defaults.RPSBLAST_CMD,
        input_database_path=rps_db,
        query_file_path=fasta_file_path,
        species=species,
        evalue=defaults.RPSBLAST_E_VALUE,
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
    blast_df['Species_Name'] = defaults.SPECIES_DICT.get(species, species)

    blast_df['Tag'] = blast_df['Query ID'].str.extract(rf'\|tag:(\w{{{defaults.RANDOM_ID_LENGTH}}})', expand=False)
    blast_df['Query ID'] = blast_df['Query ID'].str.replace(rf'\|tag:\w{{{defaults.RANDOM_ID_LENGTH}}}', '', regex=True)

    return blast_df


# -----------------------------------------------------------------------------
# Main Workflow
# -----------------------------------------------------------------------------

def main() -> None:
    """
    Executes RPS-BLAST for a single species and writes per-species output.
    """
    colored_logging(log_file_name='rpsblast.txt')

    parser = argparse.ArgumentParser(description='Run RPS-BLAST for a single species.')
    parser.add_argument('--species', type=str, required=True,
                        help='Species name to process.')
    parser.add_argument('--db', type=str, default=None,
                        help='Path to RPS-BLAST database (default: full CDD).')
    args = parser.parse_args()

    defaults.PATH_DICT['ASN_RPSBLAST_DIR'].mkdir(parents=True, exist_ok=True)

    db_path = Path(args.db) if args.db else None
    result_df = process_rps_species(args.species, db_override=db_path)

    output_dir = defaults.PATH_DICT['RPSBLAST_SPECIES_DIR']
    output_dir.mkdir(parents=True, exist_ok=True)
    parquet_path = output_dir / f'{args.species}.parquet'
    asn_path = defaults.PATH_DICT['ASN_RPSBLAST_DIR'] / f'{args.species}.asn'

    if result_df is not None and not result_df.empty:
        result_df.to_parquet(parquet_path, index=False)
        logging.info(f'Saved {len(result_df)} rows for {args.species}')
    else:
        # Write an empty parquet with the correct schema so downstream rules don't crash.
        empty_columns = [
            'Query ID', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches',
            'Gap Openings', 'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value',
            'Bit Score', 'Subject Title', 'Species', 'Species_Name', 'Tag'
        ]
        pd.DataFrame(columns=empty_columns).to_parquet(parquet_path, index=False)
        # Touch the ASN file so rpsbproc_species input is satisfied.
        asn_path.touch()
        logging.warning(f'No RPS-BLAST results for {args.species} — empty outputs written.')


# -----------------------------------------------------------------------------
# Entry Point
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    main()
