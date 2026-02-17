"""
    Runs a multi-species BLAST workflow using `tblastn`, collects results, and exports to CSV and Parquet.

    Performs parallel BLAST searches across species-specific databases, formats the results,
    and saves combined outputs to disk.

    Requires definitions from a `defaults` module and `colored_logging` for log setup.
"""

import argparse
import logging
import random
import string
import subprocess
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd

import defaults
from colored_logging import colored_logging


# -----------------------------------------------------------------------------
# BLAST Runner
# -----------------------------------------------------------------------------

def blaster(
    command: str,
    input_database_path: Path,
    query_file_path: Path,
    subject: str,
    evalue: float,
    asn_output_path: Path = None,
) -> Tuple[Optional[str], Optional[Path]]:
    """
    Runs a BLAST search for a given subject species against its genome database.

        Parameters
        ----------
            :param command: The command to run BLAST (e.g., 'tblastn').
            :param input_database_path: Path to the directory containing subject databases.
            :param query_file_path: Path to the query file (e.g., a protein FASTA).
            :param subject: Name of the species to be queried.
            :param evalue: E-value threshold for filtering BLAST results.

        Returns
        -------
            :returns: A tuple containing the BLAST output (tabular string with sequence)
                      and the path to the resulting ASN.1 archive file.

        Raises
        ------
            :raises subprocess.CalledProcessError: If either BLAST or formatter fails.
            :raises Exception: For any other unexpected runtime issues.
    """
    input_path = input_database_path / subject / subject
    try:
        if asn_output_path is None:
            asn_file_path = defaults.PATH_DICT['ASN_TBLASTN_DIR'] / f'{subject}.asn'
        else:
            asn_file_path = asn_output_path
            asn_file_path.parent.mkdir(parents=True, exist_ok=True)

        # Construct the BLAST command
        blast_command = [
            command,
            '-db', str(input_path),
            '-query', str(query_file_path),
            '-evalue', str(evalue),
            '-outfmt', '11',
            '-out', str(asn_file_path)
        ]

        result = subprocess.run(blast_command, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f'Error running {command} for {subject}: {result.stderr}')
            return None, None

        # Format ASN.1 to tabular output with sequence data
        blast_formatter_command = [
            defaults.BLAST_FORMATTER_CMD,
            '-archive', str(asn_file_path),
            '-outfmt',
            '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq'
        ]

        formatter_result = subprocess.run(blast_formatter_command, capture_output=True, text=True)

        if formatter_result.returncode != 0:
            logging.error(f'Error running blast_formatter for {subject}: {formatter_result.stderr}')
            return None, None

        return formatter_result.stdout, asn_file_path

    except Exception as e:
        logging.error(f'An exception occurred while running BLAST for {subject}: {str(e)}')
        return None, None


# -----------------------------------------------------------------------------
# BLAST Output Parser
# -----------------------------------------------------------------------------

def random_string_generator(length: int) -> str:
    """
    Generates a random string of the specified length.

        Parameters
        ----------
        :param length: The length of the string to generate.

        Returns
        -------
        :returns: The generated random string.

    """
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))


def parse_blast_output(blast_output: str) -> pd.DataFrame:
    """
    Parses the combined BLAST output (including subject sequences) into a DataFrame.

        Parameters
        ----------
            :param blast_output: A string representing tabular BLAST results with sequences.

        Returns
        -------
            :returns: A pandas DataFrame containing BLAST hits.

        Raises
        ------
            :raises ValueError: If output contains an unexpected number of columns.
    """
    if not blast_output:
        return pd.DataFrame()

    rows = [line.split('\t') for line in blast_output.strip().split('\n') if line]
    if not rows:
        return pd.DataFrame()

    # Determine column schema
    if len(rows[0]) == 12:
        columns = [
            'Query ID', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches', 'Gap Openings',
            'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value', 'Bit Score'
        ]
    elif len(rows[0]) == 13:
        columns = [
            'Query ID', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches', 'Gap Openings',
            'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value', 'Bit Score', 'Subject Sequence'
        ]
    else:
        raise ValueError('Unexpected number of columns in BLAST output')

    return pd.DataFrame(rows, columns=columns)


# -----------------------------------------------------------------------------
# Per-Species BLAST Processing
# -----------------------------------------------------------------------------

def process_species(species: str, query_file: Path = None, query_accession: str = '',
                    asn_output_path: Path = None) -> Optional[pd.DataFrame]:
    """
    Executes and parses BLAST for a single species.

        Parameters
        ----------
            :param species: Species name to process.
            :param query_file: Path to the query FASTA file.
            :param query_accession: Query accession ID for provenance tracking.
            :param asn_output_path: Explicit path for ASN output (overrides default).

        Returns
        -------
            :returns: A DataFrame with parsed BLAST results or None if unsuccessful.

        Raises
        ------
            :raises Exception: Propagates exceptions from BLAST or parsing steps.
    """
    if query_file is None:
        query_file = defaults.QUERY_FILE

    blast_output, asn_file_name = blaster(
        command=defaults.TBLASTN_CMD,
        input_database_path=defaults.PATH_DICT['SPECIES_DB'],
        query_file_path=query_file,
        subject=species,
        evalue=defaults.E_VALUE_THRESHOLD,
        asn_output_path=asn_output_path,
    )

    if not blast_output:
        logging.warning(f'No BLAST output for species {species}.')
        return None

    blast_df = parse_blast_output(blast_output)

    if blast_df.empty:
        logging.warning(f'No data returned from BLAST for species {species}.')
        return None

    # Convert numeric columns
    numeric_columns = {
        'Pct Identity': float,
        'E-value': float,
        'Alignment Length': int,
        'Bit Score': float,
        'S. Start': int,
        'S. End': int
    }
    blast_df = blast_df.astype(numeric_columns)

    if blast_df.empty:
        logging.warning(f'No hits after filtering for species {species}.')
        return None

    blast_df['Species'] = species
    blast_df['Species_Name'] = defaults.SPECIES_DICT.get(species, species)
    blast_df['Query_Accession'] = query_accession

    # Add a random string column for each row
    blast_df['Tag'] = [random_string_generator(defaults.RANDOM_ID_LENGTH) for _ in range(len(blast_df))]

    # Remove duplicates
    blast_df = blast_df.drop_duplicates().reset_index(drop=True)

    return blast_df


# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------

def main() -> None:
    """
    Executes tBLASTn for a single species and writes per-species outputs.
    """
    colored_logging(log_file_name='blast.txt')

    parser = argparse.ArgumentParser()
    parser.add_argument('--species', required=True, help='Single species to process')
    parser.add_argument('--query-file', type=str, default=None,
                        help='Path to query FASTA file (default: from config).')
    parser.add_argument('--query-accession', type=str, default='',
                        help='Query accession ID for provenance tracking.')
    parser.add_argument('--output-parquet', type=str, default=None,
                        help='Output parquet path.')
    parser.add_argument('--output-csv', type=str, default=None,
                        help='Output CSV path.')
    parser.add_argument('--output-asn', type=str, default=None,
                        help='Output ASN path.')
    args = parser.parse_args()

    query_file = Path(args.query_file) if args.query_file else defaults.QUERY_FILE
    asn_out = Path(args.output_asn) if args.output_asn else None
    blast_df = process_species(args.species, query_file=query_file,
                               query_accession=args.query_accession,
                               asn_output_path=asn_out)

    # Resolve output paths (explicit CLI args or legacy defaults)
    if args.output_parquet:
        parquet_path = Path(args.output_parquet)
    else:
        parquet_path = defaults.PATH_DICT['BLAST_SPECIES_DIR'] / f'{args.species}.parquet'

    if args.output_csv:
        csv_path = Path(args.output_csv)
    else:
        csv_path = defaults.PATH_DICT['BLAST_SPECIES_DIR'] / f'{args.species}.csv'

    if args.output_asn:
        asn_path = Path(args.output_asn)
    else:
        asn_path = defaults.PATH_DICT['ASN_TBLASTN_DIR'] / f'{args.species}.asn'

    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    asn_path.parent.mkdir(parents=True, exist_ok=True)

    if blast_df is not None and not blast_df.empty:
        blast_df.to_parquet(parquet_path, index=False)
        blast_df.to_csv(csv_path, index=False)
        logging.info(f'Saved {len(blast_df)} rows for {args.species}')
    else:
        # Write an empty parquet with the correct schema so downstream rules don't crash.
        empty_columns = [
            'Query ID', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches',
            'Gap Openings', 'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value',
            'Bit Score', 'Subject Sequence', 'Species', 'Species_Name', 'Query_Accession', 'Tag'
        ]
        pd.DataFrame(columns=empty_columns).to_parquet(parquet_path, index=False)
        csv_path.touch()
        # Touch the ASN only if tBLASTn did not already write it (paths 2 & 3).
        if not asn_path.exists():
            asn_path.touch()
        logging.warning(f'No BLAST results for {args.species} — empty outputs written.')


# -----------------------------------------------------------------------------
# Entry Point
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    main()
