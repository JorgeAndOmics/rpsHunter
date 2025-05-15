import concurrent.futures
import logging
import os
import subprocess
from typing import List, Optional, Tuple

import pandas as pd
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
        A tuple containing the combined BLAST output (with sequence) and the path to the ASN.1 file.
    """
    input_path = os.path.join(input_database_path, subject, subject)
    try:
        # Create ASN.1 file path
        asn_file_name = os.path.join(defaults.PATH_DICT[ASN_TBLASTN_DIR], f"{subject}.asn")

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

        # Run blast_formatter to get combined output (with sequence)
        blast_formatter_command = [
            'blast_formatter',
            '-archive',
            asn_file_name,
            '-outfmt',
            '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq'
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
    """Parses the combined BLAST output (with sequence) into a DataFrame."""
    if not blast_output:
        return pd.DataFrame()
    rows = [line.split('\t') for line in blast_output.strip().split('\n') if line]
    if not rows:
        return pd.DataFrame()
    # Determine the correct set of column names based on the number of fields
    if len(rows[0]) == 12:
        columns = ['Query ID', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches', 'Gap Openings',
                   'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value', 'Bit Score']
    elif len(rows[0]) == 13:
        columns = ['Query ID', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches', 'Gap Openings',
                   'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value', 'Bit Score', 'Subject Sequence']
    else:
        raise ValueError("Unexpected number of columns in BLAST output")
    df = pd.DataFrame(rows, columns=columns)
    return df


def process_species(species: str) -> Tuple[Optional[pd.DataFrame], List[SeqRecord]]:
    """Processes BLAST results for a single species."""
    blast_output, asn_file_name = blaster(
        command='tblastn',
        input_database_path=defaults.PATH_DICT[SPECIES_DB],
        query_file_path=defaults.QUERY_FILE,
        subject=species,
        evalue=defaults.E_VALUE_THRESHOLD
    )

    if not blast_output:
        logging.warning(f'No BLAST output for species {species}.')
        return None, []

    # Parse the combined BLAST output that includes sequence information
    blast_df = parse_blast_output(blast_output)
    if blast_df.empty:
        logging.warning(f'No data returned from BLAST for species {species}.')
        return None, []

    # Convert appropriate columns to numeric types
    blast_df = blast_df.astype({
        'Pct Identity': float,
        'E-value': float,
        'Alignment Length': int,
        'Bit Score': float,
        'S. Start': int,
        'S. End': int
    })

    # Clean up the sequence column if it exists
    if 'Subject Sequence' in blast_df.columns:
        blast_df['Subject Sequence'] = blast_df['Subject Sequence'].str.replace('-', '').str.replace('*', '')


    if blast_df.empty:
        logging.warning(f'No hits after filtering for species {species}.')
        return None, []

    blast_df['Species'] = species

    return blast_df


def main():
    # Set up logging
    colored_logging(log_file_name='blast.txt')

    dfs = []
    species_list = defaults.SPECIES

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {executor.submit(process_species, species): species for species in species_list}

        with tqdm(total=len(futures), desc = 'Running tBLASTn against select species...') as pbar:
            for future in concurrent.futures.as_completed(futures):
                species = futures[future]
                try:
                    blast_df = future.result()
                except Exception as exc:
                    logging.error(f"Exception processing {species}: {exc}")
                    pbar.update(1)
                    continue

                if blast_df is not None and not blast_df.empty:
                    dfs.append(blast_df)

                pbar.update(1)

    if dfs:
        df = pd.concat(dfs, axis=0, ignore_index=True)

        output_csv_path = os.path.join(defaults.PATH_DICT[TABLE_OUTPUT_DIR], 'blast.csv')
        output_parquet_path = os.path.join(defaults.PATH_DICT[TABLE_OUTPUT_DIR], 'blast.parquet')

        df.to_csv(output_csv_path, index=False)
        logging.info(f'Saved DataFrame as CSV to {output_csv_path}')

        df.to_parquet(output_parquet_path, index=False)
        logging.info(f'Saved DataFrame as Parquet to {output_parquet_path}')
    else:
        logging.warning('No data frames to concatenate. No output files were created.')


if __name__ == '__main__':
    main()
