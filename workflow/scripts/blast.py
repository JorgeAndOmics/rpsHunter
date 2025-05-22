"""
    Runs a multi-species BLAST workflow using `tblastn`, collects results, and exports to CSV and Parquet.

    Performs parallel BLAST searches across species-specific databases, formats the results,
    and saves combined outputs to disk.

    Requires definitions from a `defaults` module and `colored_logging` for log setup.
"""

import concurrent.futures
import logging
import os
import subprocess
from typing import List, Optional, Tuple

import pandas as pd
from tqdm import tqdm

import defaults
from colored_logging import colored_logging


# -----------------------------------------------------------------------------
# BLAST Runner
# -----------------------------------------------------------------------------

def blaster(
    command: str,
    input_database_path: str,
    query_file_path: str,
    subject: str,
    evalue: float,
) -> Tuple[Optional[str], Optional[str]]:
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
    input_path = os.path.join(input_database_path, subject, subject)
    try:
        asn_file_name = os.path.join(defaults.PATH_DICT['ASN_TBLASTN_DIR'], f'{subject}.asn')

        # Construct the BLAST command
        blast_command = [
            command,
            '-db', input_path,
            '-query', query_file_path,
            '-evalue', str(evalue),
            '-outfmt', '11',
            '-out', asn_file_name
        ]

        result = subprocess.run(blast_command, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f'Error running {command} for {subject}: {result.stderr}')
            return None, None

        # Format ASN.1 to tabular output with sequence data
        blast_formatter_command = [
            'blast_formatter',
            '-archive', asn_file_name,
            '-outfmt',
            '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq'
        ]

        formatter_result = subprocess.run(blast_formatter_command, capture_output=True, text=True)

        if formatter_result.returncode != 0:
            logging.error(f'Error running blast_formatter for {subject}: {formatter_result.stderr}')
            return None, None

        return formatter_result.stdout, asn_file_name

    except Exception as e:
        logging.error(f'An exception occurred while running BLAST for {subject}: {str(e)}')
        return None, None


# -----------------------------------------------------------------------------
# BLAST Output Parser
# -----------------------------------------------------------------------------

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

def process_species(species: str) -> Optional[pd.DataFrame]:
    """
    Executes and parses BLAST for a single species.

        Parameters
        ----------
            :param species: Species name to process.

        Returns
        -------
            :returns: A DataFrame with parsed BLAST results or None if unsuccessful.

        Raises
        ------
            :raises Exception: Propagates exceptions from BLAST or parsing steps.
    """
    blast_output, asn_file_name = blaster(
        command='tblastn',
        input_database_path=defaults.PATH_DICT['SPECIES_DB'],
        query_file_path=defaults.QUERY_FILE,
        subject=species,
        evalue=defaults.E_VALUE_THRESHOLD
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

    # Clean up sequences if applicable
    if 'Subject Sequence' in blast_df.columns:
        blast_df['Subject Sequence'] = blast_df['Subject Sequence'].str.replace('-', '').str.replace('*', '')

    if blast_df.empty:
        logging.warning(f'No hits after filtering for species {species}.')
        return None

    blast_df['Species'] = species
    return blast_df


# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------

def main() -> None:
    """
    Executes the multi-threaded BLAST workflow and saves results to disk.

        Parameters
        ----------
            :param None

        Returns
        -------
            :returns: None. Writes results to CSV and Parquet files.

        Raises
        ------
            :raises IOError: If the result files cannot be written to disk.
    """
    colored_logging(log_file_name='blast.txt')
    dfs: List[pd.DataFrame] = []
    species_list: List[str] = defaults.SPECIES

    # Execute BLASTs in parallel
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(process_species, species): species
            for species in species_list
        }

        with tqdm(total=len(futures), desc='Running tBLASTn against select species...') as pbar:
            for future in concurrent.futures.as_completed(futures):
                species = futures[future]
                try:
                    blast_df = future.result()
                except Exception as exc:
                    logging.error(f'Exception processing {species}: {exc}')
                    pbar.update(1)
                    continue

                if blast_df is not None and not blast_df.empty:
                    dfs.append(blast_df)

                pbar.update(1)

    # Save outputs if any data collected
    if dfs:
        df = pd.concat(dfs, axis=0, ignore_index=True)

        output_csv_path = os.path.join(defaults.PATH_DICT['TABLE_OUTPUT_DIR'], 'blast.csv')
        output_parquet_path = os.path.join(defaults.PATH_DICT['TABLE_OUTPUT_DIR'], 'blast.parquet')

        df.to_csv(output_csv_path, index=False)
        logging.info(f'Saved DataFrame as CSV to {output_csv_path}')

        df.to_parquet(output_parquet_path, index=False)
        logging.info(f'Saved DataFrame as Parquet to {output_parquet_path}')
    else:
        logging.warning('No data frames to concatenate. No output files were created.')


# -----------------------------------------------------------------------------
# Entry Point
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    main()
