import logging
import os
import subprocess
import argparse
from typing import List, Optional, Tuple, Dict, Any
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
    """
    try:
        asn_file_name: str = os.path.join(defaults.PATH_DICT[ASN_RPSBLAST_DIR], f"{species}.asn")

        rpsblast_command: List[str] = [
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

        result = subprocess.run(rpsblast_command, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f'Error running {command} for {species}: {result.stderr}')
            return None, None

        blast_formatter_command: List[str] = [
            'blast_formatter',
            '-archive',
            asn_file_name,
            '-outfmt',
            '6 std stitle',
        ]

        formatter_result: subprocess.CompletedProcess = subprocess.run(blast_formatter_command, capture_output=True, text=True)
        if formatter_result.returncode != 0:
            logging.error(f'Error running blast_formatter for {species}: {formatter_result.stderr}')
            return None, None

        blast_output: str = formatter_result.stdout

        return blast_output, asn_file_name

    except Exception as e:
        logging.error(f'An exception occurred while running RPS-BLAST for {species}: {str(e)}')
        return None, None


def parse_blast_output(blast_output: str) -> pd.DataFrame:
    """Parses the BLAST output into a DataFrame."""
    columns: List[str] = ['Query ID', 'Subject ID', 'Pct Identity', 'Alignment Length', 'Mismatches', 'Gap Openings',
                          'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value', 'Bit Score', 'Subject Title']
    if not blast_output:
        return pd.DataFrame(columns=columns)

    rows: List[List[str]] = [line.split('\t') for line in blast_output.strip().split('\n')]
    df: pd.DataFrame = pd.DataFrame(rows, columns=columns)
    return df


def process_rps_species(species: str) -> Optional[pd.DataFrame]:
    """Processes RPS-BLAST results for a single species."""
    fasta_file_path: str = os.path.join(defaults.PATH_DICT[FASTA_OUTPUT_DIR], f'{species}.fa')

    if not os.path.exists(fasta_file_path):
        logging.warning(f'FASTA file for species {species} does not exist at {fasta_file_path}.')
        return None

    blast_output: Optional[str]
    asn_file_name: Optional[str]
    blast_output, asn_file_name = rpsblaster(
        command='rpsblast',
        input_database_path=defaults.PATH_DICT[RPS_DB],
        query_file_path=fasta_file_path,
        species=species,
        evalue=defaults.E_VALUE_THRESHOLD
    )

    if not blast_output:
        logging.warning(f'No RPS-BLAST output for species {species}.')
        return None

    blast_df: pd.DataFrame = parse_blast_output(blast_output)
    if blast_df.empty:
        logging.warning(f'No data returned from RPS-BLAST for species {species}.')
        return None

    numeric_columns: List[str] = ['Pct Identity', 'E-value', 'Alignment Length', 'S. Start', 'S. End', 'Bit Score']
    col: str
    for col in numeric_columns:
        blast_df[col] = pd.to_numeric(blast_df[col], errors='coerce')

    blast_df.dropna(subset=numeric_columns, inplace=True)

    blast_df['Species'] = species

    return blast_df


def main(species_list: List[str]) -> None:
    colored_logging(log_file_name='rpsblast.txt')

    dfs: List[pd.DataFrame] = []

    with tqdm(total=len(species_list), desc='Running RPSBLAST against species...') as pbar:
        with concurrent.futures.ThreadPoolExecutor(max_workers=defaults.MAX_THREADPOOL_WORKERS) as executor:
            future_to_species: Dict[concurrent.futures.Future, str] = {
                executor.submit(process_rps_species, species): species
                for species in species_list
            }
            future: concurrent.futures.Future
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
        output_csv_path: str = os.path.join(defaults.PATH_DICT[TABLE_OUTPUT_DIR], 'rpsblast.csv')
        output_parquet_path: str = os.path.join(defaults.PATH_DICT[TABLE_OUTPUT_DIR], 'rpsblast.parquet')

        df.to_csv(output_csv_path, index=False)
        logging.info(f'Saved DataFrame as CSV to {output_csv_path}')

        df.to_parquet(output_parquet_path, index=False)
        logging.info(f'Saved DataFrame as Parquet to {output_parquet_path}')
    else:
        logging.warning('No data frames to concatenate. No output files were created.')


if __name__ == '__main__':
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description='Run RPS-BLAST for a list of species.')

    parser.add_argument(
        '--species_list',
        type=str,
        nargs='+',
        required=True,
        help='List of species to process.'
    )

    args: argparse.Namespace = parser.parse_args()

    species_list: List[str] = args.species_list

    main(species_list=species_list)
