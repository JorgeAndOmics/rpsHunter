"""
    CLI entrypoint for rpsHunter: a pipeline for gene domain identification.

    Provides interface for executing discrete Snakemake workflow components,
    standardizing FASTA filenames, and validating required input files.
"""

# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------

import os
import re
import sys
import logging
import subprocess
import argparse
from pathlib import Path
from typing import List

import defaults
import colored_logging
from validator import validation_run


# -----------------------------------------------------------------------------
# FASTA EXTENSION STANDARDIZATION
# -----------------------------------------------------------------------------

def standardize_fasta_extensions(fasta_dir_path: str) -> None:
    """
    Standardize extensions of all FASTA files in the provided directory to .fa.

        Parameters
        ----------
            :param fasta_dir_path: Path to the directory containing FASTA files
                                   with extensions like .fasta, .fna, or .fas.

        Returns
        -------
            :returns: None. Files are renamed in-place.

        Raises
        ------
            :raises FileNotFoundError: If the specified directory does not exist.
            :raises OSError: If renaming a file fails.
    """
    if not os.path.isdir(fasta_dir_path):
        raise FileNotFoundError(f'Directory not found: {fasta_dir_path}')

    pattern = re.compile(r'\.(fasta|fna|fas)$', re.IGNORECASE)

    for file in Path(fasta_dir_path).iterdir():
        if file.is_file() and pattern.search(file.name):
            new_name = file.with_name(f'{file.stem}.fa')
            logging.debug(f'Renaming: {file.name} -> {new_name.name}')
            file.rename(new_name)


# -----------------------------------------------------------------------------
# SNAKEMAKE RULE EXECUTION
# -----------------------------------------------------------------------------

def run_snakemake_rule(rule: str, num_cores: int, display_info: bool, snakemake_flags: List[str] = None) -> None:
    """
    Execute a Snakemake rule with specified options.

        Parameters
        ----------
            :param rule: Name of the Snakemake rule to execute.
            :param num_cores: Number of cores to allocate for the rule.
            :param display_info: Whether to display detailed Snakemake command output.
            :param snakemake_flags: Optional list of extra flags to pass to Snakemake.

        Returns
        -------
            :returns: None. Exits on error.

        Raises
        ------
            :raises subprocess.CalledProcessError: If the Snakemake command fails.
    """
    if snakemake_flags is None:
        snakemake_flags = []

    shell_cmd: List[str] = [
        'snakemake',
        rule,
        '--cores', str(num_cores),
        '--rerun-incomplete',
        *snakemake_flags
    ]

    if not display_info:
        shell_cmd.append('-q')

    try:
        result = subprocess.run(shell_cmd)
        if result.returncode != 0:
            logging.error(f'Snakemake rule "{rule}" failed with code {result.returncode}')
            sys.exit(result.returncode)
    except subprocess.CalledProcessError as e:
        logging.error(f'Snakemake rule execution error: {e.stderr.decode()}')
        sys.exit(1)


# -----------------------------------------------------------------------------
# CLI ENTRYPOINT
# -----------------------------------------------------------------------------

def cli_entry() -> None:
    """
    Main entrypoint for the rpsHunter CLI.

        Parameters
        ----------
            :param None

        Returns
        -------
            :returns: None. Executes selected Snakemake rules.

        Raises
        ------
            :raises SystemExit: If validation fails or CLI arguments are missing.
    """
    colored_logging.colored_logging(log_file_name='rpsHunter_main.log')

    # -------------------------------------------------------------------------
    # Argument Parsing
    # -------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description='rpsHunter: A tool for gene domain identification and analysis.'
    )

    parser.add_argument('--download_genomes', action='store_true',
                        help='Download genomes listed in the .yaml file.')

    parser.add_argument('--download_query', action='store_true',
                        help='Download the query protein sequence in FASTA format.')

    parser.add_argument('--blast_dbs', action='store_true',
                        help='Generate BLAST databases from downloaded genomes.')

    parser.add_argument('--blast', action='store_true',
                        help='Run tBLASTn of the query sequence against genome databases.')

    parser.add_argument('--rpsblast', action='store_true',
                        help='Run RPS-BLAST against a domain database.')

    parser.add_argument('--rpsbproc', action='store_true',
                        help='Process RPS-BLAST ASN files with rpsbproc.')

    parser.add_argument('--rpsbproc_parser', action='store_true',
                        help='Parse rpsbproc output into tabular format.')

    parser.add_argument('--contingency_parser', action='store_true',
                        help='Generate contingency tables and visual plots.')

    parser.add_argument('--pair_detection', action='store_true',
                        help='Detect conserved domain co-occurrence patterns.')

    parser.add_argument('--skip_validation', '-skp', action='store_true',
                        help='Skip input validation checks.')

    args, unknown = parser.parse_known_args()

    # -------------------------------------------------------------------------
    # Standardize FASTA File Extensions
    # -------------------------------------------------------------------------
    standardize_fasta_extensions(defaults.PATH_DICT['SPECIES_DB'])

    # -------------------------------------------------------------------------
    # Validate Input Files
    # -------------------------------------------------------------------------
    species_paths: List[str] = [
        os.path.join(defaults.PATH_DICT['SPECIES_DB'], f'{species}.fa')
        for species in defaults.SPECIES
    ]

    if not any(vars(args).values()):
        parser.print_help()
        sys.exit(0)

    validation: bool = True if args.skip_validation else validation_run(species_paths)

    if not validation:
        sys.exit(1)

    # -------------------------------------------------------------------------
    # Dispatch Snakemake Rules
    # -------------------------------------------------------------------------
    if args.download_genomes:
        run_snakemake_rule(
            rule='genome_downloader',
            num_cores=defaults.NUM_CORES,
            display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
            snakemake_flags=unknown
        )

    if args.download_query:
        run_snakemake_rule(
            rule='query_downloader',
            num_cores=defaults.NUM_CORES,
            display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
            snakemake_flags=unknown
        )

    if args.blast_dbs:
        run_snakemake_rule(
            rule='blast_db_generator',
            num_cores=defaults.NUM_CORES,
            display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
            snakemake_flags=unknown
        )

    if args.blast:
        run_snakemake_rule(
            rule='blaster_parser',
            num_cores=defaults.NUM_CORES,
            display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
            snakemake_flags=unknown
        )

    if args.rpsblast:
        run_snakemake_rule(
            rule='rpsblaster',
            num_cores=defaults.NUM_CORES,
            display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
            snakemake_flags=unknown
        )

    if args.rpsbproc:
        run_snakemake_rule(
            rule='rpsbproc',
            num_cores=defaults.NUM_CORES,
            display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
            snakemake_flags=unknown
        )

    if args.rpsbproc_parser:
        run_snakemake_rule(
            rule='rpsbproc_parser',
            num_cores=defaults.NUM_CORES,
            display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
            snakemake_flags=unknown
        )

    if args.contingency_parser:
        run_snakemake_rule(
            rule='contingency_sorter',
            num_cores=defaults.NUM_CORES,
            display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
            snakemake_flags=unknown
        )

    if args.pair_detection:
        run_snakemake_rule(
            rule='pair_detector',
            num_cores=defaults.NUM_CORES,
            display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
            snakemake_flags=unknown
        )

    logging.info('Finished execution.')


# -----------------------------------------------------------------------------
# ENTRY POINT
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    cli_entry()
