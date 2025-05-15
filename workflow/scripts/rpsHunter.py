# -------------------
# DEPENDENCIES
# -------------------

import os
import re
import sys
import logging
import subprocess
import argparse
from pathlib import Path

from validator import validation_run
import defaults
import colored_logging


# -----------------------------
# FASTA EXTENSION STANDARDIZATION
# -----------------------------

def standardize_fasta_extensions(fasta_dir_path: str) -> None:
    """
    Standardize extensions of all FASTA files in the provided directory to .fa.

    Parameters
    ----------
    fasta_dir_path : str
        Path to the directory containing FASTA files with various extensions (.fasta, .fna, .fas).
    """
    pattern = re.compile(r'\.(fasta|fna|fas)$', re.IGNORECASE)

    for file in Path(fasta_dir_path).iterdir():
        if file.is_file() and pattern.search(file.name):
            new_name: str = file.with_name(f'{file.stem}.fa')
            logging.debug(f"Renaming: {file.name} -> {new_name.name}")
            file.rename(new_name)


# -----------------------------
# SNAKEMAKE RULE EXECUTION
# -----------------------------

def run_snakemake_rule(rule: str, num_cores: int, display_info: bool, snakemake_flags: list[str] = None) -> None:
    """
    Execute a Snakemake rule with specified options.

    Parameters
    ----------
        :param rule : str
        Name of the Snakemake rule to execute.
        :param num_cores : int
        Number of cores to allocate for the rule.
        :param display_info : bool
        Whether to display detailed Snakemake command output.
        :param snakemake_flags:
    """
    if snakemake_flags is None:
        snakemake_flags = []
    shell_cmd: list = [
        "snakemake",
        rule,
        "--cores",
        str(num_cores),
        "--rerun-incomplete",
        *snakemake_flags
    ]

    if not display_info:
        shell_cmd.append("-q")

    try:
        result = subprocess.run(shell_cmd)
        if result.returncode != 0:
            logging.error(f"Snakemake rule '{rule}' failed with code {result.returncode}")
            logging.error(result.stderr)
            sys.exit(result.returncode)
        else:
            return
    except subprocess.CalledProcessError as e:
        logging.error(f"Snakemake rule failed with error: {e.stderr.decode()}")
        sys.exit(1)


# -----------------------------
# CLI ENTRYPOINT
# -----------------------------

def cli_entry() -> None:
    """
    Main entrypoint for RetroSeek CLI.

    Handles argument parsing, input validation, and Snakemake rule dispatch.
    """
    colored_logging.colored_logging(log_file_name='RetroSeek_main.log')

    # -----------------------------
    # PARSE COMMAND-LINE ARGUMENTS
    # -----------------------------
    parser = argparse.ArgumentParser(
        description='rpsHunter: A tool for gene domain identification and analysis.'
    )

    parser.add_argument(
        '--download_genomes',
        action='store_true',
        help="Downloads genomes specified in .yaml file. "
             "Will skip download if the same files exist currently in the target directory."
    )

    parser.add_argument(
        '--download_query',
        action='store_true',
        help="Downloads query protein sequence in fasta format."
             " Will skip download if the same file exists currently in the target directory."
    )

    parser.add_argument(
        '--blast_dbs',
        action='store_true',
        help="Generates BLAST databases from genomes."
    )

    parser.add_argument(
        '--blast',
        action='store_true',
        help="Runs tBLASTn of provided protein accession sequence against genome databases."
    )

    parser.add_argument(
        '--rpsblast',
        action='store_true',
        help="Runs rpsblast of provided protein accession sequence against domain database."
    )

    parser.add_argument(
        '--rpsbproc',
        action='store_true',
        help="Runs rpsbproc to process the rpsblast output."
    )

    parser.add_argument(
        '--rpsbproc_parser',
        action='store_true',
        help="Parses the rpsbproc output to table format."
    )

    parser.add_argument(
        '--contingency_parser',
        action='store_true',
        help="Generates contingency tables and plots."
    )

    parser.add_argument(
        '--skip_validation', '-skp',
        action='store_true',
        help='Skip input validation.'
    )

    args, unknown = parser.parse_known_args()

    # -----------------------------
    # STANDARDIZE FASTA EXTENSIONS
    # -----------------------------
    standardize_fasta_extensions(defaults.PATH_DICT["SPECIES_DB"])

    # -----------------------------
    # VALIDATE INPUT & EXECUTE RULES
    # -----------------------------
    species_paths: list = [
        os.path.join(defaults.PATH_DICT["SPECIES_DB"], f'{species}.fa')
        for species in defaults.SPECIES
    ]

    # Print help message if no arguments or wrong arguments are provided
    if not any(vars(args).values()):
        parser.print_help()
        sys.exit(0)


    validation = True if args.skip_validation else validation_run(species_paths)
    if validation:
        if args.download_genomes:
            run_snakemake_rule('genome_downloader',
                               num_cores=defaults.NUM_CORES,
                               display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                               snakemake_flags=unknown)

        if args.download_query:
            run_snakemake_rule('query_downloader',
                               num_cores=defaults.NUM_CORES,
                               display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                               snakemake_flags=unknown)

        if args.blast_dbs:
            run_snakemake_rule('blast_db_generator',
                               num_cores=defaults.NUM_CORES,
                               display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                               snakemake_flags=unknown)

        if args.blast:
            run_snakemake_rule('blaster_parser',
                               num_cores=defaults.NUM_CORES,
                               display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                               snakemake_flags=unknown)

        if args.rpsblast:
            run_snakemake_rule('rpsblaster',
                               num_cores=defaults.NUM_CORES,
                               display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                               snakemake_flags=unknown)

        if args.rpsbproc:
            run_snakemake_rule('rpsbproc',
                               num_cores=defaults.NUM_CORES,
                               display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                               snakemake_flags=unknown)

        if args.rpsbproc_parser:
            run_snakemake_rule('rpsbproc_parser',
                               num_cores=defaults.NUM_CORES,
                               display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                               snakemake_flags=unknown)

        if args.contingency_parser:
            run_snakemake_rule('contingency_sorter',
                               num_cores=defaults.NUM_CORES,
                               display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                               snakemake_flags=unknown)

        if args.pair_detection:
            run_snakemake_rule('pair_detector',
                               num_cores=defaults.NUM_CORES,
                               display_info=defaults.DISPLAY_SNAKEMAKE_INFO,
                               snakemake_flags=unknown)

    else:
        sys.exit(1)

    logging.info("Finished execution.")


# -----------------------------
# EXECUTION TRIGGER
# -----------------------------

if __name__ == "__main__":
    cli_entry()
