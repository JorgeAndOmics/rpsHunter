"""
    Validation utilities for the rpsHunter pipeline.

    This module includes YAML config validation, query accession checks, FASTA content validation,
    dependency availability checks, and NCBI API key enforcement.
"""

# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------

import os
import re
import sys
import time
import argparse
import subprocess
import logging
from typing import List, Optional

import yamale
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import defaults
import validator
from colored_logging import colored_logging


# -----------------------------------------------------------------------------
# YAML VALIDATION
# -----------------------------------------------------------------------------

def yaml_validator(yaml_file: str, yaml_schema: str) -> bool:
    """
    Validates the YAML configuration file against a Yamale schema.

        Parameters
        ----------
            :param yaml_file: Path to the YAML configuration file.
            :param yaml_schema: Path to the YAML schema file.

        Returns
        -------
            :returns: True if the YAML is valid, False otherwise.

        Raises
        ------
            :raises yamale.YamaleError: If YAML validation fails.
    """
    try:
        schema = yamale.make_schema(yaml_schema)
        data = yamale.make_data(yaml_file)
        yamale.validate(schema, data)
        logging.info('YAML configuration file is valid.')
        return True
    except yamale.YamaleError as e:
        for results in e.results:
            for error in results.errors:
                logging.warning(f'Configuration error: {error}')
        return False


# -----------------------------------------------------------------------------
# NCBI QUERY ACCESSION VALIDATION
# -----------------------------------------------------------------------------

def query_validator(query_accession: str) -> bool:
    """
    Validates the query for NCBI accession validity.

        Parameters
        ----------
            :param query_accession: Query accession ID.

        Returns
        -------
            :returns: True if the accession is valid in NCBI, False otherwise.

        Raises
        ------
            :raises Exception: If NCBI API call fails.
    """
    logging.debug('Checking NCBI entry...')
    time.sleep(0.3)
    Entrez.email = defaults.ENTREZ_EMAIL
    try:
        with Entrez.esearch(db='protein', term=query_accession) as handle:
            record = Entrez.read(handle)
            return int(record['Count']) > 0
    except Exception:
        return False


# -----------------------------------------------------------------------------
# FASTA VALIDATION
# -----------------------------------------------------------------------------

def fasta_validator(fasta_file: str) -> bool:
    """
    Validates a FASTA file for content and headers.

        Parameters
        ----------
            :param fasta_file: Path to the FASTA file.

        Returns
        -------
            :returns: True if the FASTA is valid and non-empty, False otherwise.

        Raises
        ------
            :raises Exception: If parsing the file fails.
    """
    if not os.path.exists(fasta_file):
        logging.warning(f'FASTA file does not exist: {fasta_file}')
        return False

    try:
        records = list(SeqIO.parse(fasta_file, 'fasta'))
        if not records:
            logging.warning('FASTA file is empty or has no valid records.')
            return False

        for i, record in enumerate(records):
            if not record.id:
                logging.warning(f'Record {i + 1} is missing a header.')
                return False

        logging.info(f'FASTA file {fasta_file} is valid.')
        return True

    except Exception as e:
        logging.warning(f'Error parsing FASTA: {e}')
        return False


# -----------------------------------------------------------------------------
# TOOL AVAILABILITY VALIDATION
# -----------------------------------------------------------------------------

def validate_programs() -> bool:
    """
    Validates required external programs are available in the system path.

        Returns
        -------
            :returns: True if all tools are available, False otherwise.
    """

    def check_version(cmd: str, version_cmd: str) -> bool:
        try:
            subprocess.run([cmd, version_cmd], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    checks = {
        'BLAST+': check_version('tblastn', '-version'),
        'rpsbproc': check_version('rpsbproc', '-version'),
        'Datasets': check_version('datasets', '--version')
    }

    for tool, status in checks.items():
        log_fn = logging.info if status else logging.warning
        log_fn(f'{tool} {"is installed." if status else "not found. Please install it."}')

    return all(checks.values())


# -----------------------------------------------------------------------------
# NCBI API KEY VALIDATION
# -----------------------------------------------------------------------------

def validate_ncbi_key() -> None:
    """
    Ensures that an NCBI API key is available in the environment.
    Prompts the user to enter one if not present.

        Returns
        -------
            :returns: None. Modifies environment if user provides input.
    """
    if 'NCBI_API_KEY' not in os.environ:
        logging.warning('No NCBI API key found. You can set it with: export NCBI_API_KEY="your_api_key".')
        if (
            api_key := input('Enter your NCBI API key [Leave empty to skip]: ')
            or None
        ):
            os.environ['NCBI_API_KEY'] = api_key
            logging.info('NCBI API key set in environment variables.')
        else:
            logging.warning('No NCBI API key provided. Expect slower NCBI retrievals.')
    else:
        logging.info('NCBI API key is set in the environment variables.')


# -----------------------------------------------------------------------------
# MASTER VALIDATION ORCHESTRATOR
# -----------------------------------------------------------------------------

def main_validator(fasta_files: Optional[List[str]]) -> bool:
    """
    Orchestrates validation for YAML, Query accession, FASTA, external tools, and API key.

        Parameters
        ----------
            :param fasta_files: List of FASTA file paths to validate.

        Returns
        -------
            :returns: True if all checks pass, False otherwise.
    """
    logging.debug('Starting input validation process...')

    yaml_ok = yaml_validator(
        yaml_schema=os.path.join(defaults.PATH_DICT['CONFIG_DIR'], 'schema.yaml'),
        yaml_file=defaults.CONFIG_FILE
    )

    query_ok = query_validator(query_accession=defaults.config['query']['accession'])

    if not defaults.USE_SPECIES_DICT and fasta_files:
        fasta_results = [fasta_validator(f) for f in fasta_files]
        fasta_ok = all(fasta_results)
    else:
        fasta_ok = True

    programs_ok = validate_programs()
    validate_ncbi_key()

    return all([yaml_ok, query_ok, fasta_ok, programs_ok])


# -----------------------------------------------------------------------------
# INTERACTIVE CONFIRMATION
# -----------------------------------------------------------------------------

def green_light(all_valid: bool) -> bool:
    """
    Asks user to confirm whether to proceed if all validations passed.

        Parameters
        ----------
            :param all_valid: Whether all validations were successful.

        Returns
        -------
            :returns: True if user wants to proceed, False otherwise.
    """
    if not all_valid:
        logging.warning(f'Settings validation failed. Check logs at {defaults.PATH_DICT["LOG_DIR"]}.')
        return False

    logging.info('All systems green, ready to rock.')
    time.sleep(0.1)

    proceed = input('Proceed [Y/n]: ') or 'Y'
    if proceed.upper() == 'Y':
        logging.info('rpsHunter started. Depending on your system, this may take some time.')
        return True
    elif proceed.upper() == 'N':
        logging.warning('Workflow aborted by user.')
        return False
    else:
        logging.warning('Invalid input. Please try again.')
        return green_light(all_valid)


# -----------------------------------------------------------------------------
# VALIDATION ENTRYPOINT
# -----------------------------------------------------------------------------

def validation_run(fasta_files: Optional[List[str]] = None) -> bool:
    """
    CLI entrypoint to trigger validation routines and prompt user to continue.

        Parameters
        ----------
            :param fasta_files: Optional list of FASTA file paths to validate.

        Returns
        -------
            :returns: True if user confirms execution after passing validation, False otherwise.
    """
    colored_logging.colored_logging(log_file_name='validator.log')

    fasta_files_list = fasta_files or []
    all_valid = main_validator(fasta_files=fasta_files_list)

    return green_light(all_valid=all_valid)
