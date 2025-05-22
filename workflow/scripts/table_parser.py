"""
Module to process ASN files using the `rpsbproc` command-line tool.

Reads `.asn` files from the specified input directory, processes them using `rpsbproc`
with a given database and output directory, and writes the results as `.txt` files.
"""

import os
import subprocess
from typing import List

import defaults


# -----------------------------------------------------------------------------
# Main Processing Function
# -----------------------------------------------------------------------------

def main(input_dir: str, output_dir: str, db_path: str, t_option: str) -> None:
    """
    Processes all `.asn` files in the input directory using `rpsbproc`.

        Parameters
        ----------
            :param input_dir: Directory containing input `.asn` files.
            :param output_dir: Directory where output `.txt` files will be written.
            :param db_path: Path to the RPSBProc database file.
            :param t_option: Option to pass to the `-t` flag of `rpsbproc`.

        Returns
        -------
            :returns: None. Writes `.txt` output files for each `.asn` input.

        Raises
        ------
            :raises FileNotFoundError: If the input directory does not exist.
            :raises subprocess.CalledProcessError: If the `rpsbproc` command fails.
    """

    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f'Input directory does not exist: {input_dir}')

    # -------------------------------------------------------------------------
    # Collect input ASN files from the specified input directory
    # -------------------------------------------------------------------------
    input_files: List[str] = [
        f for f in os.listdir(input_dir) if f.endswith('.asn')
    ]

    # -------------------------------------------------------------------------
    # Process each input file using the rpsbproc command-line tool
    # -------------------------------------------------------------------------
    for filename in input_files:
        input_file: str = os.path.join(input_dir, filename)
        output_filename: str = f'{os.path.splitext(filename)[0]}.txt'
        output_file: str = os.path.join(output_dir, output_filename)

        cmd: List[str] = [
            'rpsbproc',
            '-i', input_file,
            '-d', db_path,
            '-t', t_option,
            '-o', output_file
        ]

        print(f'Processing {input_file} -> {output_file}')
        subprocess.run(cmd, check=True)


# -----------------------------------------------------------------------------
# Script Entry Point
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    main(
        input_dir=defaults.ASN_RPSBLAST_DIR,
        output_dir=defaults.RPSBPROC_OUTPUT_DIR,
        db_path=defaults.RPSBPROC_DB,
        t_option='both'
    )
