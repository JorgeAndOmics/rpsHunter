"""
    Converts RPS-BLAST ASN.1 output files to readable text using the `rpsbproc` tool.

    Iterates over all `.asn` files in the input directory, processes them with `rpsbproc`,
    and writes the resulting `.txt` files to the output directory.
"""

import os
import subprocess
from typing import List

import defaults


# -----------------------------------------------------------------------------
# RPSBPROC Conversion Function
# -----------------------------------------------------------------------------

def main(input_dir: str, output_dir: str, db_path: str, t_option: str) -> None:
    """
    Processes all ASN.1 files in the input directory using `rpsbproc` and writes text outputs.

        Parameters
        ----------
            :param input_dir: Directory containing `.asn` files to be processed.
            :param output_dir: Destination directory for output `.txt` files.
            :param db_path: Path to the conserved domain database for `rpsbproc`.
            :param t_option: Option for the `-t` flag in `rpsbproc` (e.g., 'both', 'domain', etc.).

        Returns
        -------
            :returns: None. Writes processed `.txt` files to disk.

        Raises
        ------
            :raises FileNotFoundError: If the input directory does not exist.
            :raises subprocess.CalledProcessError: If the `rpsbproc` command fails.
    """
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f'Input directory does not exist: {input_dir}')

    input_files: List[str] = [
        f for f in os.listdir(input_dir) if f.endswith('.asn')
    ]

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
# Entry Point
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    main(
        input_dir=defaults.ASN_RPSBLAST_DIR,
        output_dir=defaults.RPSBPROC_OUTPUT_DIR,
        db_path=defaults.RPSBPROC_DB,
        t_option='both'
    )
