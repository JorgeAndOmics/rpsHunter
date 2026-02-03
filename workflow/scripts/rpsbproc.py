"""
    Converts RPS-BLAST ASN.1 output files to readable text using the `rpsbproc` tool.

    Iterates over all `.asn` files in the input directory, processes them with `rpsbproc`,
    and writes the resulting `.txt` files to the output directory.
"""

import subprocess
from pathlib import Path
from typing import List

import defaults


# -----------------------------------------------------------------------------
# RPSBPROC Conversion Function
# -----------------------------------------------------------------------------

def main(input_dir: Path, output_dir: Path, db_path: Path, t_option: str) -> None:
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
    if not input_dir.is_dir():
        raise FileNotFoundError(f'Input directory does not exist: {input_dir}')

    input_files: List[Path] = list(input_dir.glob('*.asn'))

    for input_file in input_files:
        output_file: Path = output_dir / f'{input_file.stem}.txt'

        cmd: List[str] = [
            'rpsbproc',
            '-i', str(input_file),
            '-d', str(db_path),
            '-t', t_option,
            '-o', str(output_file)
        ]

        print(f'Processing {input_file} -> {output_file}')
        subprocess.run(cmd, check=True)


# -----------------------------------------------------------------------------
# Entry Point
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Process a single RPS-BLAST ASN.1 file with rpsbproc.'
    )

    parser.add_argument(
        '--species',
        type=str,
        required=True,
        help='Species name to process.'
    )

    args = parser.parse_args()

    input_file  = defaults.PATH_DICT['ASN_RPSBLAST_DIR'] / f'{args.species}.asn'
    output_dir  = defaults.PATH_DICT['RPSBPROC_OUTPUT_DIR']
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f'{args.species}.txt'

    cmd: List[str] = [
        'rpsbproc',
        '-i', str(input_file),
        '-d', str(defaults.PATH_DICT['RPSBPROC_DB']),
        '-t', 'both',
        '-o', str(output_file)
    ]

    print(f'Processing {input_file} -> {output_file}')
    subprocess.run(cmd, check=True)
