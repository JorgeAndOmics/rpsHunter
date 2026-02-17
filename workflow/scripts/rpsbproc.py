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
            defaults.RPSBPROC_CMD,
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
        '--input-asn',
        type=str,
        required=True,
        help='Path to the input ASN.1 file.'
    )

    parser.add_argument(
        '--output-txt',
        type=str,
        required=True,
        help='Path to the output text file.'
    )

    args = parser.parse_args()

    input_file = Path(args.input_asn)
    output_file = Path(args.output_txt)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Empty ASN (species had 0 sequences after gating) — rpsbproc would crash on it.
    # Touch the output so the Snakemake rule is satisfied; rpsbproc_parser already
    # handles empty .txt files and produces a 0-row parquet.
    if input_file.stat().st_size == 0:
        output_file.touch()
        print(f'Empty ASN — skipping rpsbproc, touched {output_file}')
    else:
        cmd: List[str] = [
            defaults.RPSBPROC_CMD,
            '-i', str(input_file),
            '-d', str(defaults.PATH_DICT['RPSBPROC_DB']),
            '-t', 'both',
            '-o', str(output_file)
        ]

        print(f'Processing {input_file} -> {output_file}')
        subprocess.run(cmd, check=True)
