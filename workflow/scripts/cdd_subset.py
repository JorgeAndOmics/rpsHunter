"""
CDD Subset Database Builder
============================

Resolves target domain names from config against cddid.tbl, collects matching
SMP files, and builds a subset RPSBLAST database using makeprofiledb.

Usage (invoked by Snakefile):
    python cdd_subset.py --cddid <path> --smp-dir <path> --output-db <path>
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import defaults
from colored_logging import colored_logging


def main() -> None:
    colored_logging(log_file_name='cdd_subset.txt')

    parser = argparse.ArgumentParser(description='Build CDD subset database for RPSBLAST.')
    parser.add_argument('--cddid', type=str, required=True,
                        help='Path to cddid.tbl.')
    parser.add_argument('--smp-dir', type=str, required=True,
                        help='Directory containing extracted SMP files.')
    parser.add_argument('--output-db', type=str, required=True,
                        help='Output database path prefix (without extension).')
    args = parser.parse_args()

    target_domains = defaults.RPSBLAST_TARGET_DOMAINS
    if not target_domains:
        logging.error('No target_domains configured — cannot build subset database.')
        sys.exit(1)

    cddid_path = Path(args.cddid)
    smp_dir = Path(args.smp_dir)
    output_db = Path(args.output_db)
    output_db.parent.mkdir(parents=True, exist_ok=True)

    # Resolve target domain names → CDD accessions
    accessions = defaults.resolve_cdd_targets(target_domains, cddid_path)
    logging.info(f'Resolved {len(accessions)} CDD accessions from {len(target_domains)} target domains.')

    # Collect SMP files
    found_smps = []
    missing_smps = []
    for acc in accessions:
        smp_file = smp_dir / f'{acc}.smp'
        if smp_file.exists():
            found_smps.append(str(smp_file))
        else:
            missing_smps.append(acc)

    if missing_smps:
        logging.warning(f'{len(missing_smps)} SMP files not found: {missing_smps[:10]}...')

    if not found_smps:
        logging.error('No SMP files found for target domains — cannot build subset database.')
        sys.exit(1)

    # Write SMP list file for makeprofiledb
    smp_list_path = output_db.parent / 'smp_list.txt'
    smp_list_path.write_text('\n'.join(found_smps) + '\n')
    logging.info(f'Collected {len(found_smps)} SMP files → {smp_list_path}')

    # Build subset database
    logging.info(f'Running {defaults.MAKEPROFILEDB_CMD}...')
    result = subprocess.run(
        [defaults.MAKEPROFILEDB_CMD,
         '-in', str(smp_list_path),
         '-out', str(output_db),
         '-title', 'rpsHunter target domains'],
        capture_output=True, text=True
    )

    if result.returncode != 0:
        logging.error(f'makeprofiledb failed: {result.stderr}')
        sys.exit(1)

    logging.info(f'Subset database built at {output_db} ({len(found_smps)} PSSMs)')


if __name__ == '__main__':
    main()
