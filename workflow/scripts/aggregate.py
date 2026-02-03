"""
Aggregates per-species parquet files into a combined output.

Usage:
    python aggregate.py <type> <output_parquet> <input1> <input2> ...

Types: blast, rpsblast, domains
"""

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='Aggregate per-species parquet files.')
    parser.add_argument('type', choices=['blast', 'rpsblast', 'domains'],
                        help='Aggregation type (used for logging).')
    parser.add_argument('output_parquet', help='Output parquet path.')
    parser.add_argument('inputs', nargs='+', help='Input parquet file paths.')
    args = parser.parse_args()

    dfs = []
    for input_path in args.inputs:
        if Path(input_path).exists():
            dfs.append(pd.read_parquet(input_path))

    output_path = Path(args.output_parquet)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
        combined.to_parquet(output_path, index=False)

        csv_path = output_path.with_suffix('.csv')
        try:
            combined.to_csv(csv_path, index=False)
        except PermissionError:
            print(f'Warning: could not write {csv_path} (permission denied); parquet written successfully')

        print(f'Aggregated {len(dfs)} {args.type} files -> {len(combined)} rows')
    else:
        print(f'No input files found for {args.type} aggregation')


if __name__ == '__main__':
    main()
