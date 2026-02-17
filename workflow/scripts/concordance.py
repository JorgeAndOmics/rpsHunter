#!/usr/bin/env python3
"""
concordance.py — Multi-method concordance analysis for rpsHunter pipeline.

Quantifies agreement between HMMER (profile HMM) and RPSBLAST/CDD (PSSM)
domain calls for each final domain annotation. Reads existing pipeline outputs
and produces concordance-enriched tables without modifying upstream files.

Usage:
    python concordance.py \
        --aggregate results/tables/aggregate.parquet \
        --domains   results/tables/domains.parquet \
        --config    data/config/config.yaml \
        --outdir    results/concordance/
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def load_config(path: Path) -> tuple[list[str], list[str], bool]:
    """Load HMMER profiles, CDD targets, and HMMER enabled flag from config."""
    with open(path) as f:
        cfg = yaml.safe_load(f)

    hmmer_cfg = cfg.get('hmmer', {})
    hmmer_enabled = hmmer_cfg.get('enabled', False)
    hmmer_profiles = hmmer_cfg.get('profiles', []) if hmmer_enabled else []
    cdd_targets = cfg.get('rpsblast', {}).get('target_domains', [])

    return hmmer_profiles, cdd_targets, hmmer_enabled


# ---------------------------------------------------------------------------
# Domain normalization (prefix matching — same convention as defaults.py)
# ---------------------------------------------------------------------------

def normalize_domain(name: str, targets: list[str]) -> str | None:
    """Normalize a domain name to its config target via prefix matching.

    A target matches if name == target exactly, or name starts with target + '_'.
    Returns the first matching target, or None if no match.
    """
    if pd.isna(name):
        return None
    for target in targets:
        if name == target or name.startswith(f'{target}_'):
            return target
    return None


# ---------------------------------------------------------------------------
# Pre-computation of per-Tag HMMER family lookups
# ---------------------------------------------------------------------------

def _to_list(value):
    """Convert ndarray/list/None to a Python list."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return list(value)
    return []


def build_tag_hmm_lookups(
    aggregate_df: pd.DataFrame,
    hmmer_profiles: list[str],
) -> tuple[dict[str, set[str]], dict[str, dict[str, float]]]:
    """Pre-compute per-Tag HMM family sets and best E-values.

    Returns:
        tag_families: Tag → set of normalized HMMER family names
        tag_evalues:  Tag → {family: best_evalue}
    """
    tag_families: dict[str, set[str]] = {}
    tag_evalues: dict[str, dict[str, float]] = {}

    has_hmm = 'HMM_Domain' in aggregate_df.columns

    for _, row in aggregate_df[['Tag'] + (['HMM_Domain', 'HMM_Evalue'] if has_hmm else [])].iterrows():
        tag = row['Tag']

        if not has_hmm:
            tag_families[tag] = set()
            tag_evalues[tag] = {}
            continue

        domains = _to_list(row['HMM_Domain'])
        evalues = _to_list(row['HMM_Evalue'])

        families: set[str] = set()
        best: dict[str, float] = {}

        for d, e in zip(domains, evalues):
            norm = normalize_domain(d, hmmer_profiles)
            if norm is not None:
                families.add(norm)
                if norm not in best or e < best[norm]:
                    best[norm] = e

        tag_families[tag] = families
        tag_evalues[tag] = best

    return tag_families, tag_evalues


# ---------------------------------------------------------------------------
# Concordance logic
# ---------------------------------------------------------------------------

def build_concordance(
    aggregate_df: pd.DataFrame,
    domains_df: pd.DataFrame,
    hmmer_profiles: list[str],
    cdd_targets: list[str],
    hmmer_enabled: bool,
) -> pd.DataFrame:
    """Build per-domain concordance labels by joining domains to aggregate via Tag."""

    df = domains_df.copy()

    # Normalize CDD domain names to config targets
    df['Domain_Family'] = df['Domain'].apply(lambda d: normalize_domain(d, cdd_targets))
    # Domains not matching any config target keep their original name
    df['Domain_Family'] = df['Domain_Family'].fillna(df['Domain'])

    # Determine which families are HMMER-checkable
    hmmer_set = set(hmmer_profiles)
    if hmmer_enabled:
        df['HMMER_Checkable'] = df['Domain_Family'].isin(hmmer_set)
    else:
        df['HMMER_Checkable'] = False

    # Build BLAST-level columns from aggregate (one row per Tag)
    blast_cols = {'Tag': 'Tag'}
    col_map = {'E-value': 'Blast_Evalue', 'Bit Score': 'Blast_Bitscore', 'Pct Identity': 'Blast_Identity'}
    for src, dst in col_map.items():
        if src in aggregate_df.columns:
            blast_cols[src] = dst

    agg_blast = aggregate_df[list(blast_cols.keys())].drop_duplicates(subset='Tag')
    agg_blast = agg_blast.rename(columns={k: v for k, v in blast_cols.items() if k != v})
    df = df.merge(agg_blast, on='Tag', how='left')

    # Pre-compute per-Tag HMMER lookups (fast dict access per domain row)
    tag_families, tag_evalues = build_tag_hmm_lookups(aggregate_df, hmmer_profiles)

    # Assign concordance labels
    labels = np.empty(len(df), dtype=object)
    hmm_best = np.empty(len(df), dtype=object)

    for i, (_, row) in enumerate(df[['Tag', 'Domain_Family', 'HMMER_Checkable']].iterrows()):
        if not hmmer_enabled or not row['HMMER_Checkable']:
            labels[i] = 'hmmer_not_searched'
            hmm_best[i] = None
            continue

        tag = row['Tag']
        family = row['Domain_Family']
        hmm_fams = tag_families.get(tag, set())

        if family in hmm_fams:
            labels[i] = 'confirmed'
            hmm_best[i] = tag_evalues.get(tag, {}).get(family)
        else:
            labels[i] = 'hmmer_unmatched'
            hmm_best[i] = None

    df['Concordance'] = labels
    df['HMM_Best_Evalue'] = hmm_best

    return df


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def summarize_by_sequence(concordance_df: pd.DataFrame, aggregate_df: pd.DataFrame) -> pd.DataFrame:
    """Per-sequence (Tag) concordance summary."""

    has_hmm = 'HMM_Domain' in aggregate_df.columns

    # Pre-compute per-Tag unique HMM family count
    tag_hmm_count: dict[str, int] = {}
    if has_hmm:
        for _, row in aggregate_df[['Tag', 'HMM_Domain']].iterrows():
            tag_hmm_count[row['Tag']] = len(set(_to_list(row['HMM_Domain'])))

    records = []
    for tag, grp in concordance_df.groupby('Tag'):
        conc = grp['Concordance']
        n_confirmed = int((conc == 'confirmed').sum())
        n_unmatched = int((conc == 'hmmer_unmatched').sum())
        n_not_searched = int((conc == 'hmmer_not_searched').sum())

        denom = n_confirmed + n_unmatched
        conc_rate = n_confirmed / denom if denom > 0 else None

        records.append({
            'Tag': tag,
            'Species': grp['Species'].iloc[0],
            'Species_Name': grp.get('Species_Name', grp['Species']).iloc[0],
            'Query_Accession': grp['Query_Accession'].iloc[0] if 'Query_Accession' in grp.columns else '',
            'N_CDD_Domains': len(grp),
            'N_CDD_Families': grp['Domain_Family'].nunique(),
            'N_HMM_Families': tag_hmm_count.get(tag, 0),
            'N_Confirmed': n_confirmed,
            'N_Unmatched': n_unmatched,
            'N_Not_Searched': n_not_searched,
            'Concordance_Rate': conc_rate,
        })

    return pd.DataFrame(records)


def summarize_by_family(concordance_df: pd.DataFrame) -> pd.DataFrame:
    """Per-domain-family concordance summary."""

    records = []
    for family, grp in concordance_df.groupby('Domain_Family'):
        conc = grp['Concordance']
        n_confirmed = int((conc == 'confirmed').sum())
        n_unmatched = int((conc == 'hmmer_unmatched').sum())
        n_not_searched = int((conc == 'hmmer_not_searched').sum())

        denom = n_confirmed + n_unmatched
        conc_rate = n_confirmed / denom if denom > 0 else None

        evalues = pd.to_numeric(grp['Evalue'], errors='coerce')
        bitscores = pd.to_numeric(grp['Bitscore'], errors='coerce')

        confirmed_hmm = pd.to_numeric(
            grp.loc[conc == 'confirmed', 'HMM_Best_Evalue'], errors='coerce'
        )

        records.append({
            'Query_Accession': grp['Query_Accession'].iloc[0] if 'Query_Accession' in grp.columns else '',
            'Domain_Family': family,
            'Total': len(grp),
            'Confirmed': n_confirmed,
            'Unmatched': n_unmatched,
            'Not_Searched': n_not_searched,
            'Concordance_Rate': conc_rate,
            'Median_CDD_Evalue': evalues.median() if not evalues.isna().all() else None,
            'Median_CDD_Bitscore': bitscores.median() if not bitscores.isna().all() else None,
            'Median_HMM_Evalue': (
                confirmed_hmm.median()
                if len(confirmed_hmm) > 0 and not confirmed_hmm.isna().all()
                else None
            ),
        })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def write_outputs(df: pd.DataFrame, outdir: Path, name: str) -> None:
    """Write a DataFrame as both Parquet and CSV."""
    df.to_parquet(outdir / f'{name}.parquet', index=False)
    df.to_csv(outdir / f'{name}.csv', index=False)
    print(f'  {name}: {len(df):,} rows')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Multi-method concordance analysis for rpsHunter pipeline outputs.'
    )
    parser.add_argument('--aggregate', required=True, type=Path,
                        help='Path to aggregate.parquet (BLAST + HMMER enriched)')
    parser.add_argument('--domains', required=True, type=Path,
                        help='Path to domains.parquet (RPSBLAST/CDD domain calls)')
    parser.add_argument('--config', required=True, type=Path,
                        help='Path to pipeline config.yaml')
    parser.add_argument('--outdir', required=True, type=Path,
                        help='Output directory for concordance tables')
    args = parser.parse_args()

    # Validate inputs
    for path, label in [(args.aggregate, 'aggregate'), (args.domains, 'domains'), (args.config, 'config')]:
        if not path.exists():
            print(f'Error: {label} file not found: {path}', file=sys.stderr)
            sys.exit(1)

    # Load config
    print('Loading config...')
    hmmer_profiles, cdd_targets, hmmer_enabled = load_config(args.config)
    print(f'  HMMER enabled: {hmmer_enabled}')
    print(f'  HMMER profiles: {hmmer_profiles}')
    print(f'  CDD targets: {cdd_targets}')

    # Load data
    print('Loading aggregate.parquet...')
    aggregate_df = pd.read_parquet(args.aggregate)
    print(f'  {len(aggregate_df):,} sequences')

    print('Loading domains.parquet...')
    domains_df = pd.read_parquet(args.domains)
    print(f'  {len(domains_df):,} domain annotations')

    # Build concordance
    print('Building per-domain concordance...')
    concordance_df = build_concordance(
        aggregate_df, domains_df, hmmer_profiles, cdd_targets, hmmer_enabled
    )

    label_counts = concordance_df['Concordance'].value_counts()
    for label, count in label_counts.items():
        print(f'  {label}: {count:,}')

    # Summaries
    print('Summarizing by sequence...')
    seq_summary = summarize_by_sequence(concordance_df, aggregate_df)

    print('Summarizing by domain family...')
    family_summary = summarize_by_family(concordance_df)

    # Write outputs
    print('Writing outputs...')
    write_outputs(concordance_df, args.outdir, 'concordance_domains')
    write_outputs(seq_summary, args.outdir, 'concordance_sequences')
    write_outputs(family_summary, args.outdir, 'concordance_summary')

    # Print family summary
    print('\n=== Domain Family Concordance Summary ===')
    print(family_summary.to_string(index=False))
    print(f'\nDone. Outputs in {args.outdir}/')


if __name__ == '__main__':
    main()
