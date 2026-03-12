#!/usr/bin/env python3
"""
hit_domain_inventory.py — Hit–Domain Inventory for rpsHunter pipeline.

Produces two tables:

1. **hit_domain_inventory** — One row per BLAST hit (Tag), enriched with all CDD
   domains found within that hit. Per-query; multiple inputs accepted for cross-query
   merging.

2. **locus_domain_inventory** — Overlapping BLAST hits from all queries on the same
   Species × Chromosome are clustered into genomic loci. One row per locus, with the
   union of all domains across all contributing hits and queries. Completeness flags
   indicate whether the locus contains the full expected domain set.

Usage (per-query):
    python hit_domain_inventory.py \
        --aggregate results/tables/{ql}/aggregate.parquet \
        --domains   results/tables/{ql}/domains.parquet \
        --config    data/config/config.yaml \
        --outdir    results/tables/{ql}/

Usage (merged / cross-query):
    python hit_domain_inventory.py \
        --aggregate results/tables/Q1/aggregate.parquet \
                    results/tables/Q2/aggregate.parquet \
        --domains   results/tables/Q1/domains.parquet \
                    results/tables/Q2/domains.parquet \
        --config    data/config/config.yaml \
        --outdir    results/tables/
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def load_config(path: Path) -> tuple[list[str], list[str], bool]:
    """Load HMMER profiles, CDD target domains, and HMMER enabled flag."""
    with open(path) as f:
        cfg = yaml.safe_load(f)

    hmmer_cfg = cfg.get('hmmer', {})
    hmmer_enabled = hmmer_cfg.get('enabled', False)
    hmmer_profiles = hmmer_cfg.get('profiles', []) if hmmer_enabled else []
    cdd_targets = cfg.get('rpsblast', {}).get('target_domains', [])

    return hmmer_profiles, cdd_targets, hmmer_enabled


# ---------------------------------------------------------------------------
# Domain matching (prefix convention — same as defaults.py / concordance.py)
# ---------------------------------------------------------------------------

def _matches_target(domain_name: str, target: str) -> bool:
    """Check if a domain name matches a config target (exact or prefix + '_')."""
    return domain_name == target or domain_name.startswith(f'{target}_')


def _find_matching_targets(domain_name: str, targets: list[str]) -> list[str]:
    """Return all config targets that a domain name matches."""
    if pd.isna(domain_name):
        return []
    return [t for t in targets if _matches_target(domain_name, t)]


# ---------------------------------------------------------------------------
# HMM domain helpers
# ---------------------------------------------------------------------------

def _to_list(value) -> list:
    """Convert ndarray/list/None to a Python list."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return list(value)
    return []


def _join_sorted(items) -> str:
    """Join sorted unique items with '; ' separator."""
    return '; '.join(sorted(set(items))) if items else ''


# ---------------------------------------------------------------------------
# Per-hit inventory
# ---------------------------------------------------------------------------

def build_inventory(
    aggregate_df: pd.DataFrame,
    domains_df: pd.DataFrame,
    hmmer_profiles: list[str],
    cdd_targets: list[str],
    hmmer_enabled: bool,
) -> pd.DataFrame:
    """Build one-row-per-BLAST-hit inventory with domain information."""

    has_hmm = 'HMM_Domain' in aggregate_df.columns and hmmer_enabled

    # --- Group domains by Tag ---
    domain_groups = {}
    for tag, grp in domains_df.groupby('Tag'):
        raw_domains = sorted(grp['Domain'].dropna().unique())

        cdd_present = set()
        for d in raw_domains:
            cdd_present.update(_find_matching_targets(d, cdd_targets))
        cdd_present = sorted(cdd_present)
        cdd_missing = sorted(set(cdd_targets) - set(cdd_present))

        detail = []
        for _, row in grp.iterrows():
            detail.append({
                'Domain': row.get('Domain'),
                'Evalue': row.get('Evalue'),
                'Bitscore': row.get('Bitscore'),
                'Incomplete': row.get('Incomplete'),
                'Hit_type': row.get('Hit_type'),
                'From': row.get('From'),
                'To': row.get('To'),
            })

        domain_groups[tag] = {
            'CDD_Domains_Found': raw_domains,
            'N_CDD_Domains': len(grp),
            'N_CDD_Families': len(raw_domains),
            'CDD_Targets_Present': cdd_present,
            'CDD_Targets_Missing': cdd_missing,
            'CDD_Domain_Detail': detail,
        }

    # --- Build per-hit records ---
    records = []
    for _, hit in aggregate_df.iterrows():
        tag = hit['Tag']
        dom_info = domain_groups.get(tag, {})

        rec = {
            'Tag': tag,
            'Species': hit.get('Species'),
            'Species_Name': hit.get('Species_Name'),
            'Query_Accession': hit.get('Query_Accession'),
            'Subject_ID': hit.get('Subject ID'),
            'S_Start': hit.get('S. Start'),
            'S_End': hit.get('S. End'),
            'E_value': hit.get('E-value'),
            'Bit_Score': hit.get('Bit Score'),
            'Pct_Identity': hit.get('Pct Identity'),
            'Alignment_Length': hit.get('Alignment Length'),
            'HMMER_Domains_Found': [],
            'HMMER_Profiles_Present': [],
            'HMMER_Profiles_Missing': hmmer_profiles.copy(),
            'Complete_HMMER': False,
            'CDD_Domains_Found': dom_info.get('CDD_Domains_Found', []),
            'N_CDD_Domains': dom_info.get('N_CDD_Domains', 0),
            'N_CDD_Families': dom_info.get('N_CDD_Families', 0),
            'CDD_Targets_Present': dom_info.get('CDD_Targets_Present', []),
            'CDD_Targets_Missing': dom_info.get('CDD_Targets_Missing', cdd_targets.copy()),
            'CDD_Domain_Detail': dom_info.get('CDD_Domain_Detail', []),
            'Complete_CDD': False,
        }

        if has_hmm:
            hmm_domains = _to_list(hit.get('HMM_Domain'))
            rec['HMMER_Domains_Found'] = sorted(set(hmm_domains)) if hmm_domains else []
            hmm_present = set()
            for d in hmm_domains:
                hmm_present.update(_find_matching_targets(d, hmmer_profiles))
            rec['HMMER_Profiles_Present'] = sorted(hmm_present)
            rec['HMMER_Profiles_Missing'] = sorted(set(hmmer_profiles) - hmm_present)
            rec['Complete_HMMER'] = len(hmm_present) == len(hmmer_profiles) and len(hmmer_profiles) > 0

        cdd_covers_hmmer = set()
        for d in rec['CDD_Domains_Found']:
            cdd_covers_hmmer.update(_find_matching_targets(d, hmmer_profiles))
        rec['Complete_CDD'] = (
            len(cdd_covers_hmmer) == len(hmmer_profiles) and len(hmmer_profiles) > 0
        )

        records.append(rec)

    df = pd.DataFrame(records)

    list_cols = [
        'HMMER_Domains_Found', 'HMMER_Profiles_Present', 'HMMER_Profiles_Missing',
        'CDD_Domains_Found', 'CDD_Targets_Present', 'CDD_Targets_Missing',
    ]
    for col in list_cols:
        if col in df.columns:
            df[col] = df[col].apply(lambda x: '; '.join(x) if x else '')

    if 'CDD_Domain_Detail' in df.columns:
        df['CDD_Domain_Detail'] = df['CDD_Domain_Detail'].apply(
            lambda x: json.dumps(x) if x else '[]'
        )

    return df


# ---------------------------------------------------------------------------
# Locus-level aggregation
# ---------------------------------------------------------------------------

def _cluster_hits(hits_df: pd.DataFrame, gap: int) -> list[list[int]]:
    """Cluster hits by overlapping/proximal genomic coordinates.

    Parameters
    ----------
    hits_df : DataFrame
        Must contain 'coord_start' and 'coord_end' columns (orientation-normalised).
    gap : int
        Maximum gap (bp) between two hits to merge into the same locus.

    Returns
    -------
    list of lists of row indices forming each cluster.
    """
    if hits_df.empty:
        return []

    sorted_df = hits_df.sort_values('coord_start')
    indices = sorted_df.index.tolist()

    clusters = []
    current = [indices[0]]
    current_end = sorted_df.loc[indices[0], 'coord_end']

    for idx in indices[1:]:
        start = sorted_df.loc[idx, 'coord_start']
        end = sorted_df.loc[idx, 'coord_end']
        if start <= current_end + gap:
            current.append(idx)
            current_end = max(current_end, end)
        else:
            clusters.append(current)
            current = [idx]
            current_end = end
    clusters.append(current)
    return clusters


def build_locus_inventory(
    aggregate_df: pd.DataFrame,
    domains_df: pd.DataFrame,
    hmmer_profiles: list[str],
    cdd_targets: list[str],
    hmmer_enabled: bool,
    gap: int = 50000,
) -> pd.DataFrame:
    """Build one-row-per-locus inventory by clustering overlapping BLAST hits.

    Parameters
    ----------
    gap : int
        Maximum inter-hit gap (bp) to bridge into the same locus. Default 50 kb
        accommodates typical mammalian introns between PRDM9-like exons.
    """

    has_hmm = 'HMM_Domain' in aggregate_df.columns and hmmer_enabled

    # Pre-compute per-Tag domain info
    tag_cdd_domains: dict[str, list[str]] = {}
    for tag, grp in domains_df.groupby('Tag'):
        tag_cdd_domains[tag] = sorted(grp['Domain'].dropna().unique())

    # Normalise coordinates
    agg = aggregate_df.copy()
    agg['coord_start'] = agg[['S. Start', 'S. End']].min(axis=1)
    agg['coord_end'] = agg[['S. Start', 'S. End']].max(axis=1)

    records = []
    locus_id = 0

    for (species, chrom), group in agg.groupby(['Species', 'Subject ID']):
        clusters = _cluster_hits(group, gap)

        for member_indices in clusters:
            members = group.loc[member_indices]
            locus_id += 1

            locus_start = int(members['coord_start'].min())
            locus_end = int(members['coord_end'].max())
            n_hits = len(members)
            query_accessions = sorted(members['Query_Accession'].dropna().unique())
            tags = members['Tag'].tolist()

            # Best BLAST hit (by bit score)
            best_idx = members['Bit Score'].idxmax()
            best_hit = members.loc[best_idx]

            # Union CDD domains across all hits in locus
            all_cdd_domains = set()
            for tag in tags:
                all_cdd_domains.update(tag_cdd_domains.get(tag, []))
            all_cdd_domains = sorted(all_cdd_domains)

            # CDD target coverage
            cdd_present = set()
            for d in all_cdd_domains:
                cdd_present.update(_find_matching_targets(d, cdd_targets))
            cdd_present = sorted(cdd_present)
            cdd_missing = sorted(set(cdd_targets) - set(cdd_present))

            # CDD completeness against HMMER expected set
            cdd_covers_hmmer = set()
            for d in all_cdd_domains:
                cdd_covers_hmmer.update(_find_matching_targets(d, hmmer_profiles))
            complete_cdd = (
                len(cdd_covers_hmmer) == len(hmmer_profiles) and len(hmmer_profiles) > 0
            )

            # Union HMMER domains across all hits in locus
            all_hmm_domains: set[str] = set()
            if has_hmm:
                for _, hit in members.iterrows():
                    all_hmm_domains.update(_to_list(hit.get('HMM_Domain')))

            hmm_present = set()
            for d in all_hmm_domains:
                hmm_present.update(_find_matching_targets(d, hmmer_profiles))
            complete_hmmer = (
                len(hmm_present) == len(hmmer_profiles) and len(hmmer_profiles) > 0
            )

            rec = {
                'Locus_ID': locus_id,
                'Species': species,
                'Species_Name': members['Species_Name'].iloc[0],
                'Subject_ID': chrom,
                'Locus_Start': locus_start,
                'Locus_End': locus_end,
                'Locus_Length': locus_end - locus_start + 1,
                'N_Hits': n_hits,
                'Query_Accessions': _join_sorted(query_accessions),
                'N_Queries': len(query_accessions),
                'Tags': '; '.join(tags),
                'Best_E_value': best_hit.get('E-value'),
                'Best_Bit_Score': best_hit.get('Bit Score'),
                'Best_Pct_Identity': best_hit.get('Pct Identity'),
                'HMMER_Domains_Found': _join_sorted(all_hmm_domains),
                'HMMER_Profiles_Present': _join_sorted(hmm_present),
                'HMMER_Profiles_Missing': _join_sorted(set(hmmer_profiles) - hmm_present),
                'Complete_HMMER': complete_hmmer,
                'CDD_Domains_Found': _join_sorted(all_cdd_domains),
                'N_CDD_Domains': sum(len(tag_cdd_domains.get(t, [])) for t in tags),
                'N_CDD_Families': len(all_cdd_domains),
                'CDD_Targets_Present': _join_sorted(cdd_present),
                'CDD_Targets_Missing': _join_sorted(cdd_missing),
                'Complete_CDD': complete_cdd,
            }
            records.append(rec)

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
        description='Hit–Domain Inventory: join BLAST hits with their CDD domain annotations.'
    )
    parser.add_argument('--aggregate', required=True, type=Path, nargs='+',
                        help='Path(s) to aggregate.parquet (BLAST hits)')
    parser.add_argument('--domains', required=True, type=Path, nargs='+',
                        help='Path(s) to domains.parquet (CDD domain annotations)')
    parser.add_argument('--config', required=True, type=Path,
                        help='Path to pipeline config.yaml')
    parser.add_argument('--outdir', required=True, type=Path,
                        help='Output directory')
    parser.add_argument('--locus-gap', type=int, default=50000,
                        help='Max gap (bp) between BLAST hits to cluster into one locus (default: 50000)')
    args = parser.parse_args()

    # Validate inputs
    for path in args.aggregate + args.domains + [args.config]:
        if not path.exists():
            print(f'Error: file not found: {path}', file=sys.stderr)
            sys.exit(1)

    # Load config
    print('Loading config...')
    hmmer_profiles, cdd_targets, hmmer_enabled = load_config(args.config)
    print(f'  HMMER enabled: {hmmer_enabled}')
    print(f'  HMMER profiles (expected domain set): {hmmer_profiles}')
    print(f'  CDD target domains: {cdd_targets}')

    # Load data (concatenate if multiple inputs)
    print('Loading aggregate parquets...')
    aggregate_df = pd.concat(
        [pd.read_parquet(p) for p in args.aggregate], ignore_index=True
    )
    print(f'  {len(aggregate_df):,} BLAST hits from {len(args.aggregate)} file(s)')

    print('Loading domain parquets...')
    domains_df = pd.concat(
        [pd.read_parquet(p) for p in args.domains], ignore_index=True
    )
    print(f'  {len(domains_df):,} domain annotations from {len(args.domains)} file(s)')

    # Build per-hit inventory
    print('Building hit–domain inventory...')
    inventory_df = build_inventory(
        aggregate_df, domains_df, hmmer_profiles, cdd_targets, hmmer_enabled
    )

    n_with_domains = (inventory_df['N_CDD_Domains'] > 0).sum()
    n_complete_cdd = inventory_df['Complete_CDD'].sum()
    n_complete_hmmer = inventory_df['Complete_HMMER'].sum()
    print(f'  Hits with CDD domains: {n_with_domains:,} / {len(inventory_df):,}')
    print(f'  Complete domain set (CDD): {n_complete_cdd:,}')
    print(f'  Complete domain set (HMMER): {n_complete_hmmer:,}')

    # Build locus-level inventory
    print(f'Building locus–domain inventory (gap={args.locus_gap:,} bp)...')
    locus_df = build_locus_inventory(
        aggregate_df, domains_df, hmmer_profiles, cdd_targets, hmmer_enabled,
        gap=args.locus_gap,
    )

    n_loci = len(locus_df)
    n_loci_with_domains = (locus_df['N_CDD_Families'] > 0).sum() if n_loci > 0 else 0
    n_loci_complete_cdd = locus_df['Complete_CDD'].sum() if n_loci > 0 else 0
    n_loci_complete_hmmer = locus_df['Complete_HMMER'].sum() if n_loci > 0 else 0
    n_multi_query = (locus_df['N_Queries'] > 1).sum() if n_loci > 0 else 0
    print(f'  Loci: {n_loci:,} ({n_multi_query:,} multi-query)')
    print(f'  Loci with CDD domains: {n_loci_with_domains:,}')
    print(f'  Complete domain set (CDD): {n_loci_complete_cdd:,}')
    print(f'  Complete domain set (HMMER): {n_loci_complete_hmmer:,}')

    # Write outputs
    print('Writing outputs...')
    write_outputs(inventory_df, args.outdir, 'hit_domain_inventory')
    write_outputs(locus_df, args.outdir, 'locus_domain_inventory')

    # Print locus species-level summary
    if n_loci > 0:
        print('\n=== Per-Species Locus Summary ===')
        species_summary = locus_df.groupby('Species_Name').agg(
            Total_Loci=('Locus_ID', 'count'),
            Multi_Query_Loci=('N_Queries', lambda x: (x > 1).sum()),
            Loci_With_Domains=('N_CDD_Families', lambda x: (x > 0).sum()),
            Complete_CDD=('Complete_CDD', 'sum'),
            Complete_HMMER=('Complete_HMMER', 'sum'),
        ).reset_index()
        print(species_summary.to_string(index=False))

    print(f'\nDone. Outputs in {args.outdir}/')


if __name__ == '__main__':
    main()
