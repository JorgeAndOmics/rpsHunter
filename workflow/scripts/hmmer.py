"""
HMMER domain filtering — single-species mode.

Reads a per-species blast parquet, extracts Subject Sequence, runs hmmsearch
against a Pfam HMM database, and enriches the parquet with HMM_* columns:
    HMM_Filtered   (bool)   – passed all HMMER quality filters
    HMM_Evalue     (float)  – best per-domain E-value
    HMM_Score      (float)  – best domain bit score
    HMM_Domain     (str)    – best matching Pfam domain name
    HMM_Coverage   (float)  – fraction of HMM model covered by best hit

Usage (invoked by Snakefile):
    python hmmer.py --species <name> <input.parquet> <Pfam-A.hmm> <output.parquet>
"""

import argparse
import logging
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

import defaults
from colored_logging import colored_logging


# -----------------------------------------------------------------------------
# FASTA helpers
# -----------------------------------------------------------------------------

def df_to_fasta(df: pd.DataFrame, fasta_path: Path) -> None:
    """Write Subject Sequence column as a FASTA file keyed by row index."""
    with open(fasta_path, 'w') as fh:
        for idx, row in df.iterrows():
            seq = str(row['Subject Sequence']).replace('-', '').replace('*', '')
            if seq:
                fh.write(f'>seq_{idx}\n{seq}\n')


# -----------------------------------------------------------------------------
# Profile extraction (hmmfetch subset)
# -----------------------------------------------------------------------------

def extract_profiles(full_hmm: Path, accessions: list, subset_hmm: Path) -> None:
    """Extract a subset of HMM profiles from a pressed database using hmmfetch.

    Accessions are bare (e.g. PF00855).  Pfam ACC fields are versioned
    (PF00855.24), so we grep to resolve each accession to its HMM NAME,
    then fetch by NAME (hmmfetch default).

    Presses full_hmm if not already pressed.  Writes the extracted profiles
    to subset_hmm and presses that file so hmmsearch can use it directly.
    """
    # Ensure the full database is pressed (one-time cost, ~30 s for Pfam).
    # hmmpress appends suffixes to the full filename (Pfam-A.hmm.h3m, etc.).
    # Check all four; if any is missing the index is incomplete.  Use -f so
    # a partially-written index from a previous interrupted run doesn't block.
    _index_exts = ('.h3f', '.h3i', '.h3m', '.h3p')
    if not all(Path(str(full_hmm) + ext).exists() for ext in _index_exts):
        logging.info('Pressing full HMM database (one-time)...')
        result = subprocess.run(['hmmpress', '-f', str(full_hmm)], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f'hmmpress failed: {result.stderr}')

    # Resolve bare accessions → HMM NAMEs via grep on the flat file.
    # Pattern matches "ACC   PF00855." (prefix + dot) so version is irrelevant.
    pattern = '|'.join(f'^ACC   {acc}\\.' for acc in accessions)
    result = subprocess.run(
        ['grep', '-B1', '-E', pattern, str(full_hmm)],
        capture_output=True, text=True
    )
    names = [line.split()[1] for line in result.stdout.split('\n') if line.startswith('NAME')]
    if not names:
        raise RuntimeError(f'No matching HMM profiles found for accessions {accessions}')
    logging.info(f'Resolved {accessions} → {names}')

    # Write NAME key file and fetch
    key_file = subset_hmm.with_suffix('.keys')
    key_file.write_text('\n'.join(names) + '\n')

    result = subprocess.run(
        ['hmmfetch', '-f', str(full_hmm), str(key_file)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f'hmmfetch failed: {result.stderr}')
    if not result.stdout.strip():
        raise RuntimeError(f'hmmfetch returned no profiles for {names}')
    subset_hmm.write_text(result.stdout)

    # Press the subset so hmmsearch can use it
    result = subprocess.run(['hmmpress', str(subset_hmm)], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f'hmmpress on subset failed: {result.stderr}')

    logging.info(f'Extracted {len(names)} profiles → {subset_hmm}')


# -----------------------------------------------------------------------------
# hmmsearch execution
# -----------------------------------------------------------------------------

def run_hmmsearch(fasta_path: Path, hmm_path: Path, tbl_path: Path, threads: int) -> None:
    """Run hmmsearch with config-driven flags; write --tblout to tbl_path."""
    cmd = ['hmmsearch']

    if defaults.HMMER_USE_GA:
        cmd.append('--cut_ga')
    else:
        cmd.extend(['--domE', str(defaults.HMMER_DOM_EVALUE)])
        cmd.extend(['-E', str(defaults.HMMER_EVALUE)])
        if defaults.HMMER_MIN_SCORE is not None:
            cmd.extend(['--score', str(defaults.HMMER_MIN_SCORE)])

    if defaults.HMMER_MAX_SENS:
        cmd.append('--max')

    if not defaults.HMMER_BIAS_FILTER:
        cmd.append('--nobias')

    cmd.extend(['--cpu', str(threads)])
    cmd.extend(['--seed', str(defaults.HMMER_SEED)])
    cmd.extend(['--tblout', str(tbl_path)])
    cmd.extend(['--domtblout', str(tbl_path.with_suffix('.domtbl'))])
    cmd.append(str(hmm_path))
    cmd.append(str(fasta_path))

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f'hmmsearch failed: {result.stderr}')
        raise RuntimeError('hmmsearch exited with non-zero status')


# -----------------------------------------------------------------------------
# Domain table parser
# -----------------------------------------------------------------------------

def parse_domtbl(domtbl_path: Path) -> pd.DataFrame:
    """
    Parse hmmsearch --domtblout into a DataFrame.

    Columns retained: target_name, target_len, query_name (Pfam accession),
    query_name_desc (Pfam name), from, to, evalue, score, cov_from, cov_to, cov_len.
    """
    rows = []
    with open(domtbl_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 23:
                continue
            rows.append({
                'target_name': parts[0],
                'target_len': int(parts[2]),
                'query_name': parts[3],       # Pfam domain NAME (e.g. zf-C2H2)
                'evalue': float(parts[12]),   # per-domain i-evalue
                'score': float(parts[13]),    # per-domain score
                'from': int(parts[17]),       # ali from (alignment start on target)
                'to': int(parts[18]),         # ali to   (alignment end on target)
                'cov_from': int(parts[15]),   # hmm from (model start)
                'cov_to': int(parts[16]),     # hmm to   (model end)
                'cov_len': int(parts[5]),     # qlen     (model length)
            })
            # Pfam domain name is everything after column 22
            if len(parts) > 22:
                rows[-1]['query_name_desc'] = ' '.join(parts[22:])
            else:
                rows[-1]['query_name_desc'] = ''
    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=['target_name', 'target_len', 'query_name', 'evalue',
                 'score', 'from', 'to', 'cov_from', 'cov_to', 'cov_len', 'query_name_desc'])


# -----------------------------------------------------------------------------
# Per-sequence best-hit selection + quality filtering
# -----------------------------------------------------------------------------

def best_hits(domtbl: pd.DataFrame) -> pd.DataFrame:
    """For each target sequence keep only the hit with lowest E-value."""
    if domtbl.empty:
        return domtbl
    return domtbl.sort_values('evalue').drop_duplicates(subset='target_name', keep='first')


def quality_filter(hits: pd.DataFrame) -> pd.DataFrame:
    """Apply min_coverage and min_alignment_length filters."""
    if hits.empty:
        return hits
    hits = hits.copy()
    hits['coverage'] = (hits['cov_to'] - hits['cov_from'] + 1) / hits['cov_len']
    hits['aln_len'] = hits['to'] - hits['from'] + 1
    return hits[
        (hits['coverage'] >= defaults.HMMER_MIN_COVERAGE) &
        (hits['aln_len'] >= defaults.HMMER_MIN_ALN_LEN)
    ]


# -----------------------------------------------------------------------------
# Enrichment
# -----------------------------------------------------------------------------

def enrich(blast_df: pd.DataFrame, hits: pd.DataFrame) -> pd.DataFrame:
    """
    Add HMM_* columns to blast_df.  hits is indexed by 'seq_{row_index}'.
    """
    # Build a lookup: seq_N → row fields
    lookup = {}
    for _, row in hits.iterrows():
        idx = int(row['target_name'].split('_')[1])
        lookup[idx] = row

    hmm_filtered = []
    hmm_evalue  = []
    hmm_score   = []
    hmm_domain  = []
    hmm_cov     = []

    for idx in range(len(blast_df)):
        if idx in lookup:
            h = lookup[idx]
            hmm_filtered.append(True)
            hmm_evalue.append(h['evalue'])
            hmm_score.append(h['score'])
            hmm_domain.append(h['query_name'])
            hmm_cov.append(h['coverage'])
        else:
            hmm_filtered.append(False)
            hmm_evalue.append(None)
            hmm_score.append(None)
            hmm_domain.append(None)
            hmm_cov.append(None)

    blast_df = blast_df.copy()
    blast_df['HMM_Filtered']  = hmm_filtered
    blast_df['HMM_Evalue']    = hmm_evalue
    blast_df['HMM_Score']     = hmm_score
    blast_df['HMM_Domain']    = hmm_domain
    blast_df['HMM_Coverage']  = hmm_cov
    return blast_df


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    colored_logging(log_file_name='hmmer.txt')

    parser = argparse.ArgumentParser(description='Run HMMER domain filtering for a single species.')
    parser.add_argument('--species', required=True, help='Species name.')
    parser.add_argument('--threads', type=int, default=1, help='CPU threads for hmmsearch (set by Snakemake).')
    parser.add_argument('input_parquet', help='Per-species blast parquet (input).')
    parser.add_argument('pfam_hmm', help='Path to Pfam-A.hmm profile database.')
    parser.add_argument('output_parquet', help='Per-species enriched parquet (output).')
    args = parser.parse_args()

    blast_df = pd.read_parquet(args.input_parquet)

    # Remove previous HMM columns if re-running
    hmm_cols = ['HMM_Filtered', 'HMM_Evalue', 'HMM_Score', 'HMM_Domain', 'HMM_Coverage']
    blast_df = blast_df.drop(columns=[c for c in hmm_cols if c in blast_df.columns])

    # Short-circuit: no searchable sequences
    if 'Subject Sequence' not in blast_df.columns:
        logging.warning(f'{args.species}: no Subject Sequence column — skipping hmmsearch')
        hits = pd.DataFrame()
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path  = Path(tmpdir) / 'seqs.fasta'
            tbl_path    = Path(tmpdir) / 'hits.tbl'

            df_to_fasta(blast_df, fasta_path)

            # Check if any non-empty sequences were written
            if fasta_path.stat().st_size == 0:
                logging.warning(f'{args.species}: all sequences empty after gap/stop stripping — skipping hmmsearch')
                hits = pd.DataFrame()
            else:
                # If specific profiles requested, extract subset first (seconds vs minutes)
                if defaults.HMMER_PROFILES:
                    subset_hmm = Path(tmpdir) / 'subset.hmm'
                    extract_profiles(Path(args.pfam_hmm), defaults.HMMER_PROFILES, subset_hmm)
                    hmm_to_search = subset_hmm
                else:
                    hmm_to_search = Path(args.pfam_hmm)

                run_hmmsearch(fasta_path, hmm_to_search, tbl_path, args.threads)

                domtbl_path = tbl_path.with_suffix('.domtbl')
                domtbl = parse_domtbl(domtbl_path)
                hits = quality_filter(best_hits(domtbl))

    blast_enriched = enrich(blast_df, hits)

    n_pass = blast_enriched['HMM_Filtered'].sum()
    logging.info(f'{args.species}: {n_pass}/{len(blast_enriched)} sequences pass HMMER filters')

    Path(args.output_parquet).parent.mkdir(parents=True, exist_ok=True)
    blast_enriched.to_parquet(args.output_parquet, index=False)


if __name__ == '__main__':
    main()
