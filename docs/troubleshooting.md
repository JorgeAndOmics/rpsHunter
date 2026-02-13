# Troubleshooting Guide

This document covers common issues encountered when running the rpsHunter pipeline, along with their solutions and workarounds.

---

## 1. Never Use `snakemake --forceall`

**Symptom:** Pipeline attempts to re-download all genomes, rebuild all BLAST databases, and rerun every rule from scratch.

**Why this is dangerous:** The `--forceall` flag forces rerun of ALL rules, including genome download. This will:

- Delete existing genome FASTA files (`.fna` files in the species database directory)
- Attempt to re-download approximately 6 GB of genome data from NCBI
- Recreate BLAST databases unnecessarily
- Take hours to complete

**Safe alternatives:**

1. **Delete specific output files** and rerun normally. Snakemake will only regenerate the missing files:
   ```bash
   rm -rf results/fastas/*.fa results/rpsblast/*.parquet results/domains/*.parquet
   conda run -n rpsHunter ./rpsHunter --blast --skip-validation
   ```

2. **Use `--forcerun` with specific (non-wildcarded) rule names** to rerun a particular stage without touching upstream outputs.

3. **Touch upstream files** to trigger dependency-based reruns (use with caution).

See also: [Configuration](configuration.md) for details on config-driven reruns.

---

## 2. Empty Species (0 Sequences After Filtering)

**Symptom:** A species produces no sequences after quality, ORF, or HMMER filtering, but the pipeline continues without error.

**This is expected behavior, not a bug.** Some species genuinely produce zero hits after filtering. The pipeline handles this gracefully through the entire downstream chain:

1. Empty FASTA is written (0 sequences)
2. RPSBLAST writes an empty Parquet with the correct 15-column schema and touches the ASN file
3. rpsbproc detects the 0-byte ASN, touches the output `.txt`, and exits early
4. rpsbproc_parser produces a 0-row Parquet from the empty `.txt`
5. Aggregate ignores empty species cleanly via `pd.concat`

**No intervention is needed.** The empty-species path has been tested end-to-end across multiple pipeline runs.

---

## 3. Windows/WSL File Locking

**Symptom:** `PermissionError` when writing CSV files, typically in `blast_parser.py`.

**Cause:** On Windows/WSL, CSV files fail to write if Excel or another process has them open. Windows file locks prevent overwriting files that are in use.

**Solution:** Close Excel (or any other application) that has the output CSV files open before rerunning the pipeline.

**Note:** The pipeline includes a `try/except` guard in `blast_parser.py` specifically for this case. The Parquet output file will still be written successfully even when the CSV write fails, so no data is lost. The CSV is a convenience copy; the Parquet is the authoritative output.

---

## 4. HMMER Profile Not Found

**Symptom:** Error message: `Domain name(s) not found in Pfam-A.hmm`

**Cause:** The `profiles` list in `config.yaml` contains a name that does not match any entry in the Pfam HMM database.

**Solution:** Verify the exact Pfam **NAME** identifier for each profile. Key rules:

- Use the Pfam **NAME**, not the accession (e.g., `zf-C2H2`, not `PF00096`) and not the description.
- Names are **case-sensitive**: use `zf-C2H2`, not `ZF-C2H2` or `Zf-c2h2`.
- Look up correct names at [InterPro](https://www.ebi.ac.uk/interpro/).

**Example of correct config:**
```yaml
hmmer:
  profiles:
    - 'KRAB'
    - 'SET'
    - 'SSXRD'
    - 'zf-C2H2'
```

---

## 5. CDD Domain Not in cddid.tbl

**Symptom:** A domain specified in `rpsblast.target_domains` is not found when building the CDD subset database.

**Cause:** The `target_domains` list uses CDD **ShortName** identifiers from `cddid.tbl`, which are different from Pfam names. The two databases use different naming conventions:

| Database | Example Name |
|----------|-------------|
| Pfam     | `zf-C2H2`  |
| CDD      | `ZnF_C2H2` |

**Solution:**

1. Check the `cddid.tbl` file in `db_root_folder/rpsbproc_dbs/` for the exact ShortName values.
2. Use the ShortName column (column 2) from that file.
3. Prefix matching is supported: specifying `KRAB` will match both `KRAB` and `KRAB_A-box`.

**Example:**
```yaml
rpsblast:
  target_domains:
    - 'KRAB'
    - 'SET'
    - 'SSXRD'
    - 'ZnF_C2H2'
    - 'COG5048'
```

---

## 6. BLAST+ Program Naming

**Symptom:** `FileNotFoundError` or `command not found` when running BLAST+ tools. Typically occurs because the system has `rpsblast+` instead of `rpsblast`.

**Cause:** Different BLAST+ installations use different binary names. The bioconda package uses `rpsblast`, but some system installations use `rpsblast+`.

**Solution:** Use the optional `programs:` section in `config.yaml` to override program names:

```yaml
programs:
  rpsblast: 'rpsblast+'
  tblastn: 'tblastn'
  blast_formatter: 'blast_formatter'
  makeblastdb: 'makeblastdb'
  rpsbproc: 'rpsbproc'
```

The pipeline defaults match bioconda package naming (`rpsblast`, `tblastn`, etc.). Override only the programs that differ on your system.

See also: [Configuration](configuration.md) for the full `programs` schema.

---

## 7. Pipeline Resumption After Interruption

**Symptom:** Pipeline was interrupted (Ctrl+C, crash, network failure) and you want to resume.

**Solution:** Simply re-run the same command. Snakemake automatically resumes from the last successful checkpoint. Completed rules will not be re-executed.

If a job was partially completed (e.g., a file was being written when the process was killed), the `--rerun-incomplete` flag handles this. This flag is already included in rpsHunter's default Snakemake invocation, so no additional action is needed.

```bash
# Just rerun the same command:
conda run -n rpsHunter ./rpsHunter --rpsblast --skip-validation
```

---

## 8. Large Genome Downloads and NCBI Rate Limits

**Symptom:** Genome downloads are extremely slow or fail with HTTP 429 (Too Many Requests) errors.

**Cause:** NCBI rate-limits unauthenticated requests to 3 requests per second.

**Solution:**

1. Set an NCBI API key via the validation prompt when running the pipeline, or set it as an environment variable:
   ```bash
   export NCBI_API_KEY="your_api_key_here"
   ```
   This increases the rate limit from 3 to 10 requests per second.

2. For very large genomes, downloads may take hours even with an API key. Plan accordingly and ensure a stable network connection.

3. If a download is interrupted, re-running the same command will resume from where it left off (see [Pipeline Resumption](#7-pipeline-resumption-after-interruption)).

---

## 9. Chromosome Name Mismatches in contingency_sorter

**Symptom:** Entire species missing from the contingency table and 3D scatter plot, despite having data in `domains.parquet`.

**Cause (fixed):** A restrictive regex (`[A-Za-z]+_?[0-9]+[._]?[0-9]*`) was applied to `Chromosome` values via `str_extract`, which silently dropped identifiers with multiple underscores (e.g., `manual_scaffold_10`), mixed letter-digit-letter patterns (e.g., `HAP1_SUPER_5`), or other non-standard scaffold names. If all rows for a species were rejected, the species vanished from downstream outputs.

**Fix:** The regex extraction was removed. Chromosome identifiers are already parsed upstream by `rpsbproc_parser.py`, so the R script now only filters out rows with empty or NA chromosome values.

---

## 10. ORF Filter Blocking Target Domains

**Symptom:** Target domains (e.g., KRAB, SET, SSXRD) are detected by HMMER but are missing from the final domain output.

**Cause:** Short protein domains may not contain open reading frames that meet the minimum ORF length threshold. The ORF filter in `blast_parser.py` blocks these sequences before they reach RPSBLAST. Domain-specific examples:

| Domain | Typical Size | ORF Detection at 200nt |
|--------|-------------|----------------------|
| KRAB   | ~65 aa      | Often fails          |
| SET    | ~130 aa     | Often fails          |
| SSXRD  | ~25 aa      | Always fails         |
| zf-C2H2 | ~28 aa    | Passes (tandem repeats extend total length) |

**Solution:** Disable ORF filtering in `config.yaml`:

```yaml
orf:
  enabled: false
```

This was the solution adopted in production for the KRAB/SET/SSXRD capture use case. Disabling the ORF filter resulted in a 2x increase in total domain annotations (10,289 to 20,973) and successfully captured all three previously-blocked domain types.

**Alternative:** Lower the ORF threshold (`min_orf_length: 50` or lower) instead of disabling the filter entirely. However, very short ORF thresholds provide minimal biological filtering value.

---

## 11. Multi-Flag CLI Invocation Fails

**Symptom:** Running `./rpsHunter --blast --rpsblast` fails with an error about unknown arguments.

**Cause:** Each CLI flag dispatches a **separate** Snakemake subprocess. When multiple flags are provided, subsequent flags are passed as unknown arguments to the first Snakemake invocation, which fails.

**Solution:** Run one flag per command, in the correct order:

```bash
conda run -n rpsHunter ./rpsHunter --blast --skip-validation
conda run -n rpsHunter ./rpsHunter --rpsblast --skip-validation
conda run -n rpsHunter ./rpsHunter --rpsbproc --skip-validation
conda run -n rpsHunter ./rpsHunter --rpsbproc-parser --skip-validation
conda run -n rpsHunter ./rpsHunter --completeness-detector --skip-validation
conda run -n rpsHunter ./rpsHunter --contingency-parser --skip-validation
```

The `--skip-validation` flag is the sole exception: it is consumed by the Python CLI dispatcher before Snakemake is invoked and can be combined with any single stage flag.

---

## 12. rpsbproc_parser Writes 0-Column Parquets

**Symptom:** Parquet files in `results/domains/` for some species have 0 rows and 0 columns.

**Cause:** When rpsbproc produces empty `.txt` output (a touched 0-byte file from the empty-species path), the parser writes a 0-row, 0-column Parquet because there is no data to infer a schema from.

**This is benign.** The `pd.concat` call in `aggregate.py` pulls the column schema from non-empty Parquet files contributed by other species. The 0-column files are silently ignored during aggregation. If any per-species files are missing entirely, `aggregate.py` now prints a warning listing the affected species.

**No action is needed** unless all species produce empty results, in which case the aggregate table will also be empty (which is correct behavior for a dataset with no domain annotations).
