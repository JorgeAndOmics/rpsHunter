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

## 2. A Species With No Results

**Symptom:** A species has empty tables, an empty FASTA or a GFF3 file with only its header, and the
pipeline finished without error.

**This is a result, not a failure.** No hit of that species passed the quality gate, or none
carried a domain. Each step writes a valid empty file with all its columns and the next step
accepts it. A tool that fails stops its job with an error instead, so an empty table always means
"searched, nothing found". The job's log under `logs/jobs/` shows what ran.

---

## 3. Windows/WSL File Locking

**Symptom:** `PermissionError` when writing CSV files, typically in `blast_parser.py`.

**Cause:** On Windows/WSL, CSV files fail to write if Excel or another process has them open. Windows file locks prevent overwriting files that are in use.

**Solution:** Close Excel (or any other application) that has the output CSV files open before rerunning the pipeline.

**Note:** The pipeline includes a `try/except` guard in `blast_parser.py` specifically for this case. The Parquet output file will still be written successfully even when the CSV write fails, so no data is lost. The CSV is a convenience copy; the Parquet is the authoritative output.

---

## 4. Pfam Accession Not Found

**Symptom:** `domains: Pfam accessions not found in the HMM library: PF.....`

**Cause:** A `pfam:` accession under `domains:` is not in the Pfam release installed under
`db_root_folder/hmm_dbs/`. Accessions are written without their version (`PF00096`, not
`PF00096.32`); a malformed one is refused earlier, when the config is read.

**Solution:** Look the family up at [InterPro](https://www.ebi.ac.uk/interpro/) and check that the
accession exists in the release you downloaded (`grep -c "^ACC   PF00096" Pfam-A.hmm`).

---

## 5. CDD Model Not Found

**Symptom:** Building the families' RPS-BLAST database fails, naming a missing `.smp` file.

**Cause:** A `cdd:` accession under `domains:` has no model file in `db_root_folder/cdd_smp/`:
it is misspelt, or it does not exist in the CDD release installed. CDD accessions look like
`pfam00096`, `cd19193`, `smart00355` or `COG5048`, and are case-sensitive.

**Solution:** Find the accession in column 2 of `db_root_folder/rpsbproc_dbs/cddid.tbl`, or search
the family at [NCBI CDD](https://www.ncbi.nlm.nih.gov/cdd). Use accessions, not short names: short
names change between releases.

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

## 9. Multi-Flag CLI Invocation Fails

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
