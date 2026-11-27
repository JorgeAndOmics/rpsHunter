# Troubleshooting Guide

This document covers common issues encountered when running the rpsHunter pipeline, along with their solutions and workarounds.

---

## 1. Redoing Part of a Run

`--forceall` (or `-F`) is refused by the launcher: it would rebuild every database and download the genomes
again. To redo a step:

- change the config value that should take effect and run again: the steps it affects are redone
  (see "What reruns when" in the [CLI reference](cli-reference.md));
- or delete the outputs you want remade and run the stage again;
- `./rpsHunter <stage> -n` shows what would run before anything does.

---

## 2. A Species With No Results

**Symptom:** A species has empty tables, an empty FASTA or a GFF3 file with only its header, and the
pipeline finished without error.

**This is a result, not a failure.** No hit of that species passed the quality gate, or none
carried a domain. Each step writes a valid empty file with all its columns and the next step
accepts it. A tool that fails stops its job with an error instead, so an empty table always means
"searched, nothing found". The job's log under `logs/jobs/` shows what ran.

---

## 3. A CSV Cannot Be Written

**Symptom:** A job fails with `PermissionError` on a `.csv` file.

**Cause:** Another program has the file open (a spreadsheet, typically, on Windows through WSL).

**Solution:** Close it and run the same command again. Tables are written as Parquet first and
CSV second; the job is redone as a whole.

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

**Symptom:** The checks stop the run with `CDD model missing for <accession>: .../cdd_smp/<accession>.smp`.

**Cause:** A `cdd:` accession under `domains:` has no model file in `db_root_folder/cdd_smp/`:
it is misspelt, or it does not exist in the CDD release installed. CDD accessions look like
`pfam00096`, `cd19193`, `smart00355` or `COG5048`, and are case-sensitive.

**Solution:** Find the accession in column 2 of `db_root_folder/rpsbproc_dbs/cddid.tbl`, or search
the family at [NCBI CDD](https://www.ncbi.nlm.nih.gov/cdd). Use accessions, not short names: short
names change between releases.

---

## 6. BLAST+ Program Naming

**Symptom:** The checks stop the run with `program not found on PATH: rpsblast`. Typically the system has `rpsblast+` instead of `rpsblast`.

**Cause:** Different BLAST+ installations use different binary names. The bioconda package uses `rpsblast`, but some system installations use `rpsblast+`.

**Solution:** Use the optional `programs:` section in `config.yaml` to override program names:

```yaml
programs:
  rpsblast: 'rpsblast+'
  tblastn: 'tblastn'
  blast_formatter: 'blast_formatter'
  makeblastdb: 'makeblastdb'
  rpsbproc: 'rpsbproc'
  makeprofiledb: 'makeprofiledb'
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
./rpsHunter --all
```

---

## 8. Large Genome Downloads and NCBI Rate Limits

**Symptom:** Genome downloads are extremely slow or fail with HTTP 429 (Too Many Requests) errors.

**Cause:** NCBI rate-limits unauthenticated requests to 3 requests per second.

**Solution:**

1. Set an NCBI API key as an environment variable (it is never written in the config):
   ```bash
   export NCBI_API_KEY="your_api_key_here"
   ```
   This increases the rate limit from 3 to 10 requests per second.

2. For very large genomes, downloads may take hours even with an API key. Plan accordingly and ensure a stable network connection.

3. If a download is interrupted, run the same command again: genomes already fetched are skipped and the interrupted one starts over.

---
