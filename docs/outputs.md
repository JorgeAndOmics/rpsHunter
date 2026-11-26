# Output files

Everything a run writes, with the columns of each table. Results go under `root.results_root_folder`,
logs under `root.logs_root_folder`. `{query}` is a query's label from the config and `{species}` a
species key. Tables are written as Parquet; the run-level tables also as CSV beside it.

Coordinates are 1-based and inclusive with `Start <= End`; strand is its own column. A hit's `Tag`
is a 12-character identifier derived from the hit itself, the same in every run, and is what joins
tables to each other.

## Layout

```
results/
  blast/{query}/{species}.parquet        every tBLASTn hit
  asn/tblastn/{query}/{species}.asn      BLAST archive of the search
  fastas/{query}/{species}.fa            translations of the selected hits
  fastas/{query}/{species}.parquet       audit: every hit with the gate's verdict
  hmmer/{query}/{species}.parquet        HMMER domains (when HMMER runs)
  rpsblast/{query}/{species}.parquet     raw RPS-BLAST alignments
  asn/rpsblast/{query}/{species}.asn     BLAST archive, input of rpsbproc
  rpsbproc/{query}/{species}.txt         rpsbproc's report
  domains/{query}/{species}.parquet      domain calls of one query and species
  ranges/{species}.gff3                  domain instances as a genome-browser track
  tables/
    blast/{query}/selected/{species}.*   selected hits
    blast/{query}/aggregate.*            every hit of a query, all species
    blast/aggregate.*                    every hit of the run
    rpsblast/...                         raw alignments, aggregated the same way
    domains/calls.*                      every domain call of the run
    domains/domains.*                    domain instances
    domains/hmmer_domains.*              HMMER domains (when HMMER runs)
    loci/loci.*                          loci
    loci/locus_families.*                what each locus carries, per family
```

The tables to start from are `tables/loci/loci`, `tables/loci/locus_families` and
`tables/domains/domains`.

## Hits

### `blast/{query}/{species}.parquet`

One row per tBLASTn hit.

| Column | Meaning |
|---|---|
| `Query ID`, `Query_Accession` | The query protein |
| `Subject ID` | Genome sequence the hit is on |
| `Pct Identity`, `Alignment Length`, `Mismatches`, `Gap Openings` | BLAST's alignment statistics |
| `Q. Start`, `Q. End` | Residues of the query aligned |
| `S. Start`, `S. End` | Genomic span as BLAST writes it: start greater than end means the minus strand |
| `E-value`, `Bit Score` | BLAST's scores |
| `Subject Sequence` | The hit's translation, with `-` for alignment gaps and `*` for stop codons |
| `Species`, `Species_Name` | Species key and display name |
| `Tag` | The hit's identifier |
| `Internal_Stop` | A stop codon lies inside the translation (before its last residue) |

### `fastas/{query}/{species}.parquet` (audit)

The same rows with `Quality_Pass` (meets the four thresholds of `blast:`) and `Selected` (goes on to
the domain searches; today identical to `Quality_Pass`).

### `tables/blast/{query}/selected/{species}.*` and `fastas/{query}/{species}.fa`

The selected hits, and their translations as searched: alignment gaps dropped, stop codons masked
with `X`. Records are named `<sequence>:<S. Start>-<S. End>|tag:<Tag>`.

## Domains

### `tables/domains/calls.*`

One row per domain call: one model reported by rpsbproc on one hit.

| Column | Meaning |
|---|---|
| `Species`, `Species_Name`, `Query_Accession`, `Tag` | Where the call comes from |
| `Sequence`, `Strand`, `Start`, `End` | The domain's own position on the genome |
| `Hit_Start`, `Hit_End` | The hit it lies in |
| `Family` | Family from the config's `domains:` block |
| `Model`, `Model_Name`, `PSSM_ID` | CDD accession, short name and numeric id of the model or superfamily reported |
| `Hit_Type` | rpsbproc's class: Specific, Non-specific or Superfamily |
| `Evalue`, `Bitscore` | On the scale of the full CDD release |
| `Evidence` | `standard` (E <= 0.01), `borderline` (E <= 1) or `rescue` (above) |
| `Residue_From`, `Residue_To` | Residues of the hit's translation covered |
| `Aligned_PSSM`, `Model_From`, `Model_To`, `Model_Length` | The model the alignment was made to, and the part of it covered |
| `Incomplete` | rpsbproc's flag: `N`, `C` or `NC` when the alignment misses that end of the model, `-` otherwise. Normal for one exon of a split domain |
| `Method_Agreement` | `both`, `cdd_only` or `not_compared` (see [filtering.md](filtering.md#method-agreement)) |

`domains/{query}/{species}.parquet` holds the same rows for one query and species, without
`Method_Agreement`.

### `tables/domains/domains.*`

One row per **domain instance**: overlapping calls of one model on one strand, counted once however
many hits and queries found them. `Start` and `End` span the merged calls; `Model_Name`, `Hit_Type`,
`Evidence`, `Evalue` and `Bitscore` are those of the best call. `N_Calls` is the number of calls
merged, `Queries` and `Tags` list where they came from.

### `tables/domains/hmmer_domains.*`

One row per HMMER domain. `Profile` and `Pfam` name the profile; `Evalue` and `Score` are the
domain's own (independent E-value and bit score); `Above_GA` says whether the score clears Pfam's
gathering threshold; `Method_Agreement` is `both`, `hmmer_only` or `not_compared`. Position and
model-coverage columns are as in the call table.

### `ranges/{species}.gff3`

The domain instances as GFF3 features of type `protein_match`, with the best E-value as score and
the family, model and evidence tier as attributes. Written for every species; a species with no
domain gets a file with only the header, so "no domains" is never confused with "not run".

## Loci

### `tables/loci/loci.*`

One row per locus: selected hits chained within `loci.gap` bases.

| Column | Meaning |
|---|---|
| `Species`, `Species_Name` | |
| `Locus` | `<sequence>:<start>-<end>`, unique within a species |
| `Sequence`, `Start`, `End`, `Length` | Extent of the locus |
| `N_Hits`, `N_Hits_With_Stop` | Hits chained, and how many have an internal stop codon |
| `Queries` | Queries that found it |
| `Best_Identity`, `Best_Bitscore` | Of its best hit |
| `Families`, `N_Families` | Families with at least one call, in config order; empty and 0 for a locus with no domain |

### `tables/loci/locus_families.*`

One row per locus and family present.

| Column | Meaning |
|---|---|
| `Species`, `Locus`, `Family` | |
| `Evidence`, `Best_Evalue` | Best tier and E-value among the family's calls in the locus |
| `Coverage` | Largest fraction of one model covered by the locus's calls together (exons added up) |
| `Best_Model` | The model that coverage refers to |
| `Copies` | Largest number of separate places one model occurs at |
| `N_Calls` | Calls of the family in the locus |

## Empty results

A species or query with nothing to report gets valid empty files: tables with all their columns and
no rows, an empty FASTA, a header-only GFF3. A tool that fails stops its job instead, so an empty
table always means "searched, nothing found".

## Logs and benchmarks

| Path | Content |
|------|---------|
| `jobs/<rule>/<query>/<species>.log` | Everything one job printed. Per-query rules use `jobs/<rule>/<query>.log`; rules that run once use `jobs/<rule>/all.log` |
| `benchmarks/<rule>/...tsv` | Snakemake's benchmark for the same job: wall-clock seconds, peak memory, I/O |

When a job fails, its log is printed to the console and its outputs are removed.
