# What is searched, kept and reported

rpsHunter finds regions of a genome that resemble the query proteins and describes the protein
domains inside them. It does not survey a domain family genome-wide: a domain in a gene that does
not resemble any query is never seen. Counts are counts of query-like regions.

This page follows a sequence through the pipeline and says, at each step, what decides its fate.
The decisions behind each step are recorded in [`adr/`](adr/README.md).

## 1. tBLASTn: finding query-like regions

Each query protein is searched against each genome, translated in six frames. One hit covers one
exon, or part of one.

The search is deliberately permissive (`blast.e_value: 10`). On the three annotated model genomes a
stricter cut-off loses real domains long before it removes noise: in human the search reaches 95
annotated KRAB domains at E <= 10, 35 at E <= 1 and 8 at E <= 0.001
([ADR-010](adr/ADR-010-scope-and-hmmer-evidence.md)). Noise is handled later, by the domain
searches, which have their own statistics.

**Choosing queries.** A query attracts everything that resembles any part of it. A query with a
repeat array (the zinc-finger array of PRDM9, for example) will hit every member of that repeat
family, thousands of regions in a mammal. That is the right behaviour when the family is the
subject. When only orthologs of one gene are wanted, remove the repeat array from the query and
keep the domains that are specific to the gene, as the published PRDM9 surveys do. The shipped
example keeps full-length queries on purpose, to show both behaviours.

## 2. The quality gate: the only filter on hits

A hit is selected when it meets all four thresholds of the `blast:` section: percent identity,
E-value, alignment length and bit score. Nothing else selects hits. Every selected hit is searched
for domains.

Each hit also carries `Internal_Stop`: whether its translation contains a stop codon before its
last residue. Such a hit is kept, because a stop does not make a region uninteresting (it may be a
pseudogene, a frameshift, or an alignment that ran past an exon end), but the flag lets you tell
intact sequence from decayed sequence ([ADR-011](adr/ADR-011-retire-orf-gate.md)).

`fastas/{query}/{species}.parquet` keeps every hit with `Quality_Pass` and `Selected`, so the
effect of the gate can always be inspected.

## 3. Two independent domain searches

The selected hits' translations are searched twice, with different methods and libraries:

| | RPS-BLAST | HMMER |
|---|---|---|
| Library | NCBI CDD models of the families (`domains: <family>: cdd:`) | Pfam profiles of the families (`domains: <family>: pfam:`) |
| Role | the domain calls | independent evidence |
| Cut-off | `rpsblast.e_value` (default 1) | `hmmer:` section |

Neither search filters the other's input ([ADR-010](adr/ADR-010-scope-and-hmmer-evidence.md)).

### RPS-BLAST E-values are on NCBI's scale

The database searched holds only the configured families' models, which keeps the search fast. An
E-value depends on database size, so the search is told the size of the whole CDD release
(`databases.cdd.residues`). E-values are then the ones NCBI's CD-Search would report and do not
change when the list of families does ([ADR-006](adr/ADR-006-rpsblast-evalue-scale-and-cutoff.md)).

Every call gets an **evidence tier** from its E-value:

| Tier | E-value | Meaning |
|---|---|---|
| `standard` | <= 0.01 | CD-Search's default cut-off |
| `borderline` | <= 1 | NCBI's "may be real" range |
| `rescue` | > 1 | present only because `rpsblast.e_value` was raised on purpose |

`rpsbproc` then keeps representative calls per region and labels each as specific, non-specific or
superfamily, using CDD's own thresholds.

### HMMER

HMMER reports its own domains at their own genomic positions. Each says whether it clears Pfam's
curated gathering threshold (`Above_GA`), which matters because the shipped search settings are
permissive. Two light filters apply to what HMMER reports: `hmmer.min_coverage` (fraction of the
profile aligned) and `hmmer.min_alignment_length`.

### Method agreement

A CDD call and a HMMER domain agree when they report the same family, on the same hit, at
overlapping genomic positions. `Method_Agreement` is `both`, `cdd_only` / `hmmer_only`, or
`not_compared` (the family has models in only one library, or HMMER is off).

## 4. Where a domain is

A call is placed on the genome from the residues it covers inside the hit, three bases per residue,
on the hit's strand ([ADR-007](adr/ADR-007-domain-coordinates.md)). Stop codons are masked, not
removed, so that later residues keep their position.

Overlapping calls of one model on one strand are one **domain instance** (`tables/domains/domains.*`):
the same domain is often found through several queries and several overlapping hits.

`Incomplete` on a call is rpsbproc's flag that the alignment misses an end of the model. For a
domain split across exons this is the normal state of every piece, so it describes the piece, not
the gene.

## 5. Loci: putting exons back together

Selected hits within `loci.gap` bases of each other (default 20,000) on a sequence are chained into
a locus, on either strand ([ADR-008](adr/ADR-008-loci.md)). 20 kb keeps the exons of a gene
together without fusing tandem copies.

For each locus and family, `tables/loci/locus_families.*` gives:

- `Evidence`: the best tier among the family's calls in the locus;
- `Coverage`: the largest fraction of any one model covered when the locus's calls are laid on the
  model's own coordinates. A SET domain split over three exons, each piece "incomplete", adds up
  to a coverage near 1 ([ADR-009](adr/ADR-009-completeness-across-exons.md));
- `Copies`: the largest number of separate places one model occurs at (the size of a repeat array).

A locus lists the families it carries and gives no verdict of "complete": which families make a
whole gene depends on the protein family studied. A locus of chance hits is listed too; it carries
no family (`N_Families` 0).

## What cannot be concluded

- **Absence.** No hit means no region resembling the queries was found in this assembly, not that
  the gene is absent from the species.
- **Function.** A locus with every family at full coverage may still be a pseudogene; check
  `N_Hits_With_Stop`.
- **Family size.** Only query-like members are counted (see section 1).
