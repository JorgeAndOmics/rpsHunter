#!/usr/bin/env bash
# Rebuild the miniature reference databases under tests/fixtures/toy/db/ from
# full, already-downloaded copies of CDD and Pfam.
#
# The end-to-end test must run without the 40 GB of real databases, so it uses
# a handful of profiles: enough to recognise every domain of the query protein.
# The files written here are committed; run this script only to refresh them
# after changing MODELS or moving to a newer CDD or Pfam release, and record the
# releases in db/PROVENANCE.
#
# Only models from sources that allow redistribution are included: NCBI-curated
# (cd) and COG models are public domain and Pfam is CC0. SMART-derived models
# (smart*) are deliberately left out because SMART restricts redistribution.
#
# Usage: make_mini_db.sh <db_root> <cdd_version> <pfam_release>
#   <db_root> holds cdd_smp/, rpsbproc_dbs/ and hmm_dbs/Pfam-A.hmm, as laid out
#   by the pipeline's database setup.
# Needs: hmmfetch (HMMER) on PATH.
set -euo pipefail

DB_ROOT="${1:?database root}"
CDD_VERSION="${2:?CDD version, e.g. 3.22}"
PFAM_RELEASE="${3:?Pfam release, e.g. 38.2}"
OUT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/db"

# CDD accessions of the profiles kept, with the family each one stands for.
MODELS=(
    pfam01352   # KRAB
    pfam09514   # SSXRD
    cd19193     # PR-SET_PRDM7_9, the SET domain of PRDM7 and PRDM9
    cd08161     # SET, the generic SET model
    pfam00096   # zf-C2H2
    pfam13465   # zf-H2C2_2
    COG5048     # ZMS1, a multi-finger zinc-finger model
)
# Pfam profile names for the HMMER side.
PFAM_NAMES=(KRAB SET SSXRD zf-C2H2)

rm -rf "$OUT"
mkdir -p "$OUT"/{cdd_smp,rpsbproc_dbs,hmm_dbs}

for model in "${MODELS[@]}"; do
    cp "$DB_ROOT/cdd_smp/$model.smp" "$OUT/cdd_smp/"
done
touch "$OUT/cdd_smp/.download_done"

# Keep the rows of each annotation table that concern the kept models and the
# superfamilies they belong to. Column 1 of family_superfamily_links is the model
# accession and column 3 its superfamily; the other tables key on column 2.
pattern="$(IFS='|'; echo "${MODELS[*]}")"
src="$DB_ROOT/rpsbproc_dbs"
awk -F'\t' -v re="^($pattern)$" '$1 ~ re' "$src/family_superfamily_links" \
    > "$OUT/rpsbproc_dbs/family_superfamily_links"
superfamilies="$(cut -f3 "$OUT/rpsbproc_dbs/family_superfamily_links" | sort -u | paste -sd'|')"
keep="^($pattern|$superfamilies)$"
for table in cddid.tbl cddannot.dat cddannot_generic.dat bitscore_specific.txt; do
    awk -F'\t' -v re="$keep" '$2 ~ re' "$src/$table" > "$OUT/rpsbproc_dbs/$table"
done
: > "$OUT/rpsbproc_dbs/cdtrack.txt"  # empty upstream as well (CDD 3.22)

printf '%s\n' "${PFAM_NAMES[@]}" > "$OUT/hmm_dbs/names.txt"
hmmfetch -f "$DB_ROOT/hmm_dbs/Pfam-A.hmm" "$OUT/hmm_dbs/names.txt" > "$OUT/hmm_dbs/Pfam-A.hmm"
rm "$OUT/hmm_dbs/names.txt"

cat > "$OUT/PROVENANCE" <<EOF
Miniature databases for the end-to-end test. Built by tests/fixtures/toy/make_mini_db.sh.
CDD version: $CDD_VERSION (https://ftp.ncbi.nih.gov/pub/mmdb/cdd/)
Pfam release: $PFAM_RELEASE (https://ftp.ebi.ac.uk/pub/databases/Pfam/)
CDD models: ${MODELS[*]}
Pfam profiles: ${PFAM_NAMES[*]}
EOF

echo "mini database written to $OUT:"
du -sh "$OUT"/*
wc -l "$OUT"/rpsbproc_dbs/*
grep -c '^NAME' "$OUT/hmm_dbs/Pfam-A.hmm"
