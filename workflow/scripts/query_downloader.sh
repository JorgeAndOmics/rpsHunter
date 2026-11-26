#!/usr/bin/env bash
# Download one query protein from NCBI by its versioned accession.
#
# The accession is a key of `queries:` in the config. It is fetched exactly:
# no search by name and no "best match", so the same config always gets the
# same sequence.
#
# Usage: query_downloader.sh <accession> <output_dir> <record_file>
# Writes: <output_dir>/<accession>.fa
# Appends to <record_file>, tab-separated:
#   accession, sha256 of the FASTA, UTC time
# Needs: efetch (entrez-direct), sha256sum. ENTREZ_EMAIL, if set, is sent to
#   NCBI as its usage policy asks; NCBI_API_KEY is used by efetch if set.
# Fails on: wrong arguments, an accession without a version, a failed download,
#   or a reply that is not the FASTA record of that accession.
set -euo pipefail

if (($# != 3)); then
    echo "usage: $0 <accession> <output_dir> <record_file>" >&2
    exit 64
fi
accession="$1"
output_dir="$2"
record_file="$3"

if [[ ! "$accession" =~ ^[A-Z]{1,3}_?[0-9]+\.[0-9]+$ ]]; then
    echo "'$accession' is not a versioned protein accession (for example NP_000000.1)" >&2
    exit 65
fi

fasta="$output_dir/$accession.fa"
partial="$(mktemp --tmpdir="$output_dir" ".download_${accession}.XXXXXX")"
trap 'rm -f -- "$partial"' EXIT

email=()
if [[ -n "${ENTREZ_EMAIL:-}" ]]; then
    email=(-email "$ENTREZ_EMAIL")
fi

efetch -db protein -id "$accession" -format fasta "${email[@]}" > "$partial"

# efetch can exit 0 with an empty reply or an error page: check what came back.
header="$(head -n 1 "$partial")"
if [[ "$header" != ">$accession"* ]]; then
    echo "$accession: NCBI did not return its FASTA record (first line: '${header:0:80}')" >&2
    exit 66
fi

mv -- "$partial" "$fasta"
checksum="$(sha256sum "$fasta" | cut -d' ' -f1)"
printf '%s\t%s\t%s\n' "$accession" "$checksum" "$(date -u +%FT%TZ)" >> "$record_file"
echo "$accession: downloaded, sha256 $checksum"
