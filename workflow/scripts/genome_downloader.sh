#!/usr/bin/env bash
# Download one genome assembly from NCBI, exactly the version asked for.
#
# The assembly is named by a versioned accession (GCA_/GCF_ with .N) taken from
# the species entry of the config. Nothing is resolved here: no taxon lookup and
# no "latest", so the same config always fetches the same sequences.
#
# Usage: genome_downloader.sh <species_key> <assembly_accession> <output_dir> <record_file>
# Writes: <output_dir>/<species_key>.fa
# Appends to <record_file>, tab-separated:
#   species key, accession, sha256 of the FASTA, UTC time, datasets version
# Needs: datasets (ncbi-datasets-cli), unzip, sha256sum. NCBI_API_KEY is used by
#   datasets if it is set in the environment.
# Fails on: wrong arguments, an accession without a version, a download or
#   archive error, or an archive that holds a different assembly.
set -euo pipefail

if (($# != 4)); then
    echo "usage: $0 <species_key> <assembly_accession> <output_dir> <record_file>" >&2
    exit 64
fi
species="$1"
accession="$2"
output_dir="$3"
record_file="$4"

if [[ -z "$accession" ]]; then
    echo "$species: the config gives no assembly accession, so it cannot be downloaded;" \
        "place its FASTA at $3/$species.fa" >&2
    exit 65
fi
if [[ ! "$accession" =~ ^GC[AF]_[0-9]{9}\.[0-9]+$ ]]; then
    echo "$species: '$accession' is not a versioned assembly accession (GCA_/GCF_ with .N)" >&2
    exit 65
fi

fasta="$output_dir/$species.fa"
work="$(mktemp --directory --tmpdir="$output_dir" ".download_${species}.XXXXXX")"
trap 'rm -rf -- "$work"' EXIT

echo "$species: downloading assembly $accession"
datasets download genome accession "$accession" \
    --include genome --no-progressbar --filename "$work/genome.zip"

# The archive names its sequence files by assembly; a different accession inside
# means NCBI answered with another version than the one requested.
# Listed into a variable first, so that an unreadable archive stops the script
# here (a command inside a process substitution does not).
listing="$(unzip -Z -1 "$work/genome.zip" 'ncbi_dataset/data/*/*.fna')"
mapfile -t members <<< "$listing"
if [[ -z "$listing" ]]; then
    echo "$species: no sequence file in the archive for $accession" >&2
    exit 1
fi
for member in "${members[@]}"; do
    if [[ "$member" != "ncbi_dataset/data/$accession/"* ]]; then
        echo "$species: archive holds '$member', not assembly $accession" >&2
        exit 1
    fi
done

# Write beside the destination and rename, so an interrupted download never
# leaves a truncated genome under the final name.
unzip -p "$work/genome.zip" "${members[@]}" > "$work/genome.fa"
checksum="$(sha256sum "$work/genome.fa" | cut -d' ' -f1)"
mv -- "$work/genome.fa" "$fasta"

datasets_version="$(datasets --version | awk '{print $NF}')"
printf '%s\t%s\t%s\t%s\t%s\n' \
    "$species" "$accession" "$checksum" "$(date --utc '+%Y-%m-%dT%H:%M:%SZ')" \
    "$datasets_version" >> "$record_file"
echo "$species: wrote $fasta (sha256 $checksum)"
