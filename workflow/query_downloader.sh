#!/bin/bash

# Ensure correct usage
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    echo "Usage: $0 <AccessionCode_or_SearchTerm> <OutputDirectory> <LogFile>"
    exit 1
fi

QUERY="$1"       # Protein Accession, Gene name, or Search Term
OUTDIR="$(realpath -m "$2")"  # Normalize and clean output directory
LOGFILE="$(realpath -m "$3")" # Log file to store results

# Check if an Entrez email is provided
if [ -z "$ENTREZ_EMAIL" ]; then
    echo "⚠ Warning: No Entrez email found. Consider setting it with: export ENTREZ_EMAIL='your_email@example.com'"
    EMAIL_FLAG=""
else
    EMAIL_FLAG="-email $ENTREZ_EMAIL"
    echo "Found Entrez email in environment: $ENTREZ_EMAIL"
fi

echo "Fetching protein for: $QUERY"

################################################################################
# Check if input is a Protein Accession (e.g. ADK09900)
################################################################################
if [[ "$QUERY" =~ ^[A-Z]{3}[0-9]{5}(\.[0-9]+)?$ ]]; then
    echo "Detected Protein Accession: $QUERY"
    BEST_ACCESSION="$QUERY"
    LEVEL="Direct_Accession"
    QUERY_SAFE="$QUERY"
else
    ################################################################################
    # Use Entrez Direct to search for the protein accession using the search term
    ################################################################################
    BEST_ACCESSION=$(esearch -db protein -query "$QUERY" $EMAIL_FLAG | efetch -format acc $EMAIL_FLAG | head -n 1)
    if [ -z "$BEST_ACCESSION" ]; then
        echo "❌ No valid protein accession found for: $QUERY"
        exit 1
    fi
    echo "✅ Found protein accession: $BEST_ACCESSION"
    LEVEL="Searched"
    QUERY_SAFE="$QUERY"
fi

# Convert spaces in the level to underscores for a filename-safe string
LEVEL_SAFE="${LEVEL// /_}"

# Normalize filename using realpath -m
FASTA_FILE="$(realpath -m "$OUTDIR/${QUERY_SAFE}.fa")"

echo "Downloading protein accession: $BEST_ACCESSION"
efetch -db protein -id "$BEST_ACCESSION" -format fasta $EMAIL_FLAG > "$FASTA_FILE"

if [ $? -ne 0 ] || [ ! -s "$FASTA_FILE" ]; then
    echo "❌ Download failed for protein accession: $BEST_ACCESSION"
    exit 1
fi

echo "✅ Download complete: $FASTA_FILE (Accession: $BEST_ACCESSION, Source: $LEVEL)"

# Append download info to the log file
echo -e "$QUERY\t$BEST_ACCESSION\t$LEVEL" >> "$LOGFILE"
echo "✅ Download logged: $LOGFILE"
