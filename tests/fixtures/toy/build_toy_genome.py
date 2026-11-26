"""Build the toy genomes used by the end-to-end test.

Three small synthetic "species", each one FASTA of random DNA into which copies
of a query protein are planted at known places. Because the planted positions
are known, the test can check what the pipeline reports against the truth
instead of only checking that it runs.

  Toyus_plenus   the whole protein, once on each strand
  Toyus_fractus  two fragments of it: the N-terminal part and the zinc fingers
  Toyus_vacuus   nothing planted (the zero-hit path)

A protein is planted by back-translating it with one fixed codon per amino acid,
so a tBLASTn search with the same protein finds it as a single exact match. Real
genes have introns; the toy does not, on purpose: it tests the plumbing, not the
biology.

Everything is derived from the seed, so the same seed gives byte-identical
files. Coordinates in the truth table are 1-based and inclusive, on the forward
strand of the scaffold, as in BLAST and GFF3.

Usage:
    python build_toy_genome.py --protein query.fa --dest DIR [--seed N]
Writes DIR/<species>.fa for each species and DIR/truth.tsv.
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from dataclasses import dataclass
from pathlib import Path

DEFAULT_SEED = 20261126
SCAFFOLD_LENGTH = 40_000  # bp; long enough to keep planted copies far apart
LINE_WIDTH = 60

# One codon per amino acid (a common human codon for each). Any fixed choice
# works: tBLASTn translates the DNA back, so the search sees the exact protein.
CODON = {
    "A": "GCC", "C": "TGC", "D": "GAC", "E": "GAG", "F": "TTC",
    "G": "GGC", "H": "CAC", "I": "ATC", "K": "AAG", "L": "CTG",
    "M": "ATG", "N": "AAC", "P": "CCC", "Q": "CAG", "R": "CGG",
    "S": "AGC", "T": "ACC", "V": "GTG", "W": "TGG", "Y": "TAC",
}  # fmt: skip
COMPLEMENT = str.maketrans("ACGT", "TGCA")

TRUTH_COLUMNS = (
    "species",
    "scaffold",
    "start",
    "end",
    "strand",
    "protein",
    "aa_start",
    "aa_end",
)


@dataclass(frozen=True)
class Plant:
    """One planted copy: residues aa_start..aa_end of the protein (1-based)."""

    scaffold: str
    offset: int  # 0-based position of the first planted base on the scaffold
    strand: str  # "+" or "-"
    aa_start: int
    aa_end: int


def layout(protein_length: int) -> dict[str, list[Plant]]:
    """Return what is planted where, for a protein of the given length.

    The fragment boundaries are fractions of the protein so that the layout
    works for any query: the first 30% covers the N-terminal domains of a
    protein like PRDM9, and the last 40% its zinc-finger array.
    """
    head_end = protein_length * 3 // 10
    tail_start = protein_length * 6 // 10 + 1
    return {
        "Toyus_plenus": [
            Plant("plenus_scaffold_1", 5_000, "+", 1, protein_length),
            Plant("plenus_scaffold_2", 12_000, "-", 1, protein_length),
        ],
        "Toyus_fractus": [
            Plant("fractus_scaffold_1", 8_000, "+", 1, head_end),
            Plant("fractus_scaffold_1", 25_000, "+", tail_start, protein_length),
        ],
        "Toyus_vacuus": [],
    }


def back_translate(protein: str) -> str:
    """Return DNA encoding the protein with one fixed codon per residue."""
    return "".join(CODON[residue] for residue in protein)


def reverse_complement(dna: str) -> str:
    """Return the reverse complement of an upper-case DNA string."""
    return dna.translate(COMPLEMENT)[::-1]


def random_dna(length: int, rng: random.Random) -> str:
    """Return random DNA of the given length drawn from the generator."""
    return "".join(rng.choices("ACGT", k=length))


def read_single_fasta(path: Path) -> tuple[str, str]:
    """Return (identifier, sequence) of a FASTA file holding one record."""
    lines = path.read_text().splitlines()
    headers = [line for line in lines if line.startswith(">")]
    if len(headers) != 1:
        raise ValueError(f"{path}: expected one FASTA record, found {len(headers)}")
    identifier = headers[0][1:].split()[0]
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    unknown = set(sequence) - set(CODON)
    if unknown:
        raise ValueError(f"{path}: residues without a codon: {sorted(unknown)}")
    return identifier, sequence


def build_species(
    plants: list[Plant], scaffolds: list[str], protein: str, rng: random.Random
) -> dict[str, str]:
    """Return {scaffold name: sequence} with every plant written into place."""
    sequences = {name: random_dna(SCAFFOLD_LENGTH, rng) for name in scaffolds}
    for plant in plants:
        coding = back_translate(protein[plant.aa_start - 1 : plant.aa_end])
        insert = coding if plant.strand == "+" else reverse_complement(coding)
        background = sequences[plant.scaffold]
        end = plant.offset + len(insert)
        if end > len(background):
            raise ValueError(f"{plant.scaffold}: plant runs past the scaffold end")
        sequences[plant.scaffold] = (
            background[: plant.offset] + insert + background[end:]
        )
    return sequences


def write_fasta(sequences: dict[str, str], destination: Path) -> None:
    """Write the scaffolds as a FASTA file wrapped at LINE_WIDTH."""
    with destination.open("w") as handle:
        for name, sequence in sequences.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), LINE_WIDTH):
                handle.write(sequence[start : start + LINE_WIDTH] + "\n")


def truth_rows(species: str, plants: list[Plant], protein_id: str) -> list[dict]:
    """Return one truth-table row per plant, in 1-based inclusive coordinates."""
    rows = []
    for plant in plants:
        length = 3 * (plant.aa_end - plant.aa_start + 1)
        rows.append(
            {
                "species": species,
                "scaffold": plant.scaffold,
                "start": plant.offset + 1,
                "end": plant.offset + length,
                "strand": plant.strand,
                "protein": protein_id,
                "aa_start": plant.aa_start,
                "aa_end": plant.aa_end,
            }
        )
    return rows


def scaffold_names(species: str, plants: list[Plant]) -> list[str]:
    """Return the species' scaffolds: those with plants, or one empty scaffold."""
    names = list(dict.fromkeys(plant.scaffold for plant in plants))
    return names or [f"{species.split('_')[1]}_scaffold_1"]


def build(protein_fasta: Path, destination: Path, seed: int = DEFAULT_SEED) -> Path:
    """Write the toy genomes and their truth table; return the truth table path."""
    protein_id, protein = read_single_fasta(protein_fasta)
    rng = random.Random(seed)  # reproducible test data, not security
    destination.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for species, plants in layout(len(protein)).items():
        scaffolds = scaffold_names(species, plants)
        sequences = build_species(plants, scaffolds, protein, rng)
        write_fasta(sequences, destination / f"{species}.fa")
        rows.extend(truth_rows(species, plants, protein_id))
    truth_path = destination / "truth.tsv"
    with truth_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRUTH_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return truth_path


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point; returns the process exit status."""
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--dest", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args(argv)
    truth_path = build(args.protein, args.dest, args.seed)
    print(f"toy genomes and {truth_path.name} written to {args.dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
