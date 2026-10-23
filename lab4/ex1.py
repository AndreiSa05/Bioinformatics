"""Create 20,000 shuffled, overlapping reads from a nucleotide FASTA genome."""

import random
from pathlib import Path


FILENAME = "GCF_002270725.1_ViralProj399942_genomic.fna"
OUTPUT = "samples.txt"


def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip().upper()
    genome = sequence
    if not genome or any(base not in "ACGT" for base in genome):
        raise ValueError("Expected a nonempty nucleotide FASTA containing only A, C, G, T")
    return genome


def main():
    genome = read_fasta(FILENAME)
    if len(genome) < 150:
        raise ValueError("Genome must have at least 150 bases")

    rng = random.Random(42)
    samples = []
    for _ in range(20000):
        length = rng.randint(100, 150)
        start = rng.randint(0, len(genome) - length)
        samples.append(genome[start:start + length])

    rng.shuffle(samples)
    Path(OUTPUT).write_text("\n".join(samples) + "\n", encoding="utf-8")
    average_coverage = sum(map(len, samples)) / len(genome)
    print(f"Genome: {len(genome)} bases")
    print(f"Saved {len(samples)} reads (100–150 bases) to {OUTPUT}")
    print(f"Average sequencing coverage: {average_coverage:.1f}x")


if __name__ == "__main__":
    main()
