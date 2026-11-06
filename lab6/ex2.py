
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path
import csv


FASTA_FILE = Path(__file__).parent / "GCF_002270725.1_ViralProj399942_genomic.fna"
OUTPUT_FILE = Path(__file__).parent / "possible_transposons.csv"


def read_fasta(filename):
    dna = ""
    with open(filename) as file:
        for line in file:
            if not line.startswith(">"):
                dna += line.strip().upper()
    return dna


def reverse_complement(s):
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(pairs[base] for base in reversed(s))


def find_transposons(dna, writer):
    counts = {}

    for k in range(6, 9):
        # The second sliding window visits all positions, indexed by sequence.
        # This avoids comparing the same window pairs millions of times.
        second_windows = defaultdict(list)
        for j in range(len(dna) - k + 1):
            window2 = dna[j:j + k]
            second_windows[window2].append(j)

        count = 0
        # Move the first sliding window through the entire genome.
        for i in range(len(dna) - k + 1):
            window1 = dna[i:i + k]
            target = reverse_complement(window1)

            # Find matching second windows strictly after the first one.
            locations = second_windows[target]
            first_valid = bisect_left(locations, i + k)
            for j in locations[first_valid:]:
                # 1-based coordinates, inclusive at both ends.
                writer.writerow([i + 1, j + k, k, window1, dna[j:j + k]])
                count += 1

        counts[k] = count

    return counts


def main():
    dna = read_fasta(FASTA_FILE)
    if not dna or any(base not in "ATGC" for base in dna):
        raise ValueError("The FASTA sequence must contain only A, T, G and C.")

    with open(OUTPUT_FILE, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["start", "end", "repeat_length", "left_ir", "right_ir"])
        counts = find_transposons(dna, writer)

    print("DNA length:", len(dna))
    for k, count in counts.items():
        print(f"Length {k}: {count} possible matches")
    print("Total:", sum(counts.values()))
    print("Results saved to:", OUTPUT_FILE)


if __name__ == "__main__":
    main()



