"""Greedy overlap assembly of the reads produced by split_genome.py.

Only suffix-to-prefix overlaps of at least 10 bases are accepted.
Reads' original genome coordinates are never used for assembly.
"""

from collections import defaultdict
import heapq
from pathlib import Path

from ex1 import FILENAME, read_fasta


def remove_contained(reads, seed_length=12):
    """Discard duplicate reads and reads wholly inside a longer read."""
    index = defaultdict(set)
    kept = []
    for read in sorted(set(reads), key=len, reverse=True):
        seed = read[:seed_length]
        if any(read in kept[i] for i in index[seed]):
            continue
        number = len(kept)
        kept.append(read)
        for offset in range(len(read) - seed_length + 1):
            index[read[offset:offset + seed_length]].add(number)
    return kept


def assemble(reads, minimum=10):
    active = {}
    prefixes = defaultdict(set)
    suffixes = defaultdict(set)
    candidates = []
    next_id = 0

    def add(read):
        nonlocal next_id
        number = next_id
        next_id += 1
        active[number] = read
        prefixes[read[:minimum]].add(number)
        suffixes[read[-minimum:]].add(number)
        return number

    def remove(number):
        read = active.pop(number)
        prefixes[read[:minimum]].discard(number)
        suffixes[read[-minimum:]].discard(number)

    def find_edges(number):
        """Index 10-base seeds, then verify the entire suffix/prefix match."""
        read = active[number]
        right_matches = {}
        left_matches = {}
        # A read/contig may be longer than the original sample length.
        for size in range(minimum, len(read) + 1):
            tail = read[-size:]
            for other in prefixes.get(tail[:minimum], ()):
                if other != number and len(active[other]) >= size and active[other].startswith(tail):
                    right_matches[other] = size
            head = read[:size]
            for other in suffixes.get(head[-minimum:], ()):
                if other != number and len(active[other]) >= size and active[other].endswith(head):
                    left_matches[other] = size
        for other, size in right_matches.items():
            heapq.heappush(candidates, (-size, number, other))
        for other, size in left_matches.items():
            heapq.heappush(candidates, (-size, other, number))

    for read in reads:
        if len(read) >= minimum:
            add(read)
    for number in list(active):
        find_edges(number)

    merges = 0
    while candidates:
        negative_size, left, right = heapq.heappop(candidates)
        if left not in active or right not in active:
            continue  # Pair was already merged; this is an outdated edge.
        size = -negative_size
        if not active[left].endswith(active[right][:size]):
            continue

        merged = active[left] + active[right][size:]
        remove(left)
        remove(right)

        # A newly merged contig might completely contain other reads/contigs.
        contained = [number for number, read in active.items() if read in merged]
        for number in contained:
            remove(number)

        new_number = add(merged)
        find_edges(new_number)
        merges += 1

    return sorted(active.values(), key=len, reverse=True), merges


def main():
    samples_file = Path("samples.txt")
    output_file = Path("reconstructed_genome.fasta")
    min_overlap = 10

    reads = [line.strip().upper() for line in samples_file.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not reads or any(any(base not in "ACGT" for base in read) for read in reads):
        raise ValueError("Samples must contain nucleotide sequences, one per line")

    filtered = remove_contained(reads, seed_length=min(12, min_overlap))
    print(f"Read {len(reads)} fragments; {len(filtered)} remain after removing duplicates/contained reads", flush=True)
    contigs, merges = assemble(filtered, min_overlap)

    with output_file.open("w", encoding="utf-8") as output:
        for index, contig in enumerate(contigs, start=1):
            output.write(f">contig_{index} length={len(contig)}\n")
            for offset in range(0, len(contig), 80):
                output.write(contig[offset:offset + 80] + "\n")

    original = read_fasta(FILENAME)
    print(f"Performed {merges} greedy merges")
    print(f"Contigs: {len(contigs)}")
    print(f"Longest contig: {max(map(len, contigs))} bases")
    print(f"Sum of contig lengths: {sum(map(len, contigs))} bases")
    print(f"Original genome length: {len(original)} bases")
    print(f"Perfect reconstruction: {len(contigs) == 1 and contigs[0] == original}")
    print(f"Saved assembly to {output_file}")


if __name__ == "__main__":
    main()
