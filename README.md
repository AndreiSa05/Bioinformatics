# Bioinformatics Labs

This repository contains Python and C++ exercises completed during the Bioinformatics laboratory sessions.

## Lab 1 – DNA Sequence Analysis

- **Exercise 1:** Count the frequency of each nucleotide (A, T, G, C) in a DNA sequence.
- **Exercise 2:** Count the frequencies of all overlapping combinations of 2 and 3 nucleotides.

## Lab 2 – Start/Stop Codons and Nested Sequences

- Convert DNA sequences to RNA.
- Find all `AUG` start codons.
- Check for in-frame stop codons (`UAA`, `UAG`, `UGA`).
- Identify complete sequences and count nested subsequences.
- Sort the sequences in descending order by their number of nested subsequences.

## Lab 3 – Viral Genome Analysis

- **Exercise 1:** Read the Abisko virus genome from a FASTA file, count the frequencies of all 64 possible nucleotide triplets using a sliding window of 3 positions, and display the 10 most frequent.
- **Exercise 2:** Use the OpenAI API to generate food and drink recommendations based on the virus identity, without deriving medical recommendations from codon frequencies.

## Lab 4 – Genome Reconstruction

- **Exercise 1:** Read the Abisko virus genome from a FASTA file, generate 20,000 random DNA fragments of 100–150 nucleotides, and save them to a file.
- **Exercise 2:** Reconstruct the original genome using a greedy algorithm that merges overlapping fragments, considering only overlaps of at least 10 nucleotides. Compare the reconstructed genome with the original to verify the result.

## Lab 5 – Gel Electrophoresis and DNA Repetition Detection

- **Exercise 1 (Python):** Read a DNA sequence of 2,000 nucleotides from a FASTA file, generate 10 random DNA fragments of different lengths, calculate their molecular weights, and simulate their migration through an electrophoresis gel using Matplotlib.
- **Exercise 2 (C++):** Read a DNA sequence from a FASTA file, generate all possible nucleotide combinations of lengths 3–5 using backtracking, and identify combinations that occur between 6 and 10 times in the sequence.

**Note:** Exercise 2 was implemented independently, without using AI (scary, I know).

## Technologies

- Python
- C++ (Lab 5, Exercise 2)
- Matplotlib (Lab 5, Exercise 1)
- OpenAI API (Lab 3, Exercise 2)