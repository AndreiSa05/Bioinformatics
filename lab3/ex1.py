from itertools import product

def read_fasta(filename):
    sequence = ""

    with open(filename, "r") as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip().upper()

    return sequence.replace("T", "U")

def get_virus_name(filename):
    with open(filename, "r") as file:
        header = file.readline().strip()

    return header

filename = "GCF_002270725.1_ViralProj399942_genomic.fna"
virus_name = get_virus_name(filename)
S = read_fasta(filename)

freq = {
    "".join(codon): 0
    for codon in product("AUGC", repeat=3)
}

for i in range(len(S) - 2):
    codon = S[i:i+3]

    if codon in freq:
        freq[codon] += 1

top10 = sorted(freq.items(), key=lambda x: x[1], reverse=True)[:10]

print("Top 10 most frequent codons:")

for codon, count in top10:
    print(codon, count)