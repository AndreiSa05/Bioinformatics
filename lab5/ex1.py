import random
import math
import matplotlib.pyplot as plt


# 1. Read DNA sequence from FASTA file
def read_fasta(path):
    with open(path, "r") as file:
        dna = "".join(line.strip() for line in file if not line.startswith(">"))
    return dna.upper()


dna = read_fasta("GCF_002270725.1_ViralProj399942_genomic.fna")
dna = dna[:2000]

print("DNA sequence length:", len(dna))


# 2. Generate 10 random DNA samples
samples = []

for i in range(10):
    length = random.randint(100, len(dna))
    start = random.randint(0, len(dna) - length)

    sample = dna[start:start + length]
    samples.append(sample)

print("\nDNA Samples:")

for i in range(10):
    print(f"Sample {i + 1}: {len(samples[i])} nucleotides")


# 3. Calculate molecular weight
def molecular_weight(length):
    return length * 617.96 + 36.04


print("\nMolecular Weights:")

for i in range(10):
    length = len(samples[i])
    weight = molecular_weight(length)

    print(f"Sample {i + 1}: {length} bp, MW = {weight / 1000:.2f} kDa")


# 4. Simulate DNA migration
gel_length = 10
run_time = 40
max_size = 20000

k = gel_length / (run_time * (math.log10(max_size) - math.log10(100)))


def migration_distance(length, time):
    distance = k * (math.log10(max_size) - math.log10(length)) * time
    return min(distance, gel_length)


lengths = [len(sample) for sample in samples]

print("\nMigration Distances:")

for i in range(10):
    distance = migration_distance(lengths[i], run_time)

    print(f"Sample {i + 1}: {lengths[i]} bp -> {distance:.2f} cm")


# 5. Visualize gel electrophoresis
ladder = [100, 200, 300, 400, 500, 700, 1000, 1500, 2000]

fig, ax = plt.subplots(figsize=(12, 8))

ax.set_facecolor("#17233b")

# DNA ladder
for length in ladder:
    y = migration_distance(length, run_time)

    ax.hlines(y, -0.35, 0.35,
              colors="#ffff80", linewidth=5)

    ax.text(-0.5, y, f"{length} bp",
            color="white", fontsize=8,
            ha="right", va="center")


# DNA samples
for i in range(10):
    x = i + 1
    y = migration_distance(lengths[i], run_time)

    ax.hlines(y, x - 0.35, x + 0.35,
              colors="#67f5d0", linewidth=6)

    ax.text(x, y + 0.3, f"{lengths[i]} bp",
            color="white", fontsize=8, ha="center")


# Draw wells
for x in range(11):
    ax.add_patch(
        plt.Rectangle(
            (x - 0.35, -0.35),
            0.7, 0.25,
            color="#555a66"
        )
    )


ax.set_xlim(-1.5, 10.7)
ax.set_ylim(10.5, -0.6)

ax.set_xticks(range(11))
ax.set_xticklabels(["L"] + [str(i) for i in range(1, 11)])

ax.set_xlabel("DNA Samples (L = DNA Ladder)")
ax.set_ylabel("Migration Distance (cm)")
ax.set_title("Gel Electrophoresis Simulation - 10 DNA Samples")

plt.tight_layout()
plt.show()