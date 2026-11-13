from itertools import combinations
import matplotlib.pyplot as plt


def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                continue
            sequence += line
    return sequence.upper()


def find_cuts(sequence, recognition, offset):
    cuts = []
    for i in range(len(sequence) - len(recognition) + 1):
        if sequence[i:i + len(recognition)] == recognition:
            cuts.append(i + offset)
    return cuts


def get_fragments(cuts, dna_length):
    positions = [0] + sorted(set(cuts)) + [dna_length]
    fragments = []
    for i in range(len(positions) - 1):
        fragments.append(positions[i + 1] - positions[i])
    return fragments


def plot_gel(results):
    fig, ax = plt.subplots(figsize=(14, 8))
    all_lengths = []
    for lane, (names, cuts, fragments) in enumerate(results, start=1):
        all_lengths.extend(fragments)
        for length in fragments:
            ax.plot([lane - 0.32, lane + 0.32], [length, length],
                    color="deepskyblue", linewidth=3)

    ax.set_yscale("log")
    ax.set_ylim(max(1, min(all_lengths) * 0.75), max(all_lengths) * 1.3)
    ax.set_xlim(0.5, len(results) + 0.5)
    ax.set_xticks(range(1, len(results) + 1))
    ax.set_xlabel("Enzyme combination (lane number)")
    ax.set_ylabel("DNA fragment length (bp, logarithmic scale)")
    ax.set_title("Simulated Gel Electrophoresis - 3 Enzymes per Lane")
    ax.grid(axis="y", alpha=0.25)

    fig.subplots_adjust(bottom=0.30)
    for index, (names, cuts, fragments) in enumerate(results):
        column = index // 5
        row = index % 5
        fig.text(0.10 + column * 0.48, 0.21 - row * 0.037,
                 f"{index + 1}: {' + '.join(names)}", fontsize=9)

    plt.savefig("gel_electrophoresis.png", dpi=200)
    plt.show()


def plot_fragment_map(names, cuts, fragments, dna_length, cuts_by_enzyme):
    # DNA fragments alternate blue and green; enzyme cuts have separate colors.
    fig, ax = plt.subplots(figsize=(15, 4))
    start = 0
    for index, length in enumerate(fragments):
        ax.barh(0, length, left=start, height=0.45,
                color="skyblue" if index % 2 == 0 else "lightgreen",
                edgecolor="white")
        # Only label fragments wide enough to fit the text.
        if length >= dna_length * 0.06:
            ax.text(start + length / 2, 0, str(length),
                    ha="center", va="center", fontsize=9)
        start += length

    colors = {
        "EcoRI": "red",
        "BamHI": "purple",
        "HindIII": "green",
        "TaqI": "orange",
        "HaeIII": "blue"
    }

    for enzyme in names:
        for cut in cuts_by_enzyme[enzyme]:
            ax.plot([cut, cut], [-0.3, 0.3],
                    color=colors[enzyme], linewidth=1.5)
        # Empty line creates a legend entry for this enzyme.
        ax.plot([], [], color=colors[enzyme], label=enzyme, linewidth=2)

    ax.legend(title="Restriction enzymes", loc="upper right")
    ax.set_xlim(0, dna_length)
    ax.set_ylim(-0.5, 0.6)
    ax.set_yticks([])
    ax.set_xlabel("Position along DNA (bp)")
    ax.set_title("DNA Fragment Map: " + " + ".join(names))
    ax.text(0.5, -0.25,
            f"{len(cuts)} cuts produce {len(fragments)} fragments (linear DNA)",
            transform=ax.transAxes, ha="center", fontsize=10)
    plt.tight_layout()
    plt.savefig("fragment_map.png", dpi=200)
    plt.show()


filename = "GCF_002270725.1_ViralProj399942_genomic.fna"
sequence = read_fasta(filename)

enzymes = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "TaqI": ("TCGA", 1),
    "HaeIII": ("GGCC", 2)
}

cuts_by_enzyme = {}
for name, (recognition, offset) in enzymes.items():
    cuts_by_enzyme[name] = find_cuts(sequence, recognition, offset)
    print(f"{name}: {len(cuts_by_enzyme[name])} cuts")

print(f"\nDNA length: {len(sequence)} bp")

results = []
for number, names in enumerate(combinations(enzymes, 3), start=1):
    all_cuts = []
    for name in names:
        all_cuts.extend(cuts_by_enzyme[name])

    all_cuts = sorted(set(all_cuts))
    fragments = get_fragments(all_cuts, len(sequence))
    results.append((names, all_cuts, fragments))

    print(f"\nCombination {number}: {' + '.join(names)}")
    print("Combined cleavage positions:", all_cuts)
    print("Resulting fragment lengths:", fragments, "bp")
    print(f"{len(all_cuts)} cuts produce {len(fragments)} fragments (linear DNA)")

plot_gel(results)

# Choose which of the 10 combinations to show as a detailed DNA map.
selected_combination = 3    # Change this for other combinations.
names, cuts, fragments = results[selected_combination - 1]
plot_fragment_map(names, cuts, fragments, len(sequence), cuts_by_enzyme)
