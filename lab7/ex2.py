import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# We use 1-based positions in the printed output.
# A right inverted repeat is the REVERSE COMPLEMENT of the left repeat.


def reverse_complement(sequence):
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    result = ""
    for base in sequence[::-1]:
        result += pairs[base]
    return result


def random_dna(length):
    return "".join(random.choice("ATGC") for _ in range(length))


def insert_transposon(dna, start, end, left_repeat):
    """Set the repeats at the ends of [start, end) (0-based coordinates)."""
    right_repeat = reverse_complement(left_repeat)
    k = len(left_repeat)

    if end - start - 2 * k < 10:
        raise ValueError("The region between repeats must be at least 10 bases.")

    dna = dna[:start] + left_repeat + dna[start + k:]
    dna = dna[:end - k] + right_repeat + dna[end:]
    return dna


def find_transposons(dna, repeats):
    """Search the DNA; no insertion coordinates are used here."""
    found = []

    for name, left_repeat in repeats.items():
        right_repeat = reverse_complement(left_repeat)
        k = len(left_repeat)

        for start in range(len(dna) - k + 1):
            if dna[start:start + k] != left_repeat:
                continue

            # A minimum of 10 letters must lie BETWEEN the repeats.
            for right_start in range(start + k + 10, len(dna) - k + 1):
                if dna[right_start:right_start + k] == right_repeat:
                    end = right_start + k
                    found.append((name, start, end, dna[start:end]))

    return sorted(found, key=lambda item: item[1])


def show_results(title, dna, repeats, expected):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)
    print("DNA length:", len(dna))
    print("DNA sequence:", dna)

    detected = find_transposons(dna, repeats)
    print("\nDetected transposons:", len(detected))

    for name, start, end, sequence in detected:
        print("\n" + name)
        print("Start:", start + 1)
        print("End:", end)
        print("Sequence:", sequence)
        print("Length:", len(sequence))

    # Compare against the known insertions, only AFTER detection.
    detected_coordinates = {(name, start, end) for name, start, end, _ in detected}
    for name, start, end in expected:
        if (name, start, end) not in detected_coordinates:
            print("WARNING: Could not detect inserted", name)

    return detected


def assign_levels(detected):
    """Put overlapping transposons on different rows in the plot."""
    levels_last_end = []
    placed = []

    for item in detected:
        name, start, end, sequence = item
        level = None

        for i in range(len(levels_last_end)):
            if start >= levels_last_end[i]:
                level = i
                levels_last_end[i] = end
                break

        if level is None:
            levels_last_end.append(end)
            level = len(levels_last_end) - 1

        placed.append((name, start, end, sequence, level))

    return placed, len(levels_last_end)


def plot_case(title, dna, detected, filename):
    colors = ["#64B5E8", "#68C6A1", "#F4A261", "#C77DFF", "#E76F51"]

    placed, number_of_levels = assign_levels(detected)

    fig, ax = plt.subplots(figsize=(12, 3 + 0.7 * number_of_levels))
    dna_length = len(dna)

    ax.hlines(0, 1, dna_length, color="gray", linewidth=2)

    for i in range(0, dna_length + 1, 25):
        x = i if i > 0 else 1
        ax.vlines(x, -0.08, 0.08, color="lightgray", linewidth=0.8)

    for index, (name, start, end, sequence, level) in enumerate(placed):
        y = number_of_levels - 1 - level
        color = colors[index % len(colors)]

        rect = Rectangle((start + 1, y - 0.18), end - start, 0.36,
                         facecolor=color, edgecolor="black", alpha=0.9)
        ax.add_patch(rect)

        ax.text((start + end + 1) / 2, y, name,
                ha="center", va="center", fontsize=10, fontweight="bold")

        ax.text(start + 1, y + 0.28, str(start + 1),
                ha="center", va="bottom", fontsize=9)
        ax.text(end, y + 0.28, str(end),
                ha="center", va="bottom", fontsize=9)

    ax.set_xlim(1, dna_length)
    ax.set_ylim(-0.6, number_of_levels - 1 + 0.8)
    ax.set_yticks([])
    ax.set_xlabel("DNA position (bp)")
    ax.set_title(title)

    for spine in ["left", "right", "top"]:
        ax.spines[spine].set_visible(False)

    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    plt.show()


def plot_all_cases(cases):
    colors = ["#64B5E8", "#68C6A1", "#F4A261", "#C77DFF", "#E76F51"]
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))

    for ax, (title, dna, detected) in zip(axes, cases):
        placed, number_of_levels = assign_levels(detected)
        dna_length = len(dna)

        ax.hlines(0, 1, dna_length, color="gray", linewidth=2)

        for i in range(0, dna_length + 1, 25):
            x = i if i > 0 else 1
            ax.vlines(x, -0.08, 0.08, color="lightgray", linewidth=0.8)

        for index, (name, start, end, sequence, level) in enumerate(placed):
            y = number_of_levels - 1 - level
            color = colors[index % len(colors)]

            rect = Rectangle((start + 1, y - 0.18), end - start, 0.36,
                             facecolor=color, edgecolor="black", alpha=0.9)
            ax.add_patch(rect)

            ax.text((start + end + 1) / 2, y, name,
                    ha="center", va="center", fontsize=10, fontweight="bold")

            ax.text(start + 1, y + 0.28, str(start + 1),
                    ha="center", va="bottom", fontsize=9)
            ax.text(end, y + 0.28, str(end),
                    ha="center", va="bottom", fontsize=9)

        ax.set_xlim(1, dna_length)
        ax.set_ylim(-0.6, number_of_levels - 1 + 0.8)
        ax.set_yticks([])
        ax.set_xlabel("DNA position (bp)")
        ax.set_title(title)

        for spine in ["left", "right", "top"]:
            ax.spines[spine].set_visible(False)

    plt.tight_layout()
    plt.savefig("all_cases.png", dpi=200)
    plt.show()


random.seed(42)

# CASE 1: normal, separate transposons
normal_repeats = {"T1": "ATGCCGTA", "T2": "CGATTCGA"}
normal = random_dna(300)
normal_positions = [("T1", 30, 70), ("T2", 110, 150)]
for name, start, end in normal_positions:
    normal = insert_transposon(normal, start, end, normal_repeats[name])

# CASE 2: overlapping
overlap_repeats = {"T1": "ACGTAGCT", "T2": "TTGACGCA"}
overlap = random_dna(300)
overlap_positions = [("T1", 40, 110), ("T2", 90, 160)]
for name, start, end in overlap_positions:
    overlap = insert_transposon(overlap, start, end, overlap_repeats[name])

# CASE 3: nested
nested_repeats = {"Outer": "GCTAACGT", "Inner": "AACGTTGC"}
nested = random_dna(300)
nested_positions = [("Outer", 30, 230), ("Inner", 90, 150)]
for name, start, end in nested_positions:
    nested = insert_transposon(nested, start, end, nested_repeats[name])

normal_detected = show_results("CASE 1: NORMAL", normal, normal_repeats, normal_positions)
overlap_detected = show_results("CASE 2: OVERLAPPING", overlap, overlap_repeats, overlap_positions)
nested_detected = show_results("CASE 3: NESTED", nested, nested_repeats, nested_positions)

plot_case("CASE 1: NORMAL", normal, normal_detected, "case1_normal.png")
plot_case("CASE 2: OVERLAPPING", overlap, overlap_detected, "case2_overlapping.png")
plot_case("CASE 3: NESTED", nested, nested_detected, "case3_nested.png")

plot_all_cases([
    ("CASE 1: NORMAL", normal, normal_detected),
    ("CASE 2: OVERLAPPING", overlap, overlap_detected),
    ("CASE 3: NESTED", nested, nested_detected)
])