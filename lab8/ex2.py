import matplotlib.pyplot as plt
from pathlib import Path

# DNA sequence from the whiteboard
S = "CGGACTGATCTATCT" + "A" * 27 + "CGTAGCATCTATCGATCTATCTAGCGATCTA"

WINDOW_SIZE = 30

# Save plots next to the Python file
OUTPUT_DIR = Path(__file__).resolve().parent


# 1. Calculate GC content
def gc_content(S):
    gc = S.count("G") + S.count("C")
    return gc / len(S) * 100


# 2. Calculate information content
# Exact formula from the whiteboard
def information_content(S):
    n = len(S)
    total = 0

    for u in range(1, n):
        matches = 0

        for i in range(n - u):
            if S[i] == S[i + u]:
                matches += 1

        total += matches / ((n - u) * 100)

    return 100 - total / (n - 1)


# 3. Generate the DNA pattern using sliding windows
def generate_pattern(S):
    gc_values = []
    ic_values = []

    for i in range(len(S) - WINDOW_SIZE + 1):
        window = S[i:i + WINDOW_SIZE]

        gc_values.append(gc_content(window))
        ic_values.append(information_content(window))

    return gc_values, ic_values


# 4. Calculate the center of weight
def center_of_weight(gc_values, ic_values):
    center_x = sum(gc_values) / len(gc_values)
    center_y = sum(ic_values) / len(ic_values)

    return center_x, center_y


# 5. Plot the DNA pattern
def plot_pattern(gc_values, ic_values, center_x, center_y):
    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(gc_values, ic_values, "o-", color="blue",
            markersize=4, label="DNA pattern")

    ax.scatter(center_x, center_y, color="red",
               s=100, marker="X", label="Center of weight")

    ax.set_xlabel("GC Content (%)")
    ax.set_ylabel("Information Content")
    ax.set_title("DNA Pattern - Sliding Window 30")
    ax.grid(True)
    ax.legend()

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "dna_pattern.png", dpi=300,
                bbox_inches="tight")


# 6. Plot the centers of the patterns
def plot_centers(centers):
    fig, ax = plt.subplots(figsize=(9, 6))

    for name, x, y in centers:
        ax.scatter(x, y, s=100, label=name)
        ax.annotate(name, (x, y), xytext=(5, 5),
                    textcoords="offset points")

    ax.set_xlabel("GC Content (%)")
    ax.set_ylabel("Information Content")
    ax.set_title("Centers of Weight - DNA Patterns")
    ax.grid(True)
    ax.legend()

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "promoter_centers.png", dpi=300,
                bbox_inches="tight")


# 7. Process one DNA sequence
def process_sequence(name, sequence):
    gc_values, ic_values = generate_pattern(sequence)

    center_x, center_y = center_of_weight(gc_values, ic_values)

    print("\nSequence:", name)
    print("Length:", len(sequence))
    print("Number of windows:", len(gc_values))
    print("GC content:", round(gc_content(sequence), 2))
    print("Information content:", round(information_content(sequence), 2))
    print("Center of weight:", round(center_x, 2),
          round(center_y, 2))

    return gc_values, ic_values, center_x, center_y


# Main
gc_values, ic_values, center_x, center_y = process_sequence(
    "Artificial DNA", S
)

# Plot the DNA pattern
plot_pattern(gc_values, ic_values, center_x, center_y)

# Plot the center of weight
centers = [("Artificial DNA", center_x, center_y)]
plot_centers(centers)

# Display both plots
plt.show()
