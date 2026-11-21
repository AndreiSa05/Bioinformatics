import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"

MATCH = 1
MISMATCH = -1
GAP = -2


def create_matrix(S1, S2):
    rows = len(S1) + 1
    cols = len(S2) + 1

    # Create a matrix filled with zeros
    matrix = [[0 for j in range(cols)] for i in range(rows)]

    # Initialize first column
    for i in range(rows):
        matrix[i][0] = i * GAP

    # Initialize first row
    for j in range(cols):
        matrix[0][j] = j * GAP

    # Fill the scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):

            if S1[i - 1] == S2[j - 1]:
                score = MATCH
            else:
                score = MISMATCH

            diagonal = matrix[i - 1][j - 1] + score
            top = matrix[i - 1][j] + GAP
            left = matrix[i][j - 1] + GAP

            matrix[i][j] = max(diagonal, top, left)

    return matrix


def traceback(S1, S2, matrix):
    aligned1 = ""
    aligned2 = ""

    i = len(S1)
    j = len(S2)

    path = [(i, j)]

    while i > 0 or j > 0:

        # Check diagonal movement
        if i > 0 and j > 0:
            if S1[i - 1] == S2[j - 1]:
                score = MATCH
            else:
                score = MISMATCH

            if matrix[i][j] == matrix[i - 1][j - 1] + score:
                aligned1 += S1[i - 1]
                aligned2 += S2[j - 1]

                i -= 1
                j -= 1

                path.append((i, j))
                continue

        # Check upward movement
        if i > 0 and matrix[i][j] == matrix[i - 1][j] + GAP:
            aligned1 += S1[i - 1]
            aligned2 += "-"

            i -= 1

        # Otherwise, move left
        else:
            aligned1 += "-"
            aligned2 += S2[j - 1]

            j -= 1

        path.append((i, j))

    return aligned1[::-1], aligned2[::-1], path


def display_results(aligned1, aligned2, matrix):
    matches = 0
    mismatches = 0
    gaps = 0
    alignment = ""

    for i in range(len(aligned1)):
        if aligned1[i] == aligned2[i]:
            alignment += "|"
            matches += 1

        elif aligned1[i] == "-" or aligned2[i] == "-":
            alignment += " "
            gaps += 1

        else:
            alignment += "X"
            mismatches += 1

    print("Aligned sequences:")
    print(aligned1)
    print(alignment)
    print(aligned2)

    print("\nMatches:", matches)
    print("Mismatches:", mismatches)
    print("Gaps:", gaps)
    print("Alignment score:", matrix[-1][-1])
    print("Identity:", matches / len(aligned1) * 100, "%")


def plot_matrix(S1, S2, matrix, path):
    fig, ax = plt.subplots(figsize=(11, 10))

    image = ax.imshow(matrix, cmap="coolwarm")

    # Display the nucleotides on the axes
    ax.set_xticks(range(len(S2) + 1))
    ax.set_yticks(range(len(S1) + 1))

    ax.set_xticklabels(["-"] + list(S2))
    ax.set_yticklabels(["-"] + list(S1))

    # Display the score inside each cell
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            ax.text(j, i, matrix[i][j],
                    ha="center", va="center", fontsize=8)

    # Highlight the traceback path
    for i, j in path:
        rectangle = Rectangle(
            (j - 0.5, i - 0.5),
            1, 1,
            fill=False,
            edgecolor="limegreen",
            linewidth=2.5
        )

        ax.add_patch(rectangle)

    ax.set_xlabel("Sequence 2")
    ax.set_ylabel("Sequence 1")
    ax.set_title("Needleman-Wunsch Scoring Matrix")

    plt.colorbar(image, label="Alignment score")
    plt.tight_layout()
    fig.savefig("scoring_matrix.png", dpi=300, bbox_inches="tight")


def plot_alignment(aligned1, aligned2):
    n = len(aligned1)

    fig, ax = plt.subplots(figsize=(14, 3))

    for i in range(n):

        # Select the color for each position
        if aligned1[i] == aligned2[i]:
            color = "lightgreen"
            symbol = "|"

        elif aligned1[i] == "-" or aligned2[i] == "-":
            color = "orange"
            symbol = " "

        else:
            color = "lightcoral"
            symbol = "X"

        # Display nucleotides from sequence 1
        ax.text(i, 2, aligned1[i],
                ha="center", va="center",
                fontsize=14,
                bbox=dict(facecolor=color, edgecolor="black"))

        # Display match or mismatch symbol
        ax.text(i, 1.5, symbol,
                ha="center", va="center", fontsize=14)

        # Display nucleotides from sequence 2
        ax.text(i, 1, aligned2[i],
                ha="center", va="center",
                fontsize=14,
                bbox=dict(facecolor=color, edgecolor="black"))

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0.5, 2.5)

    ax.set_xticks(range(n))
    ax.set_xticklabels(range(1, n + 1))

    ax.set_yticks([1, 2])
    ax.set_yticklabels(["S2", "S1"])

    ax.set_xlabel("Position")
    ax.set_title("Needleman-Wunsch DNA Alignment")

    plt.tight_layout()
    fig.savefig("dna_alignment.png", dpi=300, bbox_inches="tight")


# Main
matrix = create_matrix(S1, S2)

aligned1, aligned2, path = traceback(S1, S2, matrix)

display_results(aligned1, aligned2, matrix)

plot_matrix(S1, S2, matrix, path)

plot_alignment(aligned1, aligned2)

plt.show()