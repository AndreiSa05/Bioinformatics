import random


random.seed(42)

BASES = "ATGC"
DNA_LENGTH = 300
TRANSPOSON_COUNT = 3
IR_LENGTH = 6
INNER_LENGTH = 18
DR_LENGTH = 5


def random_dna(length):
    return "".join(random.choice(BASES) for _ in range(length))


def reverse_complement(sequence):
    pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(pairs[base] for base in reversed(sequence))


def create_transposon():
    left_ir = random_dna(IR_LENGTH)
    middle = random_dna(INNER_LENGTH)
    right_ir = reverse_complement(left_ir)
    return left_ir + middle + right_ir


def create_sequence():
    # Each insertion adds one transposon and one extra copy of a direct repeat.
    element_length = 2 * IR_LENGTH + INNER_LENGTH
    host_length = DNA_LENGTH - TRANSPOSON_COUNT * (element_length + DR_LENGTH)
    host = random_dna(host_length)
    insertion_points = [35, 95, 155]

    dna = ""
    expected_positions = []
    previous = 0

    for point in insertion_points:
        dna += host[previous:point]

        # Copy the host's target site, giving DR + transposon + DR.
        direct_repeat = host[point:point + DR_LENGTH]
        transposon = create_transposon()
        dna += direct_repeat

        start = len(dna) + 1  # Use 1-based positions for the lab output.
        dna += transposon
        end = len(dna)
        expected_positions.append((start, end))

        previous = point  # The original target site stays in the host sequence.

    dna += host[previous:]
    return dna, expected_positions


def detect_transposons(dna):
    positions = []
    element_length = 2 * IR_LENGTH + INNER_LENGTH

    for start in range(DR_LENGTH, len(dna) - element_length - DR_LENGTH + 1):
        end = start + element_length

        left_ir = dna[start:start + IR_LENGTH]
        right_ir = dna[end - IR_LENGTH:end]
        left_dr = dna[start - DR_LENGTH:start]
        right_dr = dna[end:end + DR_LENGTH]

        if reverse_complement(left_ir) == right_ir and left_dr == right_dr:
            positions.append((start + 1, end))

    return positions


dna, expected = create_sequence()
found = detect_transposons(dna)

print("DNA sequence:")
print(dna)
print("\nDNA length:", len(dna))

print("\nInserted transposons (start, end):")
for start, end in expected:
    print(start, end)

print("\nDetected transposons (start, end):")
for start, end in found:
    print(start, end)

print("\nDetection correct:", found == expected)
