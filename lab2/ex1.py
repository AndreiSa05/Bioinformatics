
S = "ATGATGATGTAAATGTAA"

# Convert DNA to RNA
S = S.upper().replace("T", "U")

stops = ["UAA", "UAG", "UGA"]
subsequences = []
count = 0

# Find all AUG start codons
for i in range(len(S) - 2):
    if S[i:i+3] == "AUG":
        count += 1

        # Find the first stop codon in this reading frame
        for j in range(i + 3, len(S) - 2, 3):
            if S[j:j+3] in stops:
                subsequences.append((i, j + 3))
                break

print("Total AUG start codons:", count)
print("Complete sequences:", len(subsequences))

# Count nested subsequences
results = []

for start, end in subsequences:
    nested = 0

    for other_start, other_end in subsequences:
        if other_start > start and other_end <= end:
            nested += 1

    results.append((nested, start, end))

# Sort descending by number of nested subsequences
results.sort(reverse=True)

print("\nSequences sorted by nested count:")

for nested, start, end in results:
    print(S[start:end], "->", nested, "nested subsequences")

if results:
    nested, start, end = results[0]
    print("\nSequence with most nested subsequences:")
    print(S[start:end])
    print("Nested count:", nested)
