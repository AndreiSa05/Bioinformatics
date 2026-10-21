S = "ATGGCCAATGAGTGGAATCGCCCCCCCCCCCCCCCCCCCC"

freq2 = {}
freq3 = {}

for i in range(len(S) - 1):
    comb = S[i:i+2]
    freq2[comb] = freq2.get(comb, 0) + 1

for i in range(len(S) - 2):
    comb = S[i:i+3]
    freq3[comb] = freq3.get(comb, 0) + 1

for comb, count in freq2.items():
    print(comb, count)

print()

for comb, count in freq3.items():
    print(comb, count)

print(len(S))