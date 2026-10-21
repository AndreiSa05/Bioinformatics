S = "ATGGCCAATGAGTGGAATCG"

freq = {}

for ch in S:
    freq[ch] = freq.get(ch, 0) + 1

print(freq)
print(len(S))


