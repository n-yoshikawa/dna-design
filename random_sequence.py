import numpy

aminos = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

seqs = []
while True:
    seqs.append(''.join(numpy.random.choice(aminos, size=3)))
    seqs = list(set(seqs))
    if len(seqs) == 4:
        break
newSeqs = [s for s in seqs]
for seq in seqs:
    newSeqs.append(numpy.random.choice(aminos)+seq[1]+seq[2])
    newSeqs.append(seq[0]+numpy.random.choice(aminos)+seq[2])
    newSeqs.append(seq[0]+seq[1]+numpy.random.choice(aminos))

for s in sorted(newSeqs):
    print(s)

