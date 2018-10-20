import itertools
import numpy as np
from collections import Counter
import sys

def Amino2Codon(amino):
    if amino == 'A':
        return ('GCA', 'GCC', 'GCG', 'GCU')
    if amino == 'R':
        return ('CGA', 'CGC', 'CGG', 'CGU', 'AGA', 'AGG')
    if amino == 'N':
        return ('AAU', 'AAC')
    if amino == 'D':
        return ('GAU', 'GAC')
    if amino == 'C':
        return ('UGU', 'UGC')
    if amino == 'Q':
        return ('CAA', 'CAG')
    if amino == 'E':
        return ('GAA', 'GAG')
    if amino == 'G':
        return ('GGA', 'GGC', 'GGG', 'GGU')
    if amino == 'H':
        return ('CAC', 'CAU')
    if amino == 'I':
        return ('AUA', 'AUC', 'AUU')
    if amino == 'L':
        return ('UUA', 'UUG', 'CUA', 'CUC', 'CUG', 'CUU')
    if amino == 'K':
        return ('AAA', 'AAG')
    if amino == 'M':
        return ('AUG')
    if amino == 'F':
        return ('UUU', 'UUC')
    if amino == 'P':
        return ('CCA', 'CCC', 'CCG', 'CCU')
    if amino == 'S':
        return ('AGC', 'AGU', 'UCA', 'UCC', 'UCG', 'UCU')
    if amino == 'T':
        return ('ACA', 'ACC', 'ACG', 'ACU')
    if amino == 'W':
        return ('UGG')
    if amino == 'Y':
        return ('UAU', 'UAC')
    if amino == 'V':
        return ('GUA', 'GUC', 'GUG', 'GUU')
    else:
        raise ValueError("Unknown amino acid")

def HammingDistance(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("The length is not the same.")
    dist = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            dist += 1
    return dist

def MinimumHammingDistance(centerBase, baseCandidates):
    minDist = sys.maxsize
    bestBase = None
    d = [HammingDistance(centerBase, base) for base in baseCandidates]
    idx = np.argmin(d)
    return d[idx], baseCandidates[idx]

def MinimumHammingDistance2(centerBase, baseCandidates):
    minDist = sys.maxsize
    bestBase = None
    return min([HammingDistance(centerBase, base) for base in baseCandidates])


def getBaseName(baseSet):
    if baseSet == set(['A']):
        return 'A'
    if baseSet == set(['C']):
        return 'C'
    if baseSet == set(['G']):
        return 'G'
    if baseSet == set(['U']):
        return 'U'
    if baseSet == set(['A', 'G']):
        return 'R'
    if baseSet == set(['C', 'U']):
        return 'Y'
    if baseSet == set(['G', 'C']):
        return 'S'
    if baseSet == set(['A', 'U']):
        return 'W'
    if baseSet == set(['G', 'U']):
        return 'K'
    if baseSet == set(['A', 'C']):
        return 'M'
    if baseSet == set(['C', 'G', 'U']):
        return 'B'
    if baseSet == set(['A', 'G', 'U']):
        return 'D'
    if baseSet == set(['A', 'C', 'U']):
        return 'H'
    if baseSet == set(['A', 'C', 'G']):
        return 'V'
    if baseSet == set(['A', 'C', 'G', 'U']):
        return 'N'
    else:
        raise ValueError("Unknown amino acid")

def split_n(text, n):
    return [text[i*n:i*n+n] for i in range(len(text)//n)]

baselist = {}
baselist['A'] = ['A']
baselist['C'] = ['C']
baselist['G'] = ['G']
baselist['U'] = ['U']
baselist['R'] = ['A', 'G']
baselist['Y'] = ['C', 'U']
baselist['S'] = ['G', 'C']
baselist['W'] = ['A', 'U']
baselist['K'] = ['G', 'U']
baselist['M'] = ['A', 'C']
baselist['B'] = ['C', 'G', 'U']
baselist['D'] = ['A', 'G', 'U']
baselist['H'] = ['A', 'C', 'U']
baselist['V'] = ['A', 'C', 'G']
baselist['N'] = ['A', 'C', 'G', 'U']


codon = {}
codon["UUU"] = 'F'
codon["UUC"] = 'F'
codon["UUA"] = 'L'
codon["UUG"] = 'L'
codon["CUU"] = 'L'
codon["CUC"] = 'L'
codon["CUA"] = 'L'
codon["CUG"] = 'L'
codon["AUU"] = 'I'
codon["AUC"] = 'I'
codon["AUA"] = 'I'
codon["AUG"] = 'M'
codon["GUU"] = 'V'
codon["GUC"] = 'V'
codon["GUA"] = 'V'
codon["GUG"] = 'V'
codon["UCU"] = 'S'
codon["UCC"] = 'S'
codon["UCA"] = 'S'
codon["UCG"] = 'S'
codon["CCU"] = 'P'
codon["CCC"] = 'P'
codon["CCA"] = 'P'
codon["CCG"] = 'P'
codon["ACU"] = 'U'
codon["ACC"] = 'U'
codon["ACA"] = 'U'
codon["ACG"] = 'U'
codon["GCU"] = 'A'
codon["GCC"] = 'A'
codon["GCA"] = 'A'
codon["GCG"] = 'A'
codon["UAU"] = 'Y'
codon["UAC"] = 'Y'
codon["UAA"] = '-'
codon["UAG"] = '-'
codon["CAU"] = 'H'
codon["CAC"] = 'H'
codon["CAA"] = 'Q'
codon["CAG"] = 'Q'
codon["AAU"] = 'N'
codon["AAC"] = 'N'
codon["AAA"] = 'K'
codon["AAG"] = 'K'
codon["GAU"] = 'D'
codon["GAC"] = 'D'
codon["GAA"] = 'E'
codon["GAG"] = 'E'
codon["UGU"] = 'C'
codon["UGC"] = 'C'
codon["UGA"] = '-'
codon["UGG"] = 'W'
codon["CGU"] = 'R'
codon["CGC"] = 'R'
codon["CGA"] = 'R'
codon["CGG"] = 'R'
codon["AGU"] = 'S'
codon["AGC"] = 'S'
codon["AGA"] = 'R'
codon["AGG"] = 'R'
codon["GGU"] = 'G'
codon["GGC"] = 'G'
codon["GGA"] = 'G'
codon["GGG"] = 'G'

class AminoAcid:
    def __init__(self, amino):
        self.amino = amino
        candidates = []
        for c in amino:
            candidates.append(Amino2Codon(c))
        self.bases = [''.join(x) for x in itertools.product(*candidates)]
        self.cluster = None
    def __repr__(self):
        return "<{}, {}, {}>".format(self.amino, self.base, self.cluster)

def clustering(aminos, k):
    # initialization
    baseLen = len(aminos[0].bases[0])
    bases = [np.random.choice(a.bases) for a in aminos]
    centroids = np.random.choice(bases, k, replace=False)

    previousClusters = [[] for _ in range(k)]
    for step in range(1000):
        # assignment step
        clusters = [[] for _ in range(k)]
        for amino in aminos:
            minDist = sys.maxsize
            minCluster = -1
            minBase = None
            for clusterIndex, centroid in enumerate(centroids):
                d, b = MinimumHammingDistance(centroid, amino.bases)
                if d < minDist:
                    minBase = b
                    minCluster = clusterIndex
                    minDist = d
            clusters[minCluster].append(minBase)

        # convergence check
        if clusters == previousClusters:
            break
        previousClusters = clusters

        # update step
        centroids = ["" for _ in range(k)]
        for i, cluster in enumerate(clusters):
            for j in range(baseLen):
                centroids[i] += Counter([base[j] for base in cluster]).most_common()[0][0]

    # generate covering primer
    primers = ["" for _ in range(k)]
    for i, cluster in enumerate(clusters):
        for j in range(len(cluster[0])):
            primers[i] += getBaseName(set([b[j] for b in cluster]))
    return primers

def generated_amino(seqs):
    cands = []
    for seq in seqs:
        cand = []
        for b1 in baselist[seq[0]]:
            for b2 in baselist[seq[1]]:
                for b3 in baselist[seq[2]]:
                    cand.append(codon[b1+b2+b3])
        cands.append(cand)

    ret = [''.join(x) for x in itertools.product(*cands)]

    return ret

def design(sequenceList, k):
    bestResultBases = None
    bestGeneratedAminos = None
    minScore = sys.maxsize
    for trial in range(10):
        resultBases = clustering([AminoAcid(s) for s in sequenceList], k)

        generatedAminos = []
        for r in resultBases:
            generatedAminos.append(generated_amino(split_n(r, 3)))

        score = sum([len(a) for a in generatedAminos])
        print(resultBases)
        print(generatedAminos)
        print(score)
        if score < minScore:
            minScore = score
            bestResultBases = resultBases
            bestGeneratedAminos = generatedAminos
        if minScore == len(sequenceList):
            break
    print("best result:")
    print(bestResultBases)
    print(bestGeneratedAminos)
    return (list(set(bestResultBases)), [list(set(a)) for a in bestGeneratedAminos])
