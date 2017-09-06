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


def dict2base(model, dic):
    A = True if model[dic['A']]==1 else False
    C = True if model[dic['C']]==1 else False
    G = True if model[dic['G']]==1 else False
    U = True if model[dic['U']]==1 else False

    if A and (not C) and (not G) and (not U):
        return 'A'
    if (not A) and C and (not G) and (not U):
        return 'C'
    if (not A) and (not C) and G and (not U):
        return 'G'
    if (not A) and (not C) and (not G) and U:
        return 'U'
    if A and (not C) and G and (not U):
        return 'R'
    if (not A) and C and (not G) and U:
        return 'Y'
    if (not A) and C and G and (not U):
        return 'S'
    if A and (not C) and (not G) and U:
        return 'W'
    if (not A) and (not C) and G and U:
        return 'K'
    if A and C and (not G) and (not U):
        return 'M'
    if (not A) and C and G and U:
        return 'B'
    if A and (not C) and G and U:
        return 'D'
    if A and C and (not G) and U:
        return 'H'
    if A and C and G and (not U):
        return 'V'
    if A and C and G and U:
        return 'N'
