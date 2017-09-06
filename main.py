from z3 import *

def primer_design(sequence_list, size):
    length = 12
    bases = ['A', 'C', 'G', 'T']

    x = [[{base: Int("x-{}-{}-{}".format(primer, pos, base)) for base in bases} for pos in range(length)] for primer in range(size)]

    opt = Optimize()

# binary restriction
    for prime in range(size):
        for pos in range(length):
            for base in bases:
                opt.add(x[prime][pos][base] >= 0)
                opt.add(x[prime][pos][base] <= 1)

# gap is not allowed
    for i in range(size):
        for j in range(length):
            opt.add(sum([x[i][j][b] for b in bases]) >= 1)


# F: TTT, TTC
    def aminoF(primer, pos):
        return Or(And(x[primer][3*pos+0]['T']==1,
                      x[primer][3*pos+1]['T']==1,
                      x[primer][3*pos+2]['T']==1),
                  And(x[primer][3*pos+0]['T']==1,
                      x[primer][3*pos+1]['T']==1,
                      x[primer][3*pos+2]['A']==1))
# L: TTA, TTG, CU*
    def aminoL(primer, pos):
        return Or(And(x[primer][3*pos+0]['T']==1,
                      x[primer][3*pos+1]['T']==1,
                      x[primer][3*pos+2]['A']==1),
                  And(x[primer][3*pos+0]['T']==1,
                      x[primer][3*pos+1]['T']==1,
                      x[primer][3*pos+2]['G']==1),
                  And(x[primer][3*pos+0]['C']==1,
                      x[primer][3*pos+1]['U']==1))

# L: ATT, ATC, ATA
    def aminoI(primer, pos):
        return Or(And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['T']==1,
                      x[primer][3*pos+2]['T']==1),
                  And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['T']==1,
                      x[primer][3*pos+2]['C']==1),
                  And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['T']==1,
                      x[primer][3*pos+2]['A']==1))

# M: ATG
    def aminoM(primer, pos):
        return Or(And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['T']==1,
                      x[primer][3*pos+2]['G']==1))
# V: GT*
    def aminoV(primer, pos):
        return And(x[primer][3*pos+0]['G']==1,
                   x[primer][3*pos+1]['T']==1)
# S: TC*
    def aminoS(primer, pos):
        return And(x[primer][3*pos+0]['T']==1,
                   x[primer][3*pos+1]['C']==1)
# S: CC*
    def aminoP(primer, pos):
        return And(x[primer][3*pos+0]['C']==1,
                   x[primer][3*pos+1]['C']==1)

# T: AC*
    def aminoT(primer, pos):
        return And(x[primer][3*pos+0]['A']==1,
                   x[primer][3*pos+1]['C']==1)
# A: GC*
    def aminoA(primer, pos):
        return And(x[primer][3*pos+0]['G']==1,
                   x[primer][3*pos+1]['C']==1)
# Y: TAT, TAC
    def aminoY(primer, pos):
        return Or(And(x[primer][3*pos+0]['T']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['T']==1),
                  And(x[primer][3*pos+0]['T']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['C']==1))
# H: CAT, CAC
    def aminoH(primer, pos):
        return Or(And(x[primer][3*pos+0]['C']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['T']==1),
                  And(x[primer][3*pos+0]['C']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['C']==1))
# Q: CAA, CAG
    def aminoQ(primer, pos):
        return Or(And(x[primer][3*pos+0]['C']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['A']==1),
                  And(x[primer][3*pos+0]['C']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['G']==1))

# N: AAT, AAC 
    def aminoN(primer, pos):
        return Or(And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['T']==1),
                  And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['C']==1))

# K: AAA, AAG
    def aminoK(primer, pos):
        return Or(And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['A']==1),
                  And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['G']==1))

# D: GAT, GAC
    def aminoD(primer, pos):
        return Or(And(x[primer][3*pos+0]['G']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['T']==1),
                  And(x[primer][3*pos+0]['G']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['C']==1))
# E: GAT, GAC
    def aminoE(primer, pos):
        return Or(And(x[primer][3*pos+0]['G']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['A']==1),
                  And(x[primer][3*pos+0]['G']==1,
                      x[primer][3*pos+1]['A']==1,
                      x[primer][3*pos+2]['G']==1))
# C: TGT, TGC
    def aminoC(primer, pos):
        return Or(And(x[primer][3*pos+0]['T']==1,
                      x[primer][3*pos+1]['G']==1,
                      x[primer][3*pos+2]['T']==1),
                  And(x[primer][3*pos+0]['T']==1,
                      x[primer][3*pos+1]['G']==1,
                      x[primer][3*pos+2]['C']==1))
# W: TGG
    def aminoW(primer, pos):
        return Or(And(x[primer][3*pos+0]['T']==1,
                      x[primer][3*pos+1]['G']==1,
                      x[primer][3*pos+2]['G']==1))
# R: CG*, AGA, AGG
    def aminoR(primer, pos):
        return Or(And(x[primer][3*pos+0]['C']==1,
                      x[primer][3*pos+1]['G']==1),
                  And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['G']==1,
                      x[primer][3*pos+2]['A']==1),
                  And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['G']==1,
                      x[primer][3*pos+2]['G']==1))
# S: AGT, AGC
    def aminoS(primer, pos):
        return Or(And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['G']==1,
                      x[primer][3*pos+2]['T']==1),
                  And(x[primer][3*pos+0]['A']==1,
                      x[primer][3*pos+1]['G']==1,
                      x[primer][3*pos+2]['C']==1))
# G: GG*
    def aminoG(primer, pos):
        return Or(And(x[primer][3*pos+0]['G']==1,
                      x[primer][3*pos+1]['G']==1))

    def dict2base(model, dic):
        if model[dic['A']]==1 and model[dic['C']]==0 and model[dic['G']]==0 and model[dic['T']]==0:
            return 'A'
        if model[dic['A']]==0 and model[dic['C']]==1 and model[dic['G']]==0 and model[dic['T']]==0:
            return 'C'
        if model[dic['A']]==0 and model[dic['C']]==0 and model[dic['G']]==1 and model[dic['T']]==0:
            return 'G'
        if model[dic['A']]==0 and model[dic['C']]==0 and model[dic['G']]==0 and model[dic['T']]==1:
            return 'T'
        if model[dic['A']]==1 and model[dic['C']]==0 and model[dic['G']]==1 and model[dic['T']]==0:
            return 'R'
        if model[dic['A']]==0 and model[dic['C']]==1 and model[dic['G']]==0 and model[dic['T']]==1:
            return 'Y'
        if model[dic['A']]==0 and model[dic['C']]==1 and model[dic['G']]==1 and model[dic['T']]==0:
            return 'S'
        if model[dic['A']]==1 and model[dic['C']]==0 and model[dic['G']]==0 and model[dic['T']]==1:
            return 'W'
        if model[dic['A']]==0 and model[dic['C']]==0 and model[dic['G']]==1 and model[dic['T']]==1:
            return 'K'
        if model[dic['A']]==1 and model[dic['C']]==1 and model[dic['G']]==0 and model[dic['T']]==0:
            return 'M'
        if model[dic['A']]==0 and model[dic['C']]==1 and model[dic['G']]==1 and model[dic['T']]==1:
            return 'B'
        if model[dic['A']]==1 and model[dic['C']]==0 and model[dic['G']]==1 and model[dic['T']]==1:
            return 'D'
        if model[dic['A']]==1 and model[dic['C']]==1 and model[dic['G']]==0 and model[dic['T']]==1:
            return 'H'
        if model[dic['A']]==1 and model[dic['C']]==1 and model[dic['G']]==1 and model[dic['T']]==0:
            return 'V'
        else:
            return 'N'

    def split_n(text, n):
        return [text[i*n:i*n+n] for i in range(len(text)/n)]

    amino2func = {}
    amino2func['A'] = aminoA
    amino2func['R'] = aminoR
    amino2func['N'] = aminoN
    amino2func['D'] = aminoD
    amino2func['C'] = aminoC
    amino2func['Q'] = aminoQ
    amino2func['E'] = aminoE
    amino2func['G'] = aminoG
    amino2func['H'] = aminoH
    amino2func['I'] = aminoI
    amino2func['L'] = aminoL
    amino2func['K'] = aminoK
    amino2func['M'] = aminoM
    amino2func['F'] = aminoF
    amino2func['P'] = aminoP
    amino2func['S'] = aminoS
    amino2func['T'] = aminoT
    amino2func['W'] = aminoW
    amino2func['Y'] = aminoY
    amino2func['V'] = aminoV

    def seq2cond(s, p):
        return And([amino2func[s[i]](p, i) for i in range(4)])



    base = {}
    base['A'] = ['A']
    base['C'] = ['C']
    base['G'] = ['G']
    base['T'] = ['T']
    base['R'] = ['A', 'G']
    base['Y'] = ['C', 'T']
    base['S'] = ['G', 'C']
    base['W'] = ['A', 'T']
    base['K'] = ['G', 'T']
    base['M'] = ['A', 'C']
    base['B'] = ['C', 'G', 'T']
    base['D'] = ['A', 'G', 'T']
    base['H'] = ['A', 'C', 'T']
    base['V'] = ['A', 'C', 'G']
    base['N'] = ['A', 'C', 'G', 'T']

    codon = {}
    codon["TTT"] = 'F'
    codon["TTC"] = 'F'
    codon["TTA"] = 'L'
    codon["TTG"] = 'L'
    codon["CTT"] = 'L'
    codon["CTC"] = 'L'
    codon["CTA"] = 'L'
    codon["CTG"] = 'L'
    codon["ATT"] = 'I'
    codon["ATC"] = 'I'
    codon["ATA"] = 'I'
    codon["ATG"] = 'M'
    codon["GTT"] = 'V'
    codon["GTC"] = 'V'
    codon["GTA"] = 'V'
    codon["GTG"] = 'V'

    codon["TCT"] = 'S'
    codon["TCC"] = 'S'
    codon["TCA"] = 'S'
    codon["TCG"] = 'S'
    codon["CCT"] = 'P'
    codon["CCC"] = 'P'
    codon["CCA"] = 'P'
    codon["CCG"] = 'P'
    codon["ACT"] = 'T'
    codon["ACC"] = 'T'
    codon["ACA"] = 'T'
    codon["ACG"] = 'T'
    codon["GCT"] = 'A'
    codon["GCC"] = 'A'
    codon["GCA"] = 'A'
    codon["GCG"] = 'A'

    codon["TAT"] = 'Y'
    codon["TAC"] = 'Y'
    codon["TAA"] = '-'
    codon["TAG"] = '-'
    codon["CAT"] = 'H'
    codon["CAC"] = 'H'
    codon["CAA"] = 'Q'
    codon["CAG"] = 'Q'
    codon["AAT"] = 'N'
    codon["AAC"] = 'N'
    codon["AAA"] = 'K'
    codon["AAG"] = 'K'
    codon["GAT"] = 'D'
    codon["GAC"] = 'D'
    codon["GAA"] = 'E'
    codon["GAG"] = 'E'

    codon["TGT"] = 'C'
    codon["TGC"] = 'C'
    codon["TGA"] = '-'
    codon["TGG"] = 'W'
    codon["CGT"] = 'R'
    codon["CGC"] = 'R'
    codon["CGA"] = 'R'
    codon["CGG"] = 'R'
    codon["AGT"] = 'S'
    codon["AGC"] = 'S'
    codon["AGA"] = 'R'
    codon["AGG"] = 'R'
    codon["GGT"] = 'G'
    codon["GGC"] = 'G'
    codon["GGA"] = 'G'
    codon["GGG"] = 'G'

    def generated_amino(seqs):
        cands = []
        for seq in seqs:
            c = []
            for b1 in base[seq[0]]:
                for b2 in base[seq[1]]:
                    for b3 in base[seq[2]]:
                        c.append(codon[b1+b2+b3])
            cands.append(c)

        ret = []
        for c1 in cands[0]:
            for c2 in cands[1]:
                for c3 in cands[2]:
                    for c4 in cands[3]:
                        ret.append(c1+c2+c3+c4)
        return ret

    for seq in sequence_list:
        opt.add(Or([seq2cond(seq, p) for p in range(size)]))

# objective function
# smaller ambiguity is better solution
    cost = Int('cost')
    opt.add(cost==sum([x[i][j][b] for b in bases for j in range(length) for i in range(size)]))

    h = opt.minimize(cost)
    if opt.check() == sat:
        model = opt.model()
        size_all = 0
        for prime in range(size):
                result = "".join([dict2base(model, x[prime][pos]) for pos in range(length)])
                print "primer:", result
                print("generated amino acids:")
                for aminos in generated_amino(split_n(result, 3)):
                    print aminos, '-', 'in list' if aminos in sequence_list else 'Not in list'
                    size_all += 1
        print 'Hit:{} / All:{}'.format(len(sequence_list), size_all)
    else:
        print 'not satisfiable'


sequence_list = ["GAYY", "GAHY", "GAFY", "GAYH", "GGYY", "GGHY", "GAHH", "GAFH", "GGFY", "GGYH", "GAYW", "GGHH", "GAHF", "GAFF"]
primer_design(sequence_list, 2)
