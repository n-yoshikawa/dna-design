from z3 import *

# A: GC*
def aminoA(x, seq, pos):
    return And(x[seq][3*pos+0]['G']==1,
               x[seq][3*pos+1]['C']==1)
# R: CG*, AGA, AGG
def aminoR(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['C']==1,
                  x[seq][3*pos+1]['G']==1),
              And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['G']==1,
                  x[seq][3*pos+2]['A']==1),
              And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['G']==1,
                  x[seq][3*pos+2]['G']==1))
# N: AAU, AAC 
def aminoN(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['U']==1),
              And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['C']==1))
# D: GAU, GAC
def aminoD(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['G']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['U']==1),
              And(x[seq][3*pos+0]['G']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['C']==1))
# C: UGU, UGC
def aminoC(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['U']==1,
                  x[seq][3*pos+1]['G']==1,
                  x[seq][3*pos+2]['U']==1),
              And(x[seq][3*pos+0]['U']==1,
                  x[seq][3*pos+1]['G']==1,
                  x[seq][3*pos+2]['C']==1))
# Q: CAA, CAG
def aminoQ(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['C']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['A']==1),
              And(x[seq][3*pos+0]['C']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['G']==1))
# E: GAT, GAG
def aminoE(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['G']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['A']==1),
              And(x[seq][3*pos+0]['G']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['G']==1))
# G: GG*
def aminoG(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['G']==1,
                  x[seq][3*pos+1]['G']==1))
# H: CAU, CAC
def aminoH(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['C']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['U']==1),
              And(x[seq][3*pos+0]['C']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['C']==1))
# L: AUU, AUC, AUA
def aminoI(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['U']==1,
                  x[seq][3*pos+2]['U']==1),
              And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['U']==1,
                  x[seq][3*pos+2]['C']==1),
              And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['U']==1,
                  x[seq][3*pos+1]['A']==1))
# L: UUA, UUG, CU*
def aminoL(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['U']==1,
                  x[seq][3*pos+1]['U']==1,
                  x[seq][3*pos+2]['A']==1),
              And(x[seq][3*pos+0]['U']==1,
                  x[seq][3*pos+1]['U']==1,
                  x[seq][3*pos+2]['G']==1),
              And(x[seq][3*pos+0]['C']==1,
                  x[seq][3*pos+1]['U']==1))
# K: AAA, AAG
def aminoK(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['A']==1),
              And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['G']==1))
# M: AUG
def aminoM(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['U']==1,
                  x[seq][3*pos+2]['G']==1))
# F: UUU, UUC
def aminoF(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['U']==1,
                  x[seq][3*pos+1]['U']==1,
                  x[seq][3*pos+2]['U']==1),
              And(x[seq][3*pos+0]['U']==1,
                  x[seq][3*pos+1]['U']==1,
                  x[seq][3*pos+2]['C']==1))
# P: CC*
def aminoP(x, seq, pos):
    return And(x[seq][3*pos+0]['C']==1,
               x[seq][3*pos+1]['C']==1)
# S: UC*, AGU, AGC
def aminoS(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['U']==1,
                  x[seq][3*pos+1]['C']==1),
              And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['G']==1,
                  x[seq][3*pos+2]['U']==1),
              And(x[seq][3*pos+0]['A']==1,
                  x[seq][3*pos+1]['G']==1,
                  x[seq][3*pos+2]['C']==1))
# T: AC*
def aminoT(x, seq, pos):
    return And(x[seq][3*pos+0]['A']==1,
               x[seq][3*pos+1]['C']==1)
# W: UGG
def aminoW(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['U']==1,
                  x[seq][3*pos+1]['G']==1,
                  x[seq][3*pos+2]['G']==1))
# Y: UAU, UAC
def aminoY(x, seq, pos):
    return Or(And(x[seq][3*pos+0]['U']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['U']==1),
              And(x[seq][3*pos+0]['U']==1,
                  x[seq][3*pos+1]['A']==1,
                  x[seq][3*pos+2]['C']==1))
# V: GU*
def aminoV(x, seq, pos):
    return And(x[seq][3*pos+0]['G']==1,
               x[seq][3*pos+1]['U']==1)

amino2constraint = {}
amino2constraint['A'] = aminoA
amino2constraint['R'] = aminoR
amino2constraint['N'] = aminoN
amino2constraint['D'] = aminoD
amino2constraint['C'] = aminoC
amino2constraint['Q'] = aminoQ
amino2constraint['E'] = aminoE
amino2constraint['G'] = aminoG
amino2constraint['H'] = aminoH
amino2constraint['I'] = aminoI
amino2constraint['L'] = aminoL
amino2constraint['K'] = aminoK
amino2constraint['M'] = aminoM
amino2constraint['F'] = aminoF
amino2constraint['P'] = aminoP
amino2constraint['S'] = aminoS
amino2constraint['T'] = aminoT
amino2constraint['W'] = aminoW
amino2constraint['Y'] = aminoY
amino2constraint['V'] = aminoV
