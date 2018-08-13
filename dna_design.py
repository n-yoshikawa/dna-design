import math
import itertools

from z3 import *
from constants import *
from constraints import *

def optimize(amino_list, num_of_seq):
    # check if all sequence length is the same
    seq_len = len(amino_list[0])
    for seq in amino_list:
        assert len(seq) == seq_len

    base_length = 3 * seq_len

    # variables to describe constraints
    bases = ['A', 'C', 'G', 'U']
    x = [[{base: Int("x-{}-{}-{}".format(seq, pos, base)) for base in bases} for pos in range(base_length)] for seq in range(num_of_seq)]


    opt = Optimize()
    # all variable is 0 or 1
    for primer in range(num_of_seq):
        for pos in range(base_length):
            for base in bases:
                opt.add(x[primer][pos][base] >= 0)
                opt.add(x[primer][pos][base] <= 1)

    # gap is not allowed
    for i in range(num_of_seq):
        for j in range(base_length):
            opt.add(sum([x[i][j][b] for b in bases]) >= 1)

    # helper function: whole amino acid sequence will be synthesized
    def seq2constraint(s, p):
        return And([amino2constraint[s[i]](x, p, i) for i in range(seq_len)])

    # all amino acid sequence will be synthesized
    for seq in amino_list:
        opt.add(Or([seq2constraint(seq, p) for p in range(num_of_seq)]))

    # object function: smaller ambiguity is better solution
    cost = Int('cost')
    opt.add(cost==sum([x[i][j][b] for b in bases for j in range(base_length) for i in range(num_of_seq)]))

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


    def split_n(text, n):
        return [text[i*n:i*n+n] for i in range(len(text)//n)]

    h = opt.minimize(cost)
    if opt.check() == sat:
        model = opt.model()
        best_cost = int("{}".format(model[cost]))
    else:
        best_cost = -100
    best_expected = float("inf")
    best_base = []
    best_aminos = []
    count = 0

    while opt.check() == sat and count < 20:
        count += 1
        model = opt.model()
        current_cost = int("{}".format(model[cost]))
        if current_cost > best_cost:
            break
        result_base = []
        generated_aminos = []
        for n in range(num_of_seq):
            result = "".join([dict2base(model, x[n][pos]) for pos in range(base_length)])
            result_base.append(result)

        for r in result_base:
            generated_aminos.append(generated_amino(split_n(r, 3)))

        print(result_base, generated_aminos)
        expected = 0
        for sequence, aminos in zip(result_base, generated_aminos):
            n = 0
            for target in amino_list:
                if target in aminos:
                    n += 1
            expected += len(aminos) * math.log(n+0.01)

        if expected < best_expected:
            best_expected = expected
            best_base = result_base
            best_aminos = generated_aminos

        # add condition to get different solution
        opt.add(Or([x[primer][pos][base] != opt.model()[x[primer][pos][base]] for base in bases for pos in range(base_length) for primer in range(num_of_seq)]))
        h = opt.minimize(cost)

    return (best_base, best_aminos)
