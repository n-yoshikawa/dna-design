import functools
import itertools
import math
import operator
import time

from z3 import *

from constants import *
from constraints import *

def prod(iterable):
    return functools.reduce(operator.mul, iterable)

def optimize(amino_list, num_of_seq):
    # check if all sequence length is the same
    seq_len = len(amino_list[0])
    for seq in amino_list:
        assert len(seq) == seq_len

    base_length = 3 * seq_len

    # variables to describe constraints
    bases = ['A', 'C', 'G', 'U']
    x = [[{base: Int("x-{}-{}-{}".format(seq, pos, base)) for base in bases} for pos in range(base_length)] for seq in range(num_of_seq)]
    epsilon = Int("epsilon")
    delta = [Int("delta-{}".format(i)) for i in range(num_of_seq)]
    c = Int("c")
    size = Int("size")


    opt = Optimize()
    # all variable is 0 or 1
    for primer in range(num_of_seq):
        for pos in range(base_length):
            for base in bases:
                opt.add(x[primer][pos][base] >= 0)
                opt.add(x[primer][pos][base] <= 1)

    # nonnegative value
    opt.add(epsilon >= 0)
    for i in range(num_of_seq):
        opt.add(delta[i] >= 0)
    opt.add(c == 1)


    # helper function: all amino acids in the given sequence will be synthesized
    def seq2constraint(s, p):
        return And([amino2constraint[s[i]](x, p, i) for i in range(seq_len)])

    # all amino acid sequence will be synthesized
    for seq in amino_list:
        opt.add(Or([seq2constraint(seq, p) for p in range(num_of_seq)]))

    # gap is not allowed
    for i in range(num_of_seq):
        for j in range(base_length):
            opt.add(sum([x[i][j][b] for b in bases]) >= 1)

    # constraint of number of random base per primer
    for i in range(num_of_seq):
        opt.add(sum([x[i][j][b] for j in range(base_length) for b in bases]) <= epsilon + delta[i])

    # object function: smaller ambiguity is better solution
    cost = Int('cost')
    opt.add(cost == epsilon + c * sum([delta[i] for i in range(num_of_seq)]))
    #opt.add(cost==sum([x[i][j][b] for b in bases for j in range(base_length) for i in range(num_of_seq)]))

    # total number of generated amino acids
    opt.add(size == sum([prod([sum([x[i][j][b] for b in bases]) for j in range(base_length)]) for i in range(num_of_seq)]))

    # Enumerate amino acid sequences from given nucleotide sequence
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

    print('Minimizing...')
    h = opt.minimize(cost)
    if opt.check() == sat:
        initial_model = opt.model()
        initial_size = int("{}".format(initial_model[size]))
        initial_cost = int("{}".format(initial_model[cost]))

    best_size = initial_size
    best_model = initial_model

    begin_time = time.time()
    count = 0
    print("count,time,cost,epsilon,delta,size")
    while opt.check() == sat:
        count += 1
        model = opt.model()
        current_size = int("{}".format(model[size]))
        if current_size < best_size:
            best_size = current_size
            best_model = model

        result_base = []
        for n in range(num_of_seq):
            result = "".join([dict2base(model, x[n][pos]) for pos in range(base_length)])
            result_base.append(result)
        generated_aminos = []
        for r in result_base:
            generated_aminos.append(generated_amino(split_n(r, 3)))

        elapsed_time = time.time() - begin_time
        print("{},{},{},{},{},{}".format(count, elapsed_time, int("{}".format(model[cost])), int("{}".format(model[epsilon])),
                                      [int("{}".format(model[delta[i]])) for i in range(num_of_seq)], int("{}".format(model[size]))))

        opt.add(Or([x[primer][pos][base] != opt.model()[x[primer][pos][base]] for base in bases for pos in range(base_length) for primer in range(num_of_seq)]))

    model = best_model
    result_base = []
    for n in range(num_of_seq):
        result = "".join([dict2base(model, x[n][pos]) for pos in range(base_length)])
        result_base.append(result)
    generated_aminos = []
    for r in result_base:
        generated_aminos.append(generated_amino(split_n(r, 3)))

    print("cost: {}".format(int("{}".format(model[cost]))))
    print("epsilon: {}".format(int("{}".format(model[epsilon]))))
    print("delta: {}".format([int("{}".format(model[delta[i]])) for i in range(num_of_seq)]))
    print("size: {}".format(int("{}".format(model[size]))))
    print(result_base, generated_aminos)

    return (result_base, generated_aminos)
