import dna_design
import constants as c

sequence_list = ["SS","SG","SR", "GS", "GG", "GR"]
result_base, generated_aminos = dna_design.optimize(sequence_list, 1)

for b, seqs in zip(result_base, generated_aminos):
    print(b)
    for s in seqs:
        print(s, "-> in list" if s in sequence_list else "-> Not in list")
