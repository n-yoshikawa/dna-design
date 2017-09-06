import dna_design
import constants as c


sequence_list = ["GAYY", "GAHY", "GAFY", "GAYH", "GGYY", "GGHY", "GAHH", "GAFH", "GGFY", "GGYH", "GAYW", "GGHH", "GAHF", "GAFF"]
result_base, generated_aminos = dna_design.optimize(sequence_list, 2)

print "Optimal sequence"
for b in result_base:
    print b

print "Generated amino sequences"
for g in generated_aminos:
    print g, "-> in list" if g in sequence_list else "-> Not in list"

