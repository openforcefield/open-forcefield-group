import pdbbuilder
from ace_x_y_nh2_parameters import sequences

cap_NT = "ACE"
cap_CT = "NH2"

for sequence in sequences:
    filename = "./pdbs/%s.pdb" % sequence
    sequence = sequence.split("_")[1]
    print(filename)
    pdbbuilder.build_pdb(sequence, filename, True, cap_CT)
