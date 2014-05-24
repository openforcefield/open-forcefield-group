import pdbbuilder

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]
sequences = ["%s" % (aa) for aa in amino_acids]
cap_NT = "ACE"
cap_CT = "NME"

for sequence in sequences:
    filename = "./pdbs/%s_%s_%s.pdb" % (cap_NT, sequence, cap_CT)
    print(filename)
    pdbbuilder.build_pdb(sequence, filename, True, cap_CT)
