import pdbbuilder

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]
amino_acids.remove("G")

sequences = ["%s%s" % ("G", aa) for aa in amino_acids]
sequences.extend(["%s%s" % (aa, "G") for aa in amino_acids])
cap_NT = "ACE"
cap_CT = "NH2"

for sequence in sequences:
    filename = "./pdbs/%s-%s-%s.pdb" % (cap_NT, sequence, cap_CT)
    print(filename)
    pdbbuilder.build_pdb(sequence, filename, True, cap_CT)



amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]
sequences = ["%s" % (aa) for aa in amino_acids]
cap_NT = "ACE"
cap_CT = "NME"

for sequence in sequences:
    filename = "./pdbs/%s-%s-%s.pdb" % (cap_NT, sequence, cap_CT)
    print(filename)
    pdbbuilder.build_pdb(sequence, filename, True, cap_CT)
