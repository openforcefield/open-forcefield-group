import pandas as pd
import mdtraj as md
from dipeptide_parameters import assign, amino_acids

water_name = "tip3p"

data = {}
for aa in amino_acids:
    t = md.load("./dcd/amber96_%s_ACE_%s_NME.dcd" % (water_name, aa), top="./pdbs/ACE_%s_NME.pdb" % (aa))
    phi = md.compute_phi(t)[1][:, 0] * 180 / np.pi
    psi = md.compute_psi(t)[1][:, 0] * 180 / np.pi
    ass = assign(phi, psi)
    populations = ass.value_counts(normalize=True)
    populations.name = aa
    data[aa] = populations

data = pd.DataFrame(data).T
reference = pd.read_csv("./population_data/baldwin_table1.csv", index_col=0)
