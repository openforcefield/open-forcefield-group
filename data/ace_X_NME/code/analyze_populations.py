import pandas as pd
import mdtraj as md
from dipeptide_parameters import *

reference = pd.read_csv("./experimental_data/baldwin_table1_populations.csv", index_col=0)

data = []
for (ff, water, seq) in products:
    try:
        aa = seq.split("_")[1]
        t = md.load("./dcd/%s_%s_%s.dcd" % (ff, water, seq), top="./pdbs/%s.pdb" % (seq))
    except:
        continue
    phi = md.compute_phi(t)[1][:, 0] * 180 / np.pi
    psi = md.compute_psi(t)[1][:, 0] * 180 / np.pi
    ass = assign(phi, psi)
    populations = pd.Series({"PPII":0.0, "beta":0.0, "alpha":0.0, "other":0.0})
    populations += ass.value_counts(normalize=True)
    data.append([ff, water, aa, populations["PPII"], populations["beta"], populations["alpha"]])

data = pd.DataFrame(data, columns=["ff", "water", "aa", "PPII", "beta", "alpha"])

