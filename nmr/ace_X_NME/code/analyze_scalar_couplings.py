import pandas as pd
import mdtraj as md
from dipeptide_parameters import *
import scalar_couplings

reference = pd.read_csv("./experimental_data/baldwin_table1_2006_couplings.csv", index_col=0)
reference = reference["coupling"]

data = []
for (ff, water, seq) in products:
    try:
        aa = seq.split("_")[1]
        t = md.load("./dcd/%s_%s_%s.dcd" % (ff, water, seq), top="./pdbs/%s.pdb" % (seq))
    except:
        continue
    phi = md.compute_phi(t)[1][:, 0] * 180 / np.pi
    psi = md.compute_psi(t)[1][:, 0] * 180 / np.pi
    J = scalar_couplings.J3_HN_HA(phi).mean()
    data.append([ff, water, aa, J])

data = pd.DataFrame(data, columns=["ff", "water", "aa", "J"])

