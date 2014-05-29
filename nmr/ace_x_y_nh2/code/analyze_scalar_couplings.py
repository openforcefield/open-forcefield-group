import pandas as pd
import mdtraj as md
from ace_x_y_nh2_parameters import *
import scalar_couplings

larger = pd.read_csv("./data/larger_couplings.csv")
smaller = pd.read_csv("./data/smaller_couplings.csv")

data = []
for (ff, water, seq) in products:
    try:
        aa0, aa1 = seq.split("_")[1]
        t = md.load("./dcd/%s_%s_%s.dcd" % (ff, water, seq), top="./pdbs/%s.pdb" % (seq))
    except:
        continue
    phi = md.compute_phi(t)[1] * 180 / np.pi
    J0, J1 = scalar_couplings.J3_HN_HA(phi).mean(0)
    data.append([ff, water, aa0, aa1, J0, J1])

data = pd.DataFrame(data, columns=["ff", "water", "AA0", "AA1", "J0", "J1"])

X = data.pivot_table(values="J", cols=["AA"], rows=["ff", "water"])
delta = X - reference
Z = (delta / 0.36)
rms_by_model = (Z ** 2.).mean(1) ** 0.5
rms_by_model
