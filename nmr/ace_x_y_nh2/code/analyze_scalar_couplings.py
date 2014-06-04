import pandas as pd
import mdtraj as md
from ace_x_y_nh2_parameters import *

larger = pd.read_csv("./data/larger_couplings.csv")
smaller = pd.read_csv("./data/smaller_couplings.csv")

reference = []
for aa in amino_acids:
    
    value = smaller.ix["G"][aa]
    xyz = ["G%s" % aa, 0, value]
    reference.append(xyz)
    
    value = larger.ix["G"][aa]
    xyz = ["G%s" % aa, 1, value]
    reference.append(xyz)
    
    value = larger.ix[aa]["G"]
    xyz = ["%sG" % aa, 0, value]
    reference.append(xyz)
    
    value = smaller.ix[aa]["G"]
    xyz = ["%sG" % aa, 1, value]
    reference.append(xyz)

reference = pd.DataFrame(reference, columns=["seq", "resSeq", "value"])
reference = reference.set_index(["seq", "resSeq"]).value
reference = reference.drop_duplicates()


data = []
for (ff, water, seq) in products:
    try:
        aa0, aa1 = seq.split("_")[1]
        aa_string = "%s%s" % (aa0, aa1)
        t = md.load("./dcd/%s_%s_%s.dcd" % (ff, water, seq), top="./pdbs/%s.pdb" % (seq))[1500:]
    except:
        continue
    #phi = md.compute_phi(t)[1] * 180 / np.pi
    #J0, J1 = scalar_couplings.J3_HN_HA(phi).mean(0)
    J0, J1 = md.compute_J3_HN_HA(t)[1].mean(0)
    data.append([ff, water, aa_string, 0, J0])
    data.append([ff, water, aa_string, 1, J1])

data = pd.DataFrame(data, columns=["ff", "water", "seq", "resSeq", "value"])


X = data.pivot_table(cols=["seq", "resSeq"], rows=["ff", "water"], values="value")
delta = X - reference
Z = (delta / 0.36)
rms_by_model = (Z ** 2.).mean(1) ** 0.5
rms_by_model

