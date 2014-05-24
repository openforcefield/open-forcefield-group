import os
import itertools
import mdtraj
import pandas as pd

cos = np.cos
sin = np.sin
ave = lambda x: x.mean(0).mean(0)
phi0 = np.deg2rad(0.0)

amino_acids = ["A" , "C" , "D" , "E" , "F" , "G" , "H" , "I" , "K" , "L" , "M" , "N" , "Q" , "R" , "S" , "T" , "V" , "W" , "Y"]
labels = ["%s%s" % (a0,a1) for (a0, a1) in itertools.product(amino_acids, repeat=2)]

bad_pairs = ["AI","AY","CV","CY","DW","EF","EW","FA","FC","FI","FM","FN","FQ","FT","FY","IF","IM","IT","IV","IW","IY","LF","LI","LV","LW","MF","MI","ML","MM","MY","NY","QF","QY","SF","SI","TF","TI","TW","VF","VI","VV","VW","WA","WC","WF","WI","WL","WN","WV","YF","YI","YL","YM","YN","YV"]
bad_pairs.extend(["HH", "TH", "KH"])
bad_pairs.extend(["DI", "DK"])
bad_pairs.extend(["AD","CD","DD","FD","KD","LD","ND","QD","RD","TD","VD","WD","YD","EC","EH","SE","LE"])
bad_pairs.extend(["G%s" % aa for aa in amino_acids])
bad_pairs.extend(["%sG" % aa for aa in amino_acids])
labels = list(set(labels).difference(set(bad_pairs)))

small = pd.read_csv("/home/kyleb/src/tjlane/scalar-couplings/kyleb/smaller_couplings.csv", index_col=0)
large = pd.read_csv("/home/kyleb/src/tjlane/scalar-couplings/kyleb/larger_couplings.csv", index_col=0)
averaged = 0.5 * (small + large)

data = pd.DataFrame(index=labels, columns=["expt", "C1", "C2", "S1", "S2", "CS"], dtype='float')

for label in labels:
    a0, a1 = label
    if not os.path.exists("/home/kyleb/dat/peptides/dcd/%s-capped.dcd" % (label)):
        continue
    traj = mdtraj.load("/home/kyleb/dat/peptides/dcd/%s-capped.dcd" % (label), top="/home/kyleb/dat/peptides/raw/%s-capped.pdb" % (label))
    rid, indices = mdtraj.geometry.atom_sequence_finder(traj, ["H","N", "CA", "HA"], residue_offsets=[0, 0, 0, 0])
    phi = mdtraj.geometry.dihedral.compute_dihedrals(traj, indices)
    phi = mdtraj.geometry.dihedral.compute_phi(traj)[1]
    data["C1"][label] = ave(cos(phi + phi0))
    data["C2"][label] = ave(cos(phi + phi0) ** 2.)
    data["S1"][label] = ave(sin(phi + phi0))
    data["S2"][label] = ave(sin(phi + phi0) ** 2.)
    data["CS"][label] = ave(sin(phi + phi0) * cos(phi + phi0))
    data["expt"][label] = averaged[a1][a0]

data = data.dropna(axis=0)

y, X = dmatrices('expt ~ C1 + C2 + S1', data=data, return_type='dataframe')

model = sm.OLS(y, X)
results = model.fit()
print results.summary()

data["yhat"] = results.predict()
data["delta"] = data.expt - data.yhat
rms = (data.delta ** 2.).mean() ** 0.5
rms

