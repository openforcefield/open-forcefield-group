import os
import itertools
import mdtraj
import pandas as pd
from patsy import dmatrices
import statsmodels.api as sm

def weighted_ave(x, tau=1.0):
    n = len(x)
    t = np.arange(n)
    p = np.exp(-t / tau)
    p /= p.sum()
    return (x.T * p).sum(1)

cos = np.cos
sin = np.sin
ave = lambda x: x.mean(0)
phi0 = np.deg2rad(0.0)

data = pd.read_csv("/home/kyleb/src/tjlane/scalar-couplings/kyleb/data/data.csv", index_col=0)
data["C1"] = np.nan
data["S1"] = np.nan
data["C2"] = np.nan
data["S2"] = np.nan
data["CS"] = np.nan

tau_dict = {"1UBQ":8000.0, "1BU5":110.}


for name in ["1UBQ", "1BU5"]:
    ave = lambda x: weighted_ave(x, tau=tau_dict[name])
    traj = mdtraj.load("/home/kyleb/src/tjlane/scalar-couplings/nmr_trajectories/%s_protein.dcd" % name, top="/home/kyleb/src/tjlane/scalar-couplings/nmr_trajectories/%s_protein.pdb" % name)
    rid, indices = mdtraj.geometry.atom_sequence_finder(traj, ["H","N", "CA", "HA"], residue_offsets=[0, 0, 0, 0])
    rid += 1  # Fix for mdtraj issues
    phi = mdtraj.geometry.dihedral.compute_dihedrals(traj, indices)
    indices = ["%s_%d" % (name, i) for i in rid]
    data = data.reindex(data.index.union(pd.Index(indices)))
    data.C1[indices] = ave(cos(phi + phi0))
    data.S1[indices] = ave(sin(phi + phi0))
    data.C2[indices] = ave(cos(phi + phi0) ** 2.)
    data.S2[indices] = ave(sin(phi + phi0) ** 2.)
    data.CS[indices] = ave(sin(phi + phi0) * cos(phi + phi0))


data = data.dropna(axis=0)

y, X = dmatrices('expt ~ C1 + C2', data=data, return_type='dataframe')

model = sm.OLS(y, X)
results = model.fit()
print results.summary()

yhat = results.predict()
delta = y.values[:,0] - yhat
rms = (delta ** 2.).mean() ** 0.5
rms
