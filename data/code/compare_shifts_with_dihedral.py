import pandas as pd
import nmrpystar
import mdtraj as md

#t = md.load("./1am7_fixed.pdb")
t = md.load(["./Trajectories/1am7_%d.dcd" % i for i in range(15)], top="./1am7_fixed.pdb")[::50]
psi = md.compute_psi(t)[1][:, 48]  * 180 / pi


#full_prediction = md.nmr.chemical_shifts_shiftx2(t)
full_prediction = md.nmr.chemical_shifts_spartaplus(t)


parsed = nmrpystar.parse(open("./16664.str").read())
print(parsed.status)

q = parsed.value.saves["assigned_chem_shift_list_1"].loops[1]
x = pd.DataFrame(q.rows, columns=q.keys)
x = x[["Atom_chem_shift.Seq_ID", "Atom_chem_shift.Atom_ID", "Atom_chem_shift.Val"]]
x.rename(columns={"Atom_chem_shift.Seq_ID":"resSeq", "Atom_chem_shift.Atom_ID":"name", "Atom_chem_shift.Val":"value"}, inplace=True)

# Need to make dtypes match to do eventual comparison.
x["resSeq"] = x["resSeq"].astype('int')
x["value"] = x["value"].astype('float')
expt = x.set_index(["resSeq", "name"]).value


cond = psi > 50
prediction0 = full_prediction[where(cond)[0]].mean(1)  # Average over time dimensions
prediction0.name = "value"

prediction1 = full_prediction[where(~cond)[0]].mean(1)  # Average over time dimensions
prediction1.name = "value"


delta0 = (expt - prediction0).dropna()
delta1 = (expt - prediction1).dropna()

rms0 = (delta0 ** 2.).reset_index().groupby("name").value.mean() ** 0.5
rms1 = (delta1 ** 2.).reset_index().groupby("name").value.mean() ** 0.5

rms0 - rms1


local_residues = np.array([46, 47, 48, 49, 50, 51, 52, 86, 42, 36, 64, 65])
residue_cond = np.where(np.in1d(full_prediction.index.get_level_values(0), local_residues))[0]

prediction0 = full_prediction[where(cond)[0]].iloc[residue_cond].mean(1)  # Average over time dimensions
prediction0.name = "value"

prediction1 = full_prediction[where(~cond)[0]].iloc[residue_cond].mean(1)  # Average over time dimensions
prediction1.name = "value"


delta0 = (expt - prediction0).dropna()
delta1 = (expt - prediction1).dropna()

rms0 = (delta0 ** 2.).reset_index().groupby("name").value.mean() ** 0.5
rms1 = (delta1 ** 2.).reset_index().groupby("name").value.mean() ** 0.5

rms0 - rms1
