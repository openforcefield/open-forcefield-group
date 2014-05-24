import pandas as pd
import nmrpystar
import mdtraj as md

stride = 100

t0 = md.load(["./Trajectories_ff99sbnmr/1am7_%d.dcd" % i for i in range(10)], top="./1am7_fixed.pdb")[::stride]
t1 = md.load(["./Trajectories/1am7_%d.dcd" % i for i in range(15)], top="./1am7_fixed.pdb")[::stride]


#full_prediction0 = md.nmr.chemical_shifts_shiftx2(t0)
#full_prediction1 = md.nmr.chemical_shifts_shiftx2(t1)

#full_prediction0 = md.nmr.chemical_shifts_spartaplus(t0)
#full_prediction1 = md.nmr.chemical_shifts_spartaplus(t1)

full_prediction0 = md.nmr.chemical_shifts_ppm(t0)
full_prediction1 = md.nmr.chemical_shifts_ppm(t1)



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



prediction0 = full_prediction0.mean(1)  # Average over time dimensions
prediction0.name = "value"

prediction1 = full_prediction1.mean(1)  # Average over time dimensions
prediction1.name = "value"

sigma = pd.Series({"C":0.8699, "CA":0.7743, "H":0.3783, "HA":0.1967, "HA2":0.1967, "HA3":0.1967, "N":2.0862})
#sigma = pd.Series({"C":0.8699, "CA":0.7743, "H":0.3783, "HA":0.1967, "N":2.0862})
sigma.name = "value"

z0 = ((expt - prediction0)).dropna()
z1 = ((expt - prediction1)).dropna()

rms0 = (z0 ** 2.).reset_index().groupby("name").value.mean() ** 0.5
rms1 = (z1 ** 2.).reset_index().groupby("name").value.mean() ** 0.5

rms0 = rms0 / sigma
rms1 = rms1 / sigma

rms0 / rms1


