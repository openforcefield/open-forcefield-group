import pandas as pd
import nmrpystar
import mdtraj as md

t0 = md.load("./Trajectories/1am7_1.dcd", top="./1am7_fixed.pdb")[0:50]
t1 = md.load("./Trajectories/1am7_1.dcd", top="./1am7_fixed.pdb")[-50:]

prediction0 = md.nmr.chemical_shifts_shiftx2(t0).mean(1)  # Average over time dimensions
prediction1 = md.nmr.chemical_shifts_shiftx2(t1).mean(1)  # Average over time dimensions

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

prediction0.name = "value"
prediction1.name = "value"

delta0 = (expt - prediction0).dropna()
rms0 = (delta0 ** 2.).reset_index().groupby("name").value.mean() ** 0.5

delta1 = (expt - prediction1).dropna()
rms1 = (delta1 ** 2.).reset_index().groupby("name").value.mean() ** 0.5
