import pandas as pd
import nmrpystar
import mdtraj as md

t = md.load("./1am7.dcd", top="./1am7_fixed.pdb")

prediction = md.nmr.chemical_shifts_shiftx2(t).mean(1)  # Average over time dimensions
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

delta = (expt - prediction).dropna()

rms = (delta ** 2.).reset_index().groupby("name").value.mean() ** 0.5
