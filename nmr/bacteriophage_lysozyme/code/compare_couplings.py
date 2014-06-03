import pandas as pd
import nmrpystar
import mdtraj as md

t = md.load("/home/kyleb/dat/tmp/pro_0.xtc", top="/home/kyleb/dat/tmp/pro.pdb")

parsed = nmrpystar.parse(open("./19127.str").read())
print(parsed.status)

q = parsed.value.saves["coupling_constant_list_1"].loops[1]
x = pd.DataFrame(q.rows, columns=q.keys)
x = x[["Coupling_constant.Seq_ID_1", "Coupling_constant.Seq_ID_2", "Coupling_constant.Val", "Coupling_constant.Val_err"]]
x.rename(columns={"Coupling_constant.Seq_ID_1":"resSeq1", 
"Coupling_constant.Seq_ID_2":"resSeq2", "Coupling_constant.Val":"value", "Coupling_constant.Val_err":"err"},
 inplace=True)

# Need to make dtypes match to do eventual comparison.
x["resSeq"] = x["resSeq"].astype('int')
x["value"] = x["value"].astype('float')

expt = x.set_index(["resSeq", "name"]).value

delta = (expt - prediction).dropna()
delta.name = "value"

rms = (delta ** 2.).reset_index().groupby("name").value.mean() ** 0.5
