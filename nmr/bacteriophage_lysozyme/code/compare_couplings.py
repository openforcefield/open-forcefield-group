import pandas as pd
import nmrpystar
import mdtraj as md

t = md.load("/home/kyleb/dat/tmp/pro_0.xtc", top="/home/kyleb/dat/tmp/pro.pdb")

parsed = nmrpystar.parse(open("./19127.str").read())
print(parsed.status)

q = parsed.value.saves["coupling_constant_list_1"].loops[1]
x = pd.DataFrame(q.rows, columns=q.keys)
x = x[["Coupling_constant.Seq_ID_1", "Coupling_constant.Val", "Coupling_constant.Val_err"]]
x.rename(columns={"Coupling_constant.Seq_ID_1":"resSeq", "Coupling_constant.Val":"value", "Coupling_constant.Val_err":"err"}, inplace=True)

# Need to make dtypes match to do eventual comparison.
x["resSeq"] = x["resSeq"].astype('int')
x["value"] = x["value"].astype('float')

expt = x.set_index(["resSeq"]).value

top, bonds = t.top.to_dataframe()
ind, values = md.compute_J3_HN_HA(t)
prediction = pd.Series(values.mean(0), top.ix[ind[:,-1]].resSeq)


delta = (expt - prediction).dropna()
delta.name = "value"

rms = (delta ** 2.).mean(0) ** 0.5
