#!/usr/bin/env python

import os, sys
import numpy as np
import pandas as pd
import nmrpystar
import mdtraj as md
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-max', type=int, default=-1, help='Go up to a maximum frame (default is whole trajectory)')
args, sys.argv = parser.parse_known_args(sys.argv)

# Run compare_shifts.py [PDB] [XTC] [12345.str]

t = md.load(sys.argv[2], top=sys.argv[1])
outfnm = sys.argv[4]
# if os.path.exists(outfnm): sys.exit()

parsed = nmrpystar.parse(open(sys.argv[3]).read())
print(parsed.status)

q = parsed.value.saves["assigned_chem_shift_list_1"].loops[1]
x = pd.DataFrame(q.rows, columns=q.keys)
x = x[["Atom_chem_shift.Seq_ID", "Atom_chem_shift.Atom_ID", "Atom_chem_shift.Val"]]
x.rename(columns={"Atom_chem_shift.Seq_ID":"resSeq", "Atom_chem_shift.Atom_ID":"name", "Atom_chem_shift.Val":"value"}, inplace=True)

# Need to make dtypes match to do eventual comparison.
x["resSeq"] = x["resSeq"].astype('int')
x["value"] = x["value"].astype('float')

expn = x.set_index(["resSeq", "name"])
expt = x.set_index(["resSeq", "name"]).value

print "Doing ShiftX2 prediction."
outd = sys.argv[2]+'.shiftx2'

if not os.path.exists(outd): 
    os.makedirs(outd)
os.chmod(outd, 0755)
if not os.path.exists(os.path.join(outd, 'expn.csv')): 
    expn.to_csv(os.path.join(outd, 'expn.csv'))
np.savetxt(os.path.join(outd, 'expt.txt'), expt)

predictions = []
for i, ti in enumerate(t):
    if i == args.max: break
    fnm = os.path.join(outd, '%i.txt' % i)
    if os.path.exists(fnm):
        os.chmod(fnm, 0644)
        print "Loading ShiftX2 prediction for frame %i." % i
        prediction_i = pd.DataFrame.from_csv(fnm, index_col=[0,1])
    else:
        print "Computing ShiftX2 prediction for frame %i." % i
        if not os.path.exists(outd): 
            os.makedirs(outd)
        os.chmod(outd, 0755)
        prediction_i = md.nmr.chemical_shifts_shiftx2(ti)
        prediction_i.to_csv(fnm)
        os.chmod(fnm, 0644)
    predictions.append(prediction_i)

prediction = pd.concat(predictions).mean(1)
print "ShiftX2 done."

delta = (expt - prediction).dropna()
delta.name = "value"

# print "ShiftX2 prediction saved to file."
# delta.to_csv('delta.csv')

rms = (delta ** 2.).reset_index().groupby("name").value.mean() ** 0.5

with open(outfnm, 'w') as f: print >> f, rms

