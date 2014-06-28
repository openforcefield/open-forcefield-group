#!/usr/bin/env python

import os, sys
import numpy as np
import pandas as pd
import nmrpystar
import mdtraj as md

# Run compare_shifts_cache.py [PDB] [XTC] [12345.str]
# The purpose of this script is to cache the NMR predictions for each frame to disk,
# in order to save on computer time.

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

if not os.path.exists(outd): os.makedirs(outd)
expn.to_csv(os.path.join(outd, 'expn.csv'))
np.savetxt(os.path.join(outd, 'expt.txt'), expt)

predictions = []
for i, ti in enumerate(t):
    fnm = os.path.join(outd, '%i.txt' % i)
    if os.path.exists(fnm):
        print "Loading ShiftX2 prediction for frame %i." % i
        prediction_i = pd.DataFrame.from_csv(fnm, index_col=[0,1])
    else:
        print "Computing ShiftX2 prediction for frame %i." % i
        if not os.path.exists(outd): os.makedirs(outd)
        prediction_i = md.nmr.chemical_shifts_shiftx2(ti)
        prediction_i.to_csv(fnm)
    predictions.append(prediction_i)

prediction = pd.concat(predictions).mean(1)
print "ShiftX2 done."

delta = (expt - prediction).dropna()
delta.name = "value"

# print "ShiftX2 prediction saved to file."
# delta.to_csv('delta.csv')

rms = (delta ** 2.).reset_index().groupby("name").value.mean() ** 0.5

with open(outfnm, 'w') as f: print >> f, rms

