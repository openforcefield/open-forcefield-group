#!/usr/bin/env python

from forcebalance.molecule import Molecule
from forcebalance.nifty import _exec
import os

np=""
if "CORES_PER_WORKER" in os.environ and int(os.environ["CORES_PER_WORKER"]) > 1:
    np=" -np %i" % int(os.environ["CORES_PER_WORKER"])

_exec("touch opt.xyz")
_exec("touch energy.txt")
_exec("rm -f qchem.out.prev")
_exec("touch qchem.out.prev")
qcin = Molecule("qchem.in", ftype="qcin")
qcin.edit_qcrems({'geom_opt_max_cycles':'100'})
qcin.write("qchem.in")
_exec("qchem42 %s qchem.in qchem.out &> qchem.err" % np)

def special_criterion():
    mk = 0
    mx = 0
    Cnvgd = {}
    for ln, line in enumerate(open("qchem.out").readlines()):
        if "Maximum optimization cycles reached" in line:
            mx = 1
        if "Maximum     Tolerance    Cnvgd?" in line:
            mk = ln
        if mk > 0 and ln > mk and ln <= mk+3:
            s = line.split()
            try:
                Cnvgd[s[0]] = float(s[-3])
            except: pass
    if mx and len(Cnvgd) == 3:
        if Cnvgd['Gradient'] < 0.001000 and Cnvgd['Displacement'] < 0.001200 and Cnvgd['Energy'] < 0.000001:
            return 1
    return 0

Attempt = 1
message = "Optimization with standard settings."
errok = []
while True:
    _exec("echo '%s' >> qchem.out" % message)
    try:
        M = Molecule("qchem.out", errok=errok)
        break
    except:
        if special_criterion():
            M = Molecule("qchem.out", errok=["Maximum optimization cycles reached"])
            break
        Attempt += 1
        qcin = Molecule("qchem.in", ftype="qcin")
        if Attempt == 2:
            qcin.edit_qcrems({'geom_opt_dmax':'100', 'geom_opt_max_cycles':'300'})
            message = "Attempt 2: Backup optimization with reduced step."
        elif Attempt == 3:
            qcin.edit_qcrems({'geom_opt_dmax':None, 'geom_opt_max_cycles':'100', 'geom_opt_update':'5'})
            message = "Attempt 3: Using Hessian update 5."
        elif Attempt == 4:
            qcin.edit_qcrems({'geom_opt_dmax':'100', 'geom_opt_max_cycles':'300'})
            message = "Attempt 4: Using Hessian update 5 with reduced step."
        elif Attempt == 5:
            qcin.edit_qcrems({'geom_opt_dmax':None, 'geom_opt_max_cycles':'100', 'geom_opt_update':'3'})
            message = "Attempt 5: Using Hessian update 3."
        elif Attempt == 6:
            qcin.edit_qcrems({'geom_opt_dmax':'100', 'geom_opt_max_cycles':'300'})
            message = "Attempt 6: Using Hessian update 3 with reduced step."
        elif Attempt == 7:
            qcin.edit_qcrems({'geom_opt_update':None, 'geom_opt_coords':'0'})
            message = "Attempt 7: Last-ditch optimization with Cartesian coordinates."
        elif Attempt == 8:
            errok = ["Maximum optimization cycles reached", "OPTIMIZE fatal error"]
            message = "Cartesian coordinate optimization has failed miserably!"
        if Attempt < 8:
            qcin.write("qchem.in")
            _exec("cat qchem.out >> qchem.out.prev")
            _exec("qchem42 %s qchem.in qchem.out 2> qchem.err" % np)
        if Attempt > 8: break # It should never do this

_exec("mv qchem.out qchem.out_")
_exec("cat qchem.out.prev qchem.out_ > qchem.out")
_exec("bzip2 qchem.out")
    
M[-1].write("opt.xyz")
E = M.qm_energies[-1]
with open("energy.txt","w") as f: print >> f, "% .10f" % E
