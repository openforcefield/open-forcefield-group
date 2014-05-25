#!/usr/bin/env python

import os
import re
from forcebalance.molecule import *
from forcebalance.nifty import _exec

AA = re.match('^[A-Z]*', os.path.split(os.getcwd())[-1]).group(0)

P = Molecule("/home/leeping/projects/VSP27-Protein/Dihedrals/ACE-X-NME-1/%s/unmangled.pdb" % AA)

M = Molecule("all.arc")

M.resid = P.resid
M.resname = P.resname
M.atomname = P.atomname

for i in range(M.na):
    if M.resname[i] in ['ACE', 'NME']:
        if M.atomname[i] == 'H1':
            M.atomname[i] = 'HH31'
        elif M.atomname[i] == 'H2':
            M.atomname[i] = 'HH32'
        elif M.atomname[i] == 'H3':
            M.atomname[i] = 'HH33'
    if M.resname[i] == 'ILE' and 'HD1' in M.atomname[i]:
        M.atomname[i] = M.atomname[i].replace('HD1','HD')
    if M.resname[i] == 'ILE' and M.atomname[i] == 'CD1':
        M.atomname[i] = 'CD'

hletts = sorted(list(set([re.sub('[0-9]$', '', i) for i in M.atomname if re.match('^H[A-Z][0-9]', i)])))
for hlett in hletts:
    if (hlett+'1') not in M.atomname and (hlett+'3') in M.atomname:
        for i in range(M.na):
            if M.resname[i] == 'TRP' and 'HZ' in M.atomname[i]:
                continue
            elif M.atomname[i] == (hlett+'2'):
                M.atomname[i] = hlett+'1'
            elif M.atomname[i] == (hlett+'3'):
                M.atomname[i] = hlett+'2'

M.boxes = [CubicLattice(30)]
M[0].write("one.pdb")
_exec("pdb2gmx -f one.pdb", stdin="5\n6\n")
N = Molecule("conf.gro")
M.reorder_according_to(N)
M.write("all.gro")

_exec("mv qdata.txt qdata0.txt")
Q = Molecule("qdata0.txt")
Q.reorder_according_to(N)
Q.write("qdata.txt")

_exec("mv topol.top topol0.top")

fout = open('topol.top','w')

do_print = 0
for line in open('topol0.top').readlines():
    if line.startswith("#"): continue
    line = re.sub(';.*$', '', line.replace('\n',''))
    match = re.match('^ *\[ *([^ ]*) *\]', line)
    if match:
        print >> fout
        content = match.group(1)
        if content == 'moleculetype':
            do_print = 1
    if do_print and len(line) > 0:
        print >> fout, line

fout.close()
