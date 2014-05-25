#!/usr/bin/env python

import os, sys, re
import argparse
from forcebalance.molecule import Molecule
from collections import defaultdict, OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument('-bd', action="store_true", help='Label bond parameters')
parser.add_argument('-an', action="store_true", help='Label angle parameters')
parser.add_argument('-bk', action="store_true", help='Label backbone torsion angles')
parser.add_argument('-sd', action="store_true", help='Label sidechain torsion angles')
parser.add_argument('-eq', action="store_true", help='Label equilibrium positions in addition to force constants')
parser.add_argument('-sel', type=str, nargs='+', help='Activate parameters for one or more selected residues')

args, sys.argv= parser.parse_known_args(sys.argv)

llist = []
if args.sel != None:
    llist += [i.upper() for i in args.sel]
if args.bd:
    llist.append('bd')
if args.an:
    llist.append('an')
if args.bk:
    llist.append('bk')
if args.sd:
    llist.append('sd')
if args.eq:
    llist.append('eq')
if len(llist) == 0:
    print 'Did not specify any parameters to fit, use -h argument for help'
if len(sys.argv) != 2:
    raise RuntimeError('Did not specify force field file, use -h argument for help')
label = '-'+'-'.join(llist)

# Highlight pertinent dihedral interactions for amino acids in AMBER force field.

AACodes = []
AADirs = []
for d in sorted(os.listdir('../targets/')):
    if 'phipsi' not in d: continue
    if os.path.isdir('../targets/%s' % d):
        AACode = re.match('^[A-Z]*', d).group(0)
        if AACode not in AACodes:
            if args.sel != None and AACode not in [i.upper() for i in args.sel]: continue
            AACodes.append(AACode)
            AADirs.append('../targets/%s' % d)

# Get a mapping from atom names to atom classes for a molecule name.
def build_atomname_to_atomclass(mnm):
    atomgrp = 0
    content = ''
    this_mnm = ''
    anac = OrderedDict()
    for line in open('aminoacids.rtp').readlines():
        match = re.match('^ *\[ *([^ ]*) *\]', line)
        if match:
            old_content = content
            content = match.group(1)
            if content == 'atoms':
                this_mnm = old_content
                atomgrp = 1
            else:
                atomgrp = 0
        elif atomgrp and mnm == this_mnm:
            s = line.split()
            anac[s[0]] = s[1]
    return anac
        # if re.match('^ *\[', line):
        #     if re.match('^ *\[ atoms \]', line):
        #         print line,
        #     else:

# Backbone goes like N, CT, C, N, CT, C...

def is_backbone(dac):
    return dac[1:3] in [['N', 'CT'], ['CT', 'C'], ['C', 'N'], ['CT', 'N'], ['C', 'CT'], ['N', 'C']]

print AACodes

allbac = []
allaac = []
alldac = []
aadac = {}
for iAA, AA in enumerate(AACodes):
    print "Detecting parameters for", AA
    anac = build_atomname_to_atomclass(AA)
    anac.update(build_atomname_to_atomclass('ACE'))
    anac.update(build_atomname_to_atomclass('NME'))
    M = Molecule(os.path.join(AADirs[iAA], 'all.gro'))
    aadac[AA] = []
    for i in M.atomname:
        if i not in anac.keys():
            print '%s not in list of atom names' % i
    for d in M.find_dihedrals():
        if all([M.atomname[i] in anac for i in d]):
            dac = [anac[M.atomname[i]] for i in d]
            if not is_backbone(dac):
                for mult in range(1, 7):
                    dstr = "#define torsion_%s_%s_%s_%s_%s_mult%i   0.0   0.0   %i" % (AA, dac[0], dac[1], dac[2], dac[3], mult, mult)
                    if dstr not in aadac[AA]:
                        aadac[AA].append(dstr)
                        print dstr
            if is_backbone(dac):
                if args.bk and dac not in alldac:
                    alldac.append(dac[:])
                    alldac.append(dac[::-1])
            elif args.sd and dac not in alldac:
                alldac.append(dac[:])
                alldac.append(dac[::-1])
        else:
            print dac, "is not parameterized"
    for a in M.find_angles():
        if all([M.atomname[i] in anac for i in a]):
            aac = [anac[M.atomname[i]] for i in a]
            if args.an and aac not in allaac:
                allaac.append(aac[:])
                allaac.append(aac[::-1])
    for b in M.bonds:
        if all([M.atomname[i] in anac for i in b]):
            bac = [anac[M.atomname[i]] for i in b]
            if args.bd and bac not in allbac:
                allbac.append(bac[:])
                allbac.append(bac[::-1])
        
mode = 'N'
foutnm = '%s%s%s' % (os.path.splitext(sys.argv[1])[0], label, os.path.splitext(sys.argv[1])[1])
print "Output file is", foutnm
fout = open(foutnm,'w')
nprm = 0

for line in open(sys.argv[1]).readlines():
    line = re.sub(';.*$', '', line.replace('\n',''))
    match = re.match('^ *\[ *([^ ]*) *\]', line)
    if match:
        print >> fout
        content = match.group(1)
        if content == 'dihedraltypes':
            mode = 'D'
        elif content == 'angletypes':
            mode = 'A'
        elif content == 'bondtypes':
            mode = 'B'
        else:
            mode = 'N'
    elif mode == 'D':
        if len(re.sub(';.*$', '', line).strip()) > 0:
            # We are now at a line of dihedral data
            s = line.split()
            if any([all([(i == j or j == 'X') for i, j in zip(dac, s[:4])]) for dac in alldac]):
                if args.eq:
                    line += ' ; PRM 5 6'
                    nprm += 2
                else:
                    line += ' ; PARM 6'
                    nprm += 1
    elif mode == 'A':
        if len(re.sub(';.*$', '', line).strip()) > 0:
            # We are now at a line of angle data
            s = line.split()
            if any([all([i == j for i, j in zip(aac, s[:3])]) for aac in allaac]):
                if args.eq:
                    line += ' ; PRM 4 5'
                    nprm += 2
                else:
                    line += ' ; PARM 5'
                    nprm += 1
    elif mode == 'B':
        if len(re.sub(';.*$', '', line).strip()) > 0:
            # We are now at a line of bond data
            s = line.split()
            if any([all([i == j for i, j in zip(bac, s[:2])]) for bac in allbac]):
                if args.eq:
                    line += ' ; PRM 3 4'
                    nprm += 2
                else:
                    line += ' ; PARM 4'
                    nprm += 1
    if len(line) > 0:
        print >> fout, line

print nprm, "total fitting parameters"

fout.close()
            # if any([all([i == j or j == 'X'] for i, j in zip(dac, s[:4])) for dac in alldac]):
            #     print s[:4]
            
    # print sorted(list(set(anac.keys())))
    # print sorted(list(set(M.atomname)))

sys.exit()
