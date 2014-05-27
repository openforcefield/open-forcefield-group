#!/usr/bin/env python

import os, sys
from forcebalance.molecule import Molecule
from collections import defaultdict

# Atom classes belonging to the backbone or termini.
exc_mid = [1, 2, 3, 7, 40]

# Atom classes belonging to methyl hydrogens.
exc_end = []#6, 9]

# Atom class names from TINKER.
acnames = ["", "Backbone Amide Nitrogen", "Glycine Alpha Carbon", "Backbone Carbonyl Carbon", 
           "Amide or Guanidinium Hydrogen", "Amide Carbonyl Oxygen", "Methine Hydrogen", 
           "Methine Carbon", "Methyl or Methylene Carbon", "Methyl or Methylene Hydrogen", 
           "Hydroxyl Oxygen", "Hydroxyl Hydrogen", "Sulfide or Disulfide Sulfur", 
           "Sulfhydryl Hydrogen", "Thiolate Sulfur", "Proline Backbone Nitrogen", 
           "Proline Ring Methylene Carbon", "Phenyl Carbon", "Phenyl Hydrogen", 
           "Phenolic Oxygen", "Phenolic Hydrogen", "Phenoxide Oxygen", "Indole Carbon", 
           "Indole CH Hydrogen", "Imidazole or Indole NH Nitrogen", 
           "Imidazole or Indole NH Hydrogen", "Imidazole C=C Carbon", "Imidazole CH Hydrogen", 
           "Imidazole N=C-N Carbon", "Imidazole C=N Nitrogen", "Carboxylate Carbon", 
           "Carboxylate Oxygen", "Carboxylic Acid Carbonyl Oxygen", 
           "Carboxylic Acid Carbonyl Oxygen", "Carboxylic Acid Hydroxyl Oxygen", 
           "Carboxylic Acid Hydrogen", "Lysine/Ornithine Gamma Carbon", "Ammonium Nitrogen", 
           "Ammonium Hydrogen", "Guanidinium Carbon", "Acetyl or NMe Methyl Carbon", 
           "N-Terminal Ammonium Nitrogen", "N-Terminal Ammonium Hydrogen"]

# Mapping from torsional atom classes -> line number in file
tdict = {}
# Mapping from torsional atom classes -> atoms in interaction
tatoms = defaultdict(list)
# Mapping from atom types -> atom class
acdict = {}
# Force field loaded into text file
fftxt = open('amoebapro13.prm').readlines()
for ln, line in enumerate(fftxt):
    s = line.split()
    # Build the torsional dihedraltype -> line number mapping
    if len(s) > 0 and s[0] == 'torsion':
        tdict[tuple(int(i) for i in s[1:5])] = ln
    # Build the atom type -> atom class mapping
    if len(s) > 0 and s[0] == 'atom':
        acdict[int(s[1])] = int(s[2])

# Canonicalize dihedral atom number (or atom class) ordering
def cand(dih):
    # Make sure two center atoms aren't out of order
    if (dih[1] > dih[2]):
        return (dih[3], dih[2], dih[1], dih[0])
    # Two center atoms are equal but 1st and 4th are out of order
    elif dih[1] == dih[2] and dih[0] > dih[3]:
        return (dih[3], dih[2], dih[1], dih[0])
    else:
        return dih

hilite = set()
for dnm in os.listdir("../targets"):
    if 'SER' not in dnm: continue
    if 'chi1chi2' in dnm: continue
    print dnm
    M = Molecule(os.path.join("../targets", dnm, 'all.arc'))
    atomtypes = [int(suf.split()[0]) for suf in M.tinkersuf]
    # for ia, a in enumerate(atomtypes):
    #     print ia, acdict[a], acnames[acdict[a]]
    set_dclass = set()
    for dih in M.find_dihedrals():
        dclass = cand(tuple(acdict[atomtypes[i]] for i in dih))
        # Do not highlight backbone torsions.
        if dclass[1] in exc_mid and dclass[2] in exc_mid: continue
        # Do not highlight methyl hydrogens.
        if dclass[0] in exc_end or dclass[3] in exc_end: continue
        hilite.add(tdict[dclass])
        tatoms[tdict[dclass]].append(dih)
        set_dclass.add(dclass)
    for d in sorted(list(set_dclass)):
        print d
    print
    
dflds = [5, 8, 11]
nprm = 0
fout = open('amoebapro13_labeled.prm','w')
for ln, line in enumerate(fftxt):
    nprmi = 0
    pfld = []
    print >> fout, line.replace('\n',''),
    if ln in hilite:
        for i in dflds:
            #if float(line.split()[i]) != 0.0:
            if 1: # Highlight ALL dihedrals
                nprmi += 1
                nprm += 1
                pfld.append(i)
        if len(pfld) > 0:
            print >> fout, "# PRM", ' '.join(['%i' % i for i in pfld]),
    print >> fout

print "%i total parameters" % nprm
