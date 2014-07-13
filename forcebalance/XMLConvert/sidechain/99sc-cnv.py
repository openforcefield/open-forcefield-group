#!/usr/bin/env python

# The purpose of this script is to port LP's optimized intramolecular parameters for AMBER99SB
# into the OpenMM XML format.

# The strategy is to:
# 1) Read the OpenMM XML file
# 2) Read the GROMACS .itp file to figure out the interactions defined using atom classes
# 3) Write the new OpenMM XML file

from simtk.openmm.app import *
import itertools
from copy import deepcopy
import lxml.etree as ET
import numpy as np
import os, sys, re
import networkx as nx
from collections import defaultdict, OrderedDict

# This version incorporates residue-specific side chain torsions.
# What a painful script to write!!

# Parse the original AMBER99SB XML file.
A99SB = ET.parse('/home/leeping/src/OpenMM/wrappers/python/simtk/openmm/app/data/amber99sb.xml')

# Use GromacsTopFile to read the optimized. ITP file.
ITP = GromacsTopFile('/home/leeping/projects/VSP27-Protein/Dihedrals/Optimize/AMBER99SC/optimize.r1.mmc3/a99sc-v2.itp')

# Read the residue definitions with specific dihedral interactions.
RTP = '/home/leeping/opt/gromacs/share/gromacs/top/a99sc.ff/aminoacids.rtp'

# Parsed [ bondtypes ], [ angletypes ] and [ dihedraltypes ] sections
BT = ITP._bondTypes
AT = ITP._angleTypes
DT = ITP._dihedralTypes

root = A99SB.getroot()

# Amino Acid Dihedral Quartets (by New Atom Classes) to Dihedral Parameters
AA_DSC = OrderedDict()
# Amino Acid Atom Names to OpenMM Atom Types
AA_OAt = OrderedDict()
# Gromacs Atom Names to Atom Class
GAnAt = {}
# OpenMM Atom Types to Atom Class
OAtAc = {}

# Manually constructed atom name replacements.
Atom_Rename = {'ILE':{'CD':'CD1'}}
for i,j in Atom_Rename.items():
    Atom_Rename['C'+i] = j
    Atom_Rename['N'+i] = j

# Mapping of amino acid atom names to new atom classes.  Mostly new
# atom classes for beta carbons but a few new gamma carbons are
# defined.  They are named using the number "6" (for carbon) and the
# one-letter amino acid code.  A few exceptions in the case of alternate
# protonation states.
NewAC = {"SER":{"CB":"6S"},
         "THR":{"CB":"6T", "CG2":"6t"},
         "LEU":{"CB":"6L"},
         "VAL":{"CB":"6V"},
         "ILE":{"CB":"6I", "CG2":"6i"},
         "ASN":{"CB":"6N"},
         "GLN":{"CB":"6Q", "CG":"6q"},
         "ARG":{"CB":"6R"},
         "HID":{"CB":"6H"},
         "HIE":{"CB":"6h"},
         "HIP":{"CB":"6+"},
         "TRP":{"CB":"6W"},
         "TYR":{"CB":"6Y"},
         "PHE":{"CB":"6F"},
         "GLU":{"CB":"6E", "CG":"6e"},
         "ASP":{"CB":"6D"},
         "LYS":{"CB":"6K"},
         "LYN":{"CB":"6k"},
         "PRO":{"CB":"6P"},
         "CYS":{"CB":"6C"},
         "CYM":{"CB":"6c"},
         "MET":{"CB":"6M"},
         "ASH":{"CB":"6d"},
         "GLH":{"CB":"6J", "CG":"6j"}}

for i in NewAC.keys():
    NewAC['C'+i] = NewAC[i]
    NewAC['N'+i] = NewAC[i]

# Obtain a canonicalized ordering of dihedral atom classes.
def get_ijkl(aci, acj, ack, acl):
    if ack > acj:
        acijkl = tuple([aci, acj, ack, acl])
    elif ack < acj:
        acijkl = tuple([acl, ack, acj, aci])
    else:
        if acl >= aci:
            acijkl = tuple([aci, acj, ack, acl])
        else:
            acijkl = tuple([acl, ack, acj, aci])
    return acijkl

# Read amino acid definitions;
# This is so we can figure out which atom class quartets receive
# the sidechain-specific dihedral parameters.
def ParseRTP(rtp):
    for line in open(rtp).readlines():
        line = line.split(';')[0].strip()
        s = line.split()
        if len(s) == 0: continue
        if re.match('\[ .* \]$', line):
            section = s[1]
            if section not in ['bondedtypes', 'atoms', 'bonds', 'dihedrals', 'impropers']:
                AA = section
                GAnAt[AA] = dict()
        elif section == 'atoms':
            # The GROMACS atom types are the same as OpenMM atom classes
            # and they should serve as the default atom class when we 
            # haven't defined a new one
            GAnAt[AA][s[0]] = s[1]
        elif section == 'dihedrals':
            dihan = tuple(s[:4])
            # Obtain the quartet of new atom classes corresponding to this particular dihedral interaction
            aci, acj, ack, acl = [NewAC.get(AA, {}).get(an,GAnAt[AA][an]) for an in dihan]
            acijkl = get_ijkl(aci, acj, ack, acl)
            # Insert the dihedral parameters into AA_DSC (keyed by the quartet of new atom classes)
            if acijkl not in AA_DSC:
                AA_DSC[acijkl] = []
            if ITP._defines[s[4]].split() not in AA_DSC[acijkl]:
                AA_DSC[acijkl].append(ITP._defines[s[4]].split())

ParseRTP(RTP)

def almostequal(i, j, tol):
    return np.abs(i-j) < tol

def get_periodn(elem, period):
    nprd = 0
    for i, j in elem.items():
        if 'periodicity' in i:
            nprd += 1
            if j == period:
                return i[-1]
    return "%i" % (nprd + 1)

NewAtAc = deepcopy(OAtAc)

# Amino Acid Graphs
AA_Gphs = OrderedDict()
# A dictionary of atom classes to copy FROM -> [TO]
FrcCopy = defaultdict(list)
# A list of bond types
NeedBT = set()
NeedAT = set()
NeedDT = set()
NeedIT = set()
HaveBT = set()
HaveAT = set()
HaveDT = set()
HaveIT = set()
for force in root:
    # Check for atom classes that are missing from the ITP file
    if force.tag == 'AtomTypes':
        for elem in force:
            OAtAc[elem.attrib['name']] = elem.attrib['class']
            if elem.attrib['class'] not in ITP._atomTypes:
                print "Atom Type", elem.attrib['name'], "Class", elem.attrib['class'], "not present in ITP file"

    # Residue processing
    if force.tag == 'Residues':
        for elem in force:
            res = elem.attrib['name']
            # Initialize NetworkX graph
            AA_Gphs[res] = nx.Graph()
            G = AA_Gphs[res]
            # List of (new) atom classes for each atom in the residue
            resAc = []
            # Number of nodes in the graph
            nn = 0
            for subelem in elem:
                # Atom tag: Create a NetworkX node
                if subelem.tag == 'Atom':
                    G.add_node(nn)
                    G.node[nn]['name'] = subelem.attrib['name']
                    G.node[nn]['type'] = subelem.attrib['type']
                    G.node[nn]['class'] = OAtAc[subelem.attrib['type']]
                    nn += 1
                    # FrcCopy is a dictionary denoting which interactions will get copied
                    # due to atom class duplication (e.g. (CT, CT) -> (CT, 6S))
                    if res in NewAC and subelem.attrib['name'] in NewAC[res]:
                        FrcCopy[OAtAc[subelem.attrib['type']]].append(NewAC[res][subelem.attrib['name']])
                    # Create OpenMM atom type dictionary.
                    if res not in AA_OAt:
                        AA_OAt[res] = OrderedDict()
                    AA_OAt[res][subelem.attrib['name']] = subelem.attrib['type']
                    # resAc will be used to create pairs, triplets and quartets of atom types in this residue
                    # for constructing the NeedBT, NeedAT, NeedDT and NeedIT dictionaries
                    resAc.append(NewAC.get(res, {}).get(subelem.attrib['name'], OAtAc[subelem.attrib['type']]))
                    # Record of new atom type -> atom class mapping; used to rewrite AtomTypes section.
                    NewAtAc[subelem.attrib['type']] = NewAC.get(res, {}).get(subelem.attrib['name'], OAtAc[subelem.attrib['type']])
                if subelem.tag == 'Bond':
                    # Add edges in the graph.
                    G.add_edge(int(subelem.attrib['from']), int(subelem.attrib['to']))
                if subelem.tag == 'ExternalBond':
                    # Add edges in the graph for bonds to the next residue.
                    ifrom = int(subelem.attrib['from'])
                    if G.node[ifrom]['class'] == 'C':
                        G.add_node(nn)
                        G.node[nn]['name'] = 'N'
                        G.node[nn]['type'] = '706'
                        G.node[nn]['class'] = 'N'
                        G.add_edge(ifrom, nn)
                        resAc.append('N')
                        nn += 1
                    elif G.node[ifrom]['class'] == 'N':
                        G.add_node(nn)
                        G.node[nn]['name'] = 'C'
                        G.node[nn]['type'] = '711'
                        G.node[nn]['class'] = 'C'
                        G.add_edge(ifrom, nn)
                        resAc.append('C')
                        nn += 1
                    # For now, don't treat nucleic acids and disulfide bonds.
                    elif G.node[ifrom]['name'] in ['SG', 'P', "O3'"]:
                        pass
                    else:
                        print G.node[ifrom]
                        raise RuntimeError('Spoo!')

            # Build the NeedBT, NeedAT, NeedDT and NeedIT sets.
            # These are the atom class combinations for the bond / angle / dihedral 
            # interactions that actually occur in the residues.
            # If we didn't have these as a filter, then one atom class that gets expanded to 30 (i.e. CT -> 6T, 6L, 6V, etc...)
            # would suddenly blow up to a humongous number of dihedral types.
            for edge in G.edges():
                # Build a list of bond types in this residue.
                NeedBT.add(tuple(sorted([resAc[edge[0]], resAc[edge[1]]])))
            for a2 in list(G.nodes()):
                # List of angle types in this residue.
                # Find all bonded neighbors to this atom
                friends = sorted(list(nx.all_neighbors(G, a2)))
                if len(friends) < 2: continue
                # Double loop over bonded neighbors
                for i, a1 in enumerate(friends):
                    for a3 in friends[i+1:]:
                        c1, c2, c3 = (resAc[k] for k in [a1, a2, a3])
                        if c3 > c1:
                            NeedAT.add(tuple((c1, c2, c3)))
                        else:
                            NeedAT.add(tuple((c3, c2, c1)))
            for edge in G.edges():
                # List of proper dihedral types.
                a2 = edge[0]
                a3 = edge[1]
                for a1 in sorted(list(nx.all_neighbors(G, a2))):
                    if a1 != a3:
                        for a4 in sorted(list(nx.all_neighbors(G, a3))):
                            if a4 != a2:
                                c1, c2, c3, c4 = (resAc[i] for i in [a1, a2, a3, a4])
                                ijkl = get_ijkl(c1, c2, c3, c4)
                                NeedDT.add((ijkl[0], ijkl[1], ijkl[2], ijkl[3]))
                                NeedDT.add((ijkl[3], ijkl[2], ijkl[1], ijkl[0]))
                                NeedDT.add(("", ijkl[1], ijkl[2], ""))
                                NeedDT.add(("", ijkl[2], ijkl[1], ""))
            for a1 in list(G.nodes()):
                # List of improper dihedral types.
                # Find all bonded neighbors to this atom
                friends = sorted(list(nx.all_neighbors(G, a1)))
                if len(friends) < 3: continue
                for a2, a3, a4 in list(itertools.permutations(friends, 3)):
                    c1, c2, c3, c4 = (resAc[i] for i in [a1, a2, a3, a4])
                    # Include wildcards.
                    NeedIT.add((c1, c2, c3, c4))
                    NeedIT.add((c1, "", c3, c4))
                    NeedIT.add((c1, "", "", c4))
        # FrcCopy is responsible for taking existing atom classes and
        # propagating the interaction types out to the copied atom
        # classes.
        for i in FrcCopy.keys():
            FrcCopy[i].append(i)

    # Harmonic bond parameters
    if force.tag == 'HarmonicBondForce':
        # List of new interaction types if needed
        newfrc = []
        for elem in force:
            att = elem.attrib
            BC = (att['class1'], att['class2'])
            BCr = (att['class2'], att['class1'])
            # Look up parameters from the parsed GROMACS ITP file
            if BC in BT.keys():
                prm = BT[BC]
            elif BCr in BT.keys():
                prm = BT[BCr]
            else:
                print BC, "has no parameters from the ITP file"
                prm = None
            # Set parameters if they differ from the OpenMM values
            if prm != None:
                if not almostequal(float(elem.attrib['length']), float(prm[3]), 1e-8):
                    elem.attrib['length'] = '%.8f' % float(prm[3])
                if not almostequal(float(elem.attrib['k']), float(prm[4]), 1e-8):
                    elem.attrib['k'] = '%.8f' % float(prm[4])
            acij = tuple(sorted(BC))
            # Add interaction type to HaveBT to avoid double counting
            HaveBT.add(acij)
            # Copy harmonic bond parameters to "copied" interaction types.
            for aci, acj in itertools.product(*[FrcCopy.get(BC[0], [BC[0]]), FrcCopy.get(BC[1], [BC[1]])]):
                acij = tuple(sorted([aci, acj]))
                if acij in NeedBT and acij not in HaveBT:
                    elem1 = deepcopy(elem)
                    elem1.attrib['class1'] = acij[0]
                    elem1.attrib['class2'] = acij[1]
                    newfrc.append(elem1)
                    HaveBT.add(acij)
        for elem in newfrc:
            force.append(elem)

    # Harmonic angle parameters.  Same as for harmonic bonds.
    if force.tag == 'HarmonicAngleForce':
        newfrc = []
        for elem in force:
            att = elem.attrib
            AC = (att['class1'], att['class2'], att['class3'])
            ACr = (att['class3'], att['class2'], att['class1'])
            if AC in AT.keys():
                prm = AT[AC]
            elif ACr in AT.keys():
                prm = AT[ACr]
            else:
                print AC, "has no parameters from the ITP file"
                prm = None
            if prm != None:
                gang = float(prm[4]) * np.pi / 180
                oang = float(elem.attrib['angle'])
                if not almostequal(gang, oang, 1e-8):
                    elem.attrib['angle'] = '%.8f' % gang
                if not almostequal(float(prm[5]), float(elem.attrib['k']), 1e-8):
                    elem.attrib['k'] = '%.8f' % float(prm[5])
                    
            if AC[2] >= AC[0]:
                acijk = tuple(AC)
            else:
                acijk = tuple(ACr)
            HaveAT.add(acijk)

            # Duplicate harmonic angle parameters for new atom classes.
            for aci, acj, ack in itertools.product(*[FrcCopy.get(AC[0], [AC[0]]), FrcCopy.get(AC[1], [AC[1]]), FrcCopy.get(AC[2], [AC[2]])]):
                if ack >= aci:
                    acijk = tuple([aci, acj, ack])
                else:
                    acijk = tuple([ack, acj, aci])
                if acijk in NeedAT and acijk not in HaveAT:
                    elem1 = deepcopy(elem)
                    elem1.attrib['class1'] = acijk[0]
                    elem1.attrib['class2'] = acijk[1]
                    elem1.attrib['class3'] = acijk[2]
                    newfrc.append(elem1)
                    HaveAT.add(acijk)
        for elem in newfrc:
            force.append(elem)

    # Periodic torsion parameters.  These are a bit of a pain.
    if force.tag == 'PeriodicTorsionForce':
        newfrc = []
        for elem in force:
            att = elem.attrib
            
            # Since I got confused about the orderings, I add both orderings to HaveDT
            if elem.tag == 'Proper':
                DC = (att['class1'], att['class2'], att['class3'], att['class4'], '9') 
                DCr = (att['class4'], att['class3'], att['class2'], att['class1'], '9')
                HaveDT.add((att['class1'], att['class2'], att['class3'], att['class4']))
                HaveDT.add((att['class4'], att['class3'], att['class2'], att['class1']))
            elif elem.tag == 'Improper':
                DC = (att['class2'], att['class3'], att['class1'], att['class4'], '4') 
                DCr = (att['class2'], att['class3'], att['class4'], att['class1'], '4')
                HaveIT.add((att['class1'], att['class2'], att['class3'], att['class4']))

            # Look up parameters from the Gromacs ITP file
            DC = tuple('X' if i == '' else i for i in DC)
            DCr = tuple('X' if i == '' else i for i in DCr)
            if DC in DT.keys():
                prms = DT[DC]
            elif DCr in DT.keys():
                prms = DT[DCr]
            else:
                print DC, "has no parameters from the ITP file"
                prms = None

            # Edit parameters in-place for the existing interaction type
            if prms != None:
                for prm in prms:
                    prd = prm[7]
                    prdn = get_periodn(elem, prd)
                    if ('periodicity' + prdn) not in elem.attrib:
                        elem.attrib['periodicity' + prdn] = prd
                        elem.attrib['phase' + prdn] = '%.8f' % (float(prm[5]) * np.pi / 180)
                        elem.attrib['k' + prdn] = '%.8f' % (float(prm[6]))
                    else:
                        if not almostequal(float(prm[5]) * np.pi / 180, float(elem.attrib['phase' + prdn]), 1e-8):
                            elem.attrib['phase' + prdn] = '%.8f' % (float(prm[5]) * np.pi / 180)
                        if not almostequal(float(prm[6]), float(elem.attrib['k' + prdn]), 1e-8):
                            elem.attrib['k' + prdn] = '%.8f' % (float(prm[6]))

            # Propagate interaction type to "copied" types
            if elem.tag == 'Improper':
                for dci, dcj, dck, dcl in itertools.product(*[FrcCopy.get(att['class1'], [att['class1']]),
                                                              FrcCopy.get(att['class2'], [att['class2']]),
                                                              FrcCopy.get(att['class3'], [att['class3']]),
                                                              FrcCopy.get(att['class4'], [att['class4']])]):
                    ijkl = (dci, dcj, dck, dcl) 
                    if ijkl in NeedIT and ijkl not in HaveIT:
                        elem1 = deepcopy(elem)
                        elem1.attrib['class1'] = dci
                        elem1.attrib['class2'] = dcj
                        elem1.attrib['class3'] = dck
                        elem1.attrib['class4'] = dcl
                        newfrc.append(elem1)
                        HaveIT.add(ijkl)

            elif elem.tag == 'Proper':
                for dci, dcj, dck, dcl in itertools.product(*[FrcCopy.get(att['class1'], [att['class1']]),
                                                              FrcCopy.get(att['class2'], [att['class2']]),
                                                              FrcCopy.get(att['class3'], [att['class3']]),
                                                              FrcCopy.get(att['class4'], [att['class4']])]):
                    ijkl = (dci, dcj, dck, dcl) 
                    lkji = ijkl[::-1]
                    if (ijkl in NeedDT or lkji in NeedDT) and ijkl not in HaveDT:
                        elem1 = deepcopy(elem)
                        elem1.attrib['class1'] = dci
                        elem1.attrib['class2'] = dcj
                        elem1.attrib['class3'] = dck
                        elem1.attrib['class4'] = dcl
                        newfrc.append(elem1)
                        HaveDT.add(ijkl)
                        HaveDT.add(lkji)

        for elem in newfrc:
            force.append(elem)
            
        # Finally get side chain parameters from the ITP file.
        newfrc = []
        # Needed interaction types are pulled from the atom class quartets constructed in ParseRTP
        for ijkl in NeedDT:
            # Look up the interaction type in the ones from the ITP file
            if ijkl in AA_DSC:
                print ijkl, "sidechain parameters:",
                prm = AA_DSC[ijkl]
                replace = False
                # Loop through the existing interaction types in the XML file.  If already exists,
                # then replace the parameters.
                for elem in force:
                    if (ijkl == (elem.attrib['class1'], elem.attrib['class2'], elem.attrib['class3'], elem.attrib['class4']) or 
                        ijkl == (elem.attrib['class4'], elem.attrib['class3'], elem.attrib['class2'], elem.attrib['class1'])):
                        print "replacing existing"
                        replace = True
                        for iixn, ixn in enumerate(prm):
                            elem.attrib["periodicity%i" % (iixn + 1)] = ixn[2]
                            elem.attrib["phase%i" % (iixn + 1)] = '%.10f' % (float(ixn[0]) * np.pi / 180)
                            elem.attrib["k%i" % (iixn + 1)] = ixn[1]
                        break
                # If the interaction type doesn't exist in the XML file, create it anew.
                if not replace:
                    print "creating new"
                    elem1 = ET.Element("Proper")
                    elem1.attrib['class1'] = ijkl[0]
                    elem1.attrib['class2'] = ijkl[1]
                    elem1.attrib['class3'] = ijkl[2]
                    elem1.attrib['class4'] = ijkl[3]
                    for iixn, ixn in enumerate(prm):
                        elem1.attrib["periodicity%i" % (iixn + 1)] = ixn[2]
                        elem1.attrib["phase%i" % (iixn + 1)] = '%.10f' % (float(ixn[0]) * np.pi / 180)
                        elem1.attrib["k%i" % (iixn + 1)] = ixn[1]
                    elem1.tail = '\n  '
                    newfrc.append(elem1)

        for elem in newfrc:
            force.append(elem)

# Finally replace the atom classes in the AtomTypes section.
for force in root:
    if force.tag == 'AtomTypes':
        for elem in force:
            if OAtAc[elem.attrib['name']] != NewAtAc[elem.attrib['name']]:
                elem.attrib['class'] = NewAtAc[elem.attrib['name']]

# Write the output. (Whew!)
with open('new.xml', 'wb') as f:
    A99SB.write(f)
