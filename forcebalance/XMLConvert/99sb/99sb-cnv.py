#!/usr/bin/env python

# The purpose of this script is to port LP's optimized intramolecular parameters for AMBER99SB
# into the OpenMM XML format.

# The strategy is to:
# 1) Read the OpenMM XML file
# 2) Read the GROMACS .itp file to figure out the interactions defined using atom classes
# 3) Write the new OpenMM XML file

# This script was easy to write as it simply replaces numbers in the XML file.

from simtk.openmm.app import *
import lxml.etree as ET
import numpy as np
import os, sys

A99SB = ET.parse('/home/leeping/src/OpenMM/wrappers/python/simtk/openmm/app/data/amber99sb.xml')
ITP = GromacsTopFile(sys.argv[1])
BT = ITP._bondTypes
AT = ITP._angleTypes
DT = ITP._dihedralTypes

root = A99SB.getroot()

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

for force in root:
    # Check for atom classes that are missing from the ITP file
    if force.tag == 'AtomTypes':
        for elem in force:
            if elem.attrib['class'] not in ITP._atomTypes:
                print "Atom Type", elem.attrib['name'], "Class", elem.attrib['class'], "not present in ITP file"

    # Copy over harmonic bond parameters
    if force.tag == 'HarmonicBondForce':
        for elem in force:
            att = elem.attrib
            BC = (att['class1'], att['class2'])
            BCr = (att['class2'], att['class1'])
            if BC in BT.keys():
                prm = BT[BC]
            elif BCr in BT.keys():
                prm = BT[BCr]
            else:
                print BC, "has no parameters from the ITP file"
                prm = None
            if prm != None:
                if not almostequal(float(elem.attrib['length']), float(prm[3]), 1e-8):
                    elem.attrib['length'] = '%.8f' % float(prm[3])
                if not almostequal(float(elem.attrib['k']), float(prm[4]), 1e-8):
                    elem.attrib['k'] = '%.8f' % float(prm[4])

    if force.tag == 'HarmonicAngleForce':
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

    if force.tag == 'PeriodicTorsionForce':
        for elem in force:
            att = elem.attrib
            if elem.tag == 'Proper':
                DC = (att['class1'], att['class2'], att['class3'], att['class4'], '9') 
                DCr = (att['class4'], att['class3'], att['class2'], att['class1'], '9')
            elif elem.tag == 'Improper':
                DC = (att['class2'], att['class3'], att['class1'], att['class4'], '4') 
                DCr = (att['class4'], att['class1'], att['class2'], att['class3'], '4')
            DC = tuple('X' if i == '' else i for i in DC)
            DCr = tuple('X' if i == '' else i for i in DCr)
            if DC in DT.keys():
                prms = DT[DC]
            elif DCr in DT.keys():
                prms = DT[DCr]
            else:
                print DC, "has no parameters from the ITP file"
                prms = None
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
                
A99SB.write('new.xml')
