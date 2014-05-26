#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
split-molecule-set.py

Split molecule set into training and test sets.

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

from optparse import OptionParser # For parsing of command line arguments

import os
import math
import numpy
import re
import string
import sys
import time

import openeye.oechem
import openeye.oequacpac
import openeye.oeiupac
import openeye.oeomega

#=============================================================================================
# MAIN
#=============================================================================================

if __name__=="__main__":

    # Create command-line argument options.
    usage_string = """\
    usage: %prog --molecules molfile --test testfile --train trainfile
    
    example: %prog --molecules datasets/neutrals.sdf --test datasets/test.sdf --train datasets/train.sdf
    
    """
    version_string = "%prog %__version__"
    parser = OptionParser(usage=usage_string, version=version_string)

    parser.add_option("-t", "--table", metavar='TABLE',
                      action="store", type="string", dest='table_filename', default='results_table_origcarbon.tex',
                      help="LaTeX table from SI.")
    parser.add_option("-o", "--output", metavar='OUTPUT',
                      action="store", type="string", dest='output_filename', default='mobley.sdf',
                      help="SDF file to write.")
    
    # Parse command-line arguments.
    (options,args) = parser.parse_args()

    # Read LaTeX table.
    infile = open(options.table_filename, 'r')
    lines = infile.readlines()
    infile.close()

    # Skip header and footer.
    lines = lines[3:]
    lines = lines[:-2]

    # Process table.
    molecules = list()
    unparsed_molecule_names = list()
    for line in lines:
        # Trim training \\\hline
        line = string.replace(line, r'\\\hline', '')
        
        # Split into elements.
        # EXAMPLE:
        # 1112\_tetrachloroethane & $-1.43\pm 0.01 $ & $1.54\pm0.02$  &  $0.11\pm0.02 $ & -1.28\\\hline
        elements = line.split('&')

        # Parse elements.
        molecule_name = elements[0]
        dg_exp = float(elements[4]) # in kcal/mol

        # Process molecule name.
        # Tackle acids.
        molecule_name = string.replace(molecule_name, r'\_acid', r' acid')
        # Change protected understrikes into dashes.
        molecule_name = string.replace(molecule_name, r'\_', r'-')
        # Eliminate double dashes.
        molecule_name = string.replace(molecule_name, r'--', r'-')
        # Change repeated numbers to have commas.
        count = 1
        while (count > 0):
            (molecule_name, count) = re.subn(r'(\d)(\d)', r'\1,\2', molecule_name)
        print "%64s %8.3f kcal/mol" % (molecule_name, dg_exp)

        try:
            # Create molecule.
            molecule = openeye.oechem.OEMol()
            status = openeye.oeiupac.OEParseIUPACName(molecule, molecule_name) 

            # Assign aromaticity.
            #OEAssignAromaticFlags(molecule, OEAroModelOpenEye)   
            
            # Add hydrogens.
            openeye.oechem.OEAddExplicitHydrogens(molecule)

            # Set title.
            molecule.SetTitle(molecule_name)
            
            # Set free energy.
            openeye.oechem.OESetSDData(molecule, 'dG(exp)', '%f' % dg_exp);

            # Add molecule to list.
            molecules.append(molecule)
        except:
            unparsed_molecule_names.append(molecule_name)

    print "%d unparsed molecule names:" % len(unparsed_molecule_names)
    for name in unparsed_molecule_names:
        print name

    # Build a conformation for all molecules with Omega.
    print "Building conformations for all molecules..."
    start_time = time.time()    
    omega = openeye.oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetStrictStereo(False)
    omega.SetFromCT(True)
    for molecule in molecules:
        #omega.SetFixMol(molecule)
        omega(molecule)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print "%.3f s elapsed" % elapsed_time

    # Write SDF file.
    print "Writing test set..."
    omolstream = openeye.oechem.oemolostream(options.output_filename)    
    for molecule in molecules:
        # Write molecule as mol2, changing molecule through normalization.    
        openeye.oechem.OEWriteMolecule(omolstream, molecule)
    omolstream.close()
