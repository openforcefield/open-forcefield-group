#!/usr/bin/python

"""
A module for the optimization of parameters

CL: call from ~/SVN/proteinparam/test
module load anaconda
python ../scripts/get_properties.py
"""

#=========================================================================
# COPYRIGHT NOTICE
#
# Written by Brittany S. Zimmerman <brittanys.zimmerman@gmail.com>
#
# Portions of this code modified from code written by Michael R. Shirts, P. Klimovich, D. Mobley, and John Chodera (alchemical-gromacs.py and pymbar.py, both available via svn checkout https://simtk.org/svn/pymbar)
# 
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the Licence, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; including but not limited to the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#=========================================================================


#=========================================================================
# REFERENCES
# [1] Wang, Junmei; Hou, Tingjun. Application of Molecular Dynamics Simulations in Molecular Property Prediction. 1. Density and Heat of Vaporization. JCTC. 2011; 7: 2151-2165
# [2] Shirts, Michael R.; Pande, Vinjay S. Solvation free energies of amino acid side chain analogs for common molecular mechanics water models. J. Chem. Phys. 122, 134508 (2005)
#=========================================================================


#=========================================================================
# TO DO
#
# HIGH PRIORITY
# 2. Doesn't work with optimizers that don't use derivatives.  Doesn't rerun in main_single_point is the issue I think.
# 3. Trim down Constants and Parameters section - as much of this as possible should be in the input option parsing section or automatically generated from the files so that users dont have to change the acutal script text
#
# LOW PRIORITY
# 1. Add a check for normality of Volume distribution for 'Generate Delta P' section.  If distribution is not normal then generate a warning to suggest running longer initial simulations.
# 2. Make the nsamples override also override the rerun length (shorten trr files and edit mdp file).
#========================================================================


#========================================================================
# VERSION CONTROL INFORMATION
#========================================================================
__version__ = '1.0 beta'
__authors__ = 'Brittany S. Zimmerman'
__license__ = 'GPL 2.0'


#========================================================================
# IMPORTS
#========================================================================
from optparse import OptionParser # for parsing command-line options
import numpy as np
import scipy
from commands import getoutput # in Python 3.x, getoutput is moved to the subprocess module
import re
import pdb
import os
import subprocess
import time
import datetime
import sys
import string
import stat
#import matplotlib.pyplot as plt
#from scipy.optimize import minimize
# IMPORT MODULES
sys.path.append('/home/bsz9ur/pymbar/trunk/pymbar')
import pymbar
import timeseries # For timeseries analysis

#=========================================================================
# INPUT OPTION PARSER and SET CONSTANTS
#=========================================================================
def parse_inputs():
    """
    Parse inputs for optimizer.py.
    Inputs:
    Outputs:
    """
    parser = OptionParser()
    parser.add_option('-o', '--optimize', help = 'Optimizer option. Default = False.  When optimize = False, only property prediction is performed. When optimize = True, property prediction and optimization are performed.', default = False, action = 'store_true')
    parser.add_option('-a', '--temperature_t4900ewcoul', help = 'Temperature of "_t4900ewcoul" suffix simulations in K. Default: 300 K.', default = 300, type = 'float')
    parser.add_option('-b', '--temperature_pure', help = 'Temperature of "_pure" suffix simulations in K. Default: 298.15 K.', default = 298.15, type = 'float')
    parser.add_option('-e', '--temperature_idh2o', help = 'Temperature of "_idh2o" suffix simulations in K. Default: 298.15 K.', default = 298.15, type = 'float')
    parser.add_option('-T', '--deltaT', help = 'Change in Temperature, default = 1 K', default = 1, type = 'float') 
    parser.add_option('-P', '--pressure', help = 'Pressure of simulations in bar, default = 1.01325 bar', default = 1.01325, type = 'float')
    parser.add_option('-n', '--nequil', help = 'Number of equilibration points, or samples at the beginning of the file that should be ignored.  Default = 0.', default = 0, type = 'int')
    parser.add_option('-v', '--verbose', help = 'Verbose option. Default = False.', default = False, action = 'store_true')
    parser.add_option('-l', '--logfile', help = 'Location of log file. Default = properties_log.txt', default = 'properties_log.txt')
    parser.add_option('-f', '--figurefile', help = 'Location of pdf file with figures.  Default = figures.pdf', default = 'figures.pdf')
    parser.add_option('-d', '--datafile', help = 'Location of txt file with data for figures.  Default = figure_data.txt', default = 'figure_data.txt')
    parser.add_option('-p', '--pbsfilename', help = 'Name of PBS job file for rerunning (not resimulating).  If running this script multiple times simultaneously, make sure that each run is using a different pbsfile name or each run will be drastically slowed down as they have to wait for every rerun job to complete. Default = 2_rerun.sh', default = '2_rerun.sh')
    parser.add_option('-s', '--single', help = 'GROMACS compiled in single precision. Default = False.', default = False, action = 'store_true')
    parser.add_option('-g', '--gromacs', help = 'Location of GROMACS binaries. Default = /h3/n1/shirtsgroup/gromacs_46/Install_v462/bin/', default = '/h3/n1/shirtsgroup/gromacs_46/Install_v462/bin/')
    parser.add_option('-c', '--cluster', help = 'Cluster option. Default = False.', default = False, action = 'store_true')
    parser.add_option('-j', '--jobfile', help = 'location of template PBS job file for resimulating. Default = /home/bsz9ur/proteinparam/scripts/pbsfile.sh', default='/home/bsz9ur/proteinparam/scripts/pbsfile.sh')
    parser.add_option('-r', '--restart', help = 'restart optimization scheme from last recorded parameters in logfile', default = False, action = 'store_true')
    parser.add_option('-w', '--directory', help = 'the directory where you would like output files (logfile, datafile, figurefile, etc...) to be saved. Default = /bigtmp/bsz9ur/proteinparam/test_lbfgsb', default = '/bigtmp/bsz9ur/proteinparam/test_lbfgsb')

    (options, args) = parser.parse_args()
    cluster = options.cluster
    optimize = options.optimize
    verbose = options.verbose
    temp_t4900ewcoul = options.temperature_t4900ewcoul
    temp_pure = options.temperature_pure
    temp_idh2o = options.temperature_idh2o
    delT = options.deltaT
    press = options.pressure * 1e5 # Convert pressure from bar to Pascals (kg m-1 s-2)
    nequil = options.nequil
    logfile = options.logfile
    figurefile = options.figurefile
    datafile = options.datafile
    pbsfilename = options.pbsfilename
    single = options.single
    gro_loc = options.gromacs
    job_script = options.jobfile
    orestart = options.restart
    directory = options.directory
    return(cluster, optimize, verbose, temp_t4900ewcoul, temp_pure, temp_idh2o, delT, press, nequil, logfile, figurefile, datafile, pbsfilename, single, gro_loc, job_script, orestart, directory)

#==============================================================================
def set_parameters(single, gro_loc):
    """
    User set parameters for optimizer.py
    Inputs:
        single : TRUE or FALSE
        gro_loc: string (full directory path of GROMACS bin)
    Outputs: prefixes, nsamples, cwd, dir_t4900ewcoul, dir_pure, dir_inputs_t4900ewcoul, dir_inputs_pure, dir_tops, g_energy, grompp, mdrun, ke_call, pe_call, v_call, n_systems, energy_suffix, weight, pbs_lines, maxcount, errortol, error, rr_states, not_rr_states
    """

    prefixes = ['sergaf_t4900ewcoul', 'thrgaf_t4900ewcoul','thr_pure', 'ser_pure']

    nsamples = 200 #This is an override that will force the program to only analyze the first nsamples points of the data! Set to 'False' if you want to use all of the data

    cwd = os.getcwd()
    cwd = cwd + '/'
 
    dir_pure = '/bigtmp/bsz9ur/proteinparam/short_pure/'
    dir_t4900ewcoul ='/bigtmp/bsz9ur/proteinparam/short_t4900ewcoul/'
    dir_idh2o = '/bigtmp/bsz9ur/proteinparam/short_idh2o/'
    dir_inputs_pure = '/bigtmp/bsz9ur/proteinparam/short_pure/'
    dir_inputs_t4900ewcoul = '/bigtmp/bsz9ur/proteinparam/short_t4900ewcoul/'
    dir_inputs_idh2o = '/bigtmp/bsz9ur/proteinparam/short_idh2o/'
    dir_tops = '/bigtmp/bsz9ur/proteinparam/inputs/topologies/'

    # Gromacs command names
    if single:
        g_energy = gro_loc + '/g_energy'
        grompp = gro_loc + '/grompp'
        mdrun = gro_loc + '/mdrun'
    else:
        g_energy = gro_loc + '/g_energy_d'
        grompp = gro_loc + '/grompp_d'
        mdrun = gro_loc + '/mdrun_d'

    # The number that is fed to g_energy to get kinetic energy
    ke_call = 'Kinetic-En.'
    v_call = 'Volume'
    pe_call = 'Potential'

    energy_suffix = '.dhdl.xvg' # Suffix of the energy file(s) created by simulation.

    weight = dict([('idChemPot', 1.), ('pChemPot', 1/5.),('dens', 1/5.), ('hCap', 1/5.), ('hVap', 1/5.), ('kappa', 1/5.), ('idh2oChemPot', 0.)])

    pbs_lines = list()
    pbs_lines.append('#!/bin/bash\n')
    #pbs_lines.append('#PBS -W group_list=shirtsgroup\n')
    #pbs_lines.append('#PBS -q shirtscluster\n')
    pbs_lines.append('\n')
    pbs_lines.append('\n')
    pbs_lines.append('#PBS -l walltime=05:00:00\n')
    pbs_lines.append('#PBS -l select=1:mpiprocs=1:ncpus=4\n')
    pbs_lines.append('#PBS -r n\n')
    #pbs_lines.append('#PBS -m abe\n')
    #pbs_lines.append('#PBS -M bsz9ur@virginia.edu\n')
    #pbs_lines.append('\n')
    pbs_lines.append('cat $PBS_NODEFILE\n')
    pbs_lines.append('LS="/jobtmp/pbstmp.$PBS_JOBID"')
    pbs_lines.append('cd $LS')
    pbs_lines.append("NP='wc -l < $PBS_NODEFILE'\n")
    pbs_lines.append("NN='sort -u $PBS_NODEFILE | wc -l'\n")
    pbs_lines.append('echo Number of nodes is $NN\n')
    pbs_lines.append('echo Number of processors is $NP\n')
    pbs_lines.append('\n')
    pbs_lines.append('export GMX_MAXBACKUP=-1\n')
    pbs_lines.append('\n')
    pbs_lines.append('export OMP_NUM_THREADS=1\n')
    pbs_lines.append('\n')
    pbs_lines.append('COMMAND GOES HERE')
    pbs_lines.append('\n')
    pbs_lines.append('sleep 2')
    pbs_lines.append('\n')
    pbs_lines.append('REMOVE COMMAND GOES HERE')
    pbs_lines.append('\n')
    pbs_lines.append('cp * %(cwd)s' % vars())


    maxcount = 100 # Maximum number of optimization loops
    errortol = 0.1 # Maximum error desired
    error = errortol + 0.1 # Initialize error (arbitrary value higher than errortol)
    objlog = list()
    objlog.append(error)

    rr_states = [0, 15]
    not_rr_states = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

    n_states = 16

    n_systems = len(prefixes)

    return(prefixes, nsamples, cwd, dir_idh2o, dir_t4900ewcoul, dir_pure, dir_inputs_idh2o, dir_inputs_t4900ewcoul, dir_inputs_pure, dir_tops, g_energy, grompp, mdrun, ke_call, pe_call, v_call, energy_suffix, weight, pbs_lines, maxcount, errortol, error, rr_states, not_rr_states, n_systems)

#===============================================================================
def set_constants():
    kB = 1.381 * 6.02214/1000.0 # Boltzmann's constant (kJ/(mol K)).
    NA = 6.022e23 # Avogadro's Constant (1/mol)
    relative_tolerance = 1e-10 # Convergence criterion of the energy estimates for BAR and MBAR.

    convertP = (((1e-9)**3) * 6.0221413e23)/1000 # convert Pa/nm^3 to kJ/mol
    convertE = 4.184 # convert kcal to kJ

    mol_mass = dict([('SOL',18.01540),('ala',16.00000),('asn',59.00000),('cys',48.00000),('gln',73.00000),('hid',82.00000),('hie',82.00000),('ile',58.00000),('leu',58.00000),('met',76.00000),('phe',92.00000),('ser',32.00000),('thr',46.00000),('trp',131.00000),('tyr',108.00000),('val',44.00000),('aceto', 58.00000), ('bdiol', 90.00000), ('bhyde', 106.00000), ('bol', 74.00000), ('bthio', 90.00000), ('cyclo', 98.00000), ('dimet', 87.00000), ('hepta', 114.00000), ('pyrid', 79.00000), ('pyrro', 67.00000), ('quino', 129.00000), ('trime', 121.00000)])

    beta_t4900ewcoul = 1./(kB * temp_t4900ewcoul) # mol/kJ
    beta_pure = 1./(kB * temp_pure) # mol/kJ
    beta_idh2o = 1./(kB * temp_idh2o) # mol/kJ

    return(kB, NA, relative_tolerance, convertP, convertE, mol_mass, beta_t4900ewcoul, beta_pure, beta_idh2o)

#=========================================================================
# EXPERIMENTAL DATA
#=========================================================================
def experimental_data(convertE, beta_idh2o):
    """Experimental data for an assortment of properties and systems
    Outputs: EXPdens, EXPdensERR, EXPhCap, EXPhCapERR, EXPkappa, EXPkappaERR, EXPidChemPot, EXPidChemPotERR, EXPhVap, EXPhVapERR, EXPpChemPot, EXPpChemPotERR"""
    # Experimental Pure Fluid density
    # Units = g/(cm3)
    # Source (all except bthio): Y. Marcus. The properties of solvents.  Wiley, 1998. 70-77
    # Source (bthio): D. R. Lide.  CRC Handbook of Chemistry and Physics 90th edition, CRC Press: Cleveland, Ohio. 2009.
    EXPdens = dict([('aceto_pure', 0.78452), ('bdiol_pure', 1.01280), ('bhyde_pure', 1.0418), ('bol_pure', 0.80600), ('bthio_pure', 0.83700), ('cyclo_pure', 0.94220), ('dimet_pure', 0.93663), ('hepta_pure', 0.81105), ('phe_pure', 0.86453), ('pyrid_pure', 0.97820), ('pyrro_pure', 0.96540), ('quino_pure', 1.08933), ('ser_pure', 0.78704), ('thr_pure', 0.78620), ('trime_pure', 0.9115)])

    EXPdensERR = dict([('aceto_pure', 0.001627), ('bdiol_pure', 0.0005), ('bhyde_pure', 0.002546), ('bol_pure', 0.000424), ('bthio_pure', 0.0005), ('cyclo_pure', 0.0005), ('dimet_pure', 0.0005), ('hepta_pure', 0.00050), ('phe_pure', 0.002554), ('pyrid_pure', 0.000566), ('pyrro_pure', 0.0005), ('quino_pure', 0.000577), ('ser_pure', 0.000359), ('thr_pure', 0.001405), ('trime_pure', 0.001556)])

    # Experimental Pure Fluid heat capacity
    # Units = J/(mol K)
    # Source (all except bthio): Y. Marcus. The properties of solvents.  Wiley, 1998. 70-77
    # Source (bthio): D. R. Lide.  CRC Handbook of Chemistry and Physics 90th edition, CRC Press: Cleveland, Ohio. 2009.
    EXPhCap = dict([('aceto_pure', 126.215), ('bdiol_pure',203.235), ('bhyde_pure',172.0), ('bol_pure',180.2825), ('bthio_pure', 172.31), ('cyclo_pure', 178.09), ('dimet_pure',176.5633), ('hepta_pure',246.425), ('phe_pure',157.323), ('pyrid_pure', 134.72), ('pyrro_pure', 127.97), ('quino_pure',192.44), ('ser_pure', 80.78), ('thr_pure', 111.948), ('trime_pure', 214.32)]) 

    EXPhCapERR = dict([('aceto_pure', 1.860), ('bdiol_pure',1.463), ('bhyde_pure',0.5), ('bol_pure',6.8797), ('bthio_pure', 1.853), ('cyclo_pure', 1.7112), ('dimet_pure',1.44015), ('hepta_pure',3.995), ('phe_pure',0.5), ('pyrid_pure', 1.2244508), ('pyrro_pure', 0.5), ('quino_pure',10.6844), ('ser_pure', 0.738), ('thr_pure', 1.813), ('trime_pure', 0.5)])

    # Experimental Pure Fluid isothermal compressibility
    # Units = 1/GPa
    # Source (all except bthio, hepta, and trime): Y. Marcus. The properties of solvents.  Wiley, 1998. 70-77
    # Source (bthio, hepta, trime): David van der Spoel, Paul J van Maaren and Carl Caleman.  GROMACS Molecule & Liquid Database. (virtual-chemistry.org)
    EXPkappa = dict([('aceto_pure', 1.324), ('bdiol_pure', 0.440), ('bhyde_pure', 0.230), ('bol_pure', 0.941), ('bthio_pure', 1.20), ('cyclo_pure', 0.662), ('dimet_pure', 0.630), ('hepta_pure', 1.06), ('phe_pure', 0.922), ('pyrid_pure', 0.715), ('pyrro_pure', 0.652), ('quino_pure', 0.440), ('ser_pure', 1.248), ('thr_pure', 1.153), ('trime_pure', 0.83), ('water', 0.45)]) 

    EXPkappaERR = dict([('aceto_pure', 0.005), ('bdiol_pure',0.005), ('bhyde_pure',0.005), ('bol_pure',0.005), ('bthio_pure', 0.005), ('cyclo_pure', 0.005), ('dimet_pure',0.005), ('hepta_pure',0.005), ('phe_pure',0.005), ('pyrid_pure', 0.005), ('pyrro_pure', 0.005), ('quino_pure',0.005), ('ser_pure', 0.005), ('thr_pure', 0.005), ('trime_pure', 0.005)])

    # Infinite Dilution Chemical Potential
    # Units = kJ/mol
    # Source: M. R. Shirts, V. S. Pande. Solvation free energies of amino acid side chain analogs for common molecular mechanics water models.J. of Chem. Phys. (2005) 122: 134508
    EXPidChemPot = dict([('alagaf_t4900ewcoul', 1.94 * convertE), ('asngaf_t4900ewcoul', -9.68 * convertE), ('cysgaf_t4900ewcoul', -1.24 * convertE), ('glngaf_t4900ewcoul', -9.38 * convertE), ('hidgaf_t4900ewcoul', -10.27 * convertE), ('hiegaf_t4900ewcoul', -10.27 * convertE), ('ilegaf_t4900ewcoul', 2.15 * convertE), ('leugaf_t4900ewcoul', 2.28 * convertE), ('metgaf_t4900ewcoul', -1.48 * convertE), ('phegaf_t4900ewcoul', -0.76 * convertE), ('sergaf_t4900ewcoul', -5.06 * convertE), ('thrgaf_t4900ewcoul', -4.88 * convertE), ('trpgaf_t4900ewcoul', -5.88 * convertE), ('tyrgaf_t4900ewcoul', -6.11 * convertE), ('valgaf_t4900ewcoul', 1.99 * convertE)])

    EXPidChemPotERR = dict([('alagaf_t4900ewcoul', 0.05 * convertE), ('asngaf_t4900ewcoul', 0.05 * convertE), ('cysgaf_t4900ewcoul', 0.052 * convertE), ('glngaf_t4900ewcoul', 0.05 * convertE), ('hidgaf_t4900ewcoul', 0.05 * convertE), ('hiegaf_t4900ewcoul', 0.05 * convertE), ('ilegaf_t4900ewcoul', 0.05 * convertE), ('leugaf_t4900ewcoul', 0.05 * convertE), ('metgaf_t4900ewcoul', 0.05 * convertE), ('phegaf_t4900ewcoul', 0.05 * convertE), ('sergaf_t4900ewcoul', 0.05 * convertE), ('thrgaf_t4900ewcoul', 0.056 * convertE), ('trpgaf_t4900ewcoul', 0.05 * convertE), ('tyrgaf_t4900ewcoul', 0.05 * convertE), ('valgaf_t4900ewcoul', 0.05 * convertE)])

    # Heat of Vaporization
    # Units = kJ/mol
    # Source (all except bthio): Y. Marcus. The properties of solvents.  Wiley, 1998. 70-77
    # Source (bthio): David van der Spoel, Paul J van Maaren and Carl Caleman.  GROMACS Molecule & Liquid Database.
    EXPhVap =  dict([('aceto_pure', 30.99), ('bdiol_pure', 76.60), ('bhyde_pure', 39.60), ('bol_pure', 52.35), ('bthio_pure', 37.14), ('cyclo_pure', 45.13), ('dimet_pure', 50.23), ('hepta_pure', 47.24), ('phe_pure', 37.99), ('pyrid_pure', 40.15), ('pyrro_pure', 45.15), ('quino_pure', 64.10), ('ser_pure', 37.43), ('thr_pure', 42.32), ('trime_pure', 50.34)]) 

    EXPhVapERR = dict([('aceto_pure', 0.05), ('bdiol_pure',0.05), ('bhyde_pure',0.05), ('bol_pure',0.05), ('bthio_pure', 0.05), ('cyclo_pure', 0.05), ('dimet_pure',0.05), ('hepta_pure',0.05), ('phe_pure',0.05), ('pyrid_pure', 0.05), ('pyrro_pure', 0.05), ('quino_pure',0.05), ('ser_pure', 0.05), ('thr_pure', 0.05), ('trime_pure', 0.05)])

    # Pure fluid chemical potential
    # Units = kJ/mol
    # Source (Vapor Pressure Eqation and Coefficients): C. L. Yaws. Yaw's Handbook of Thermodynamic and Physical Properties of Chemical Compounds. Knovel: Beaumont, Texas. 2003.
    # Source (Equation relating vapor pressure to free energy of solvation): P. Winget, G. D. Hawkins, C. J. Cramer, D. G. Truhlar. Prediction of Vapor Pressures from Self-Solvation Free Energies Calculated by the SM5 Series of Universal Solvation Models
    EXPpChemPot = dict([('aceto_pure', -17.3490), ('bdiol_pure', -41.6733), ('bhyde_pure', -29.6180), ('bol_pure', -25.4437), ('bthio_pure', -20.4281), ('cyclo_pure', -26.8316), ('dimet_pure', -28.5365), ('hepta_pure', -25.9244), ('phe_pure', -21.6207), ('pyrid_pure', -23.0807), ('pyrro_pure', -25.7154), ('quino_pure', -36.6153), ('ser_pure', -20.3143), ('thr_pure', -21.2850), ('trime_pure', -27.6789)]) 

    EXPpChemPotERR = dict([('aceto_pure', 0.012667), ('bdiol_pure', 0.003338), ('bhyde_pure', 0.005784), ('bol_pure', 0.007793), ('bthio_pure', 0.004267), ('cyclo_pure', 0.010101), ('dimet_pure', 0.003212), ('hepta_pure', 0.089138), ('phe_pure', 0.026587), ('pyrid_pure', 0.004198), ('pyrro_pure', 0.029554), ('quino_pure', 0.011536), ('ser_pure', 0.011523), ('thr_pure', 0.005415), ('trime_pure', 0.011541)])

    # Experimental activity coefficient of water at infinite dilution in solvent X
    # Units = unitless
    # Source: DeChema
    # aceto measurement is at 307.85K, ser and thr measurements are at 298.15K
    EXPactcoeff = dict([('aceto_idh2o', 6.02), ('ser_idh2o',4.18), ('thr_idh2o',3.28)]) 
    EXPactcoeffERR = dict([('aceto_idh2o', 0.005), ('ser_idh2o',0.005), ('thr_idh2o',0.005)])

    #Calculate Experimental ID free energy of solvation from experimental activity coefficient
    # Units = kJ/mol
    EXPidh2oChemPot=dict()
    EXPidh2oChemPotERR=dict()
    for molecule in EXPactcoeff:
        EXPidh2oChemPot[molecule] = np.log(EXPactcoeff[molecule])/beta_idh2o
        EXPidh2oChemPotERR[molecule] = np.abs(EXPactcoeffERR[molecule]/(EXPactcoeff[molecule]*beta_idh2o))
    return(EXPdens, EXPdensERR, EXPhCap, EXPhCapERR, EXPkappa, EXPkappaERR, EXPidChemPot, EXPidChemPotERR, EXPhVap, EXPhVapERR, EXPpChemPot, EXPpChemPotERR, EXPidh2oChemPot, EXPidh2oChemPotERR)

#=========================================================================
# CHECK FOR REQUIRED FILES
#=========================================================================
def file_check(prefixes, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, dir_tops, optimize):
    """
    Check that all required file types are present before beginning optimization.  This saves time in debugging.
    If not optimizing then the following file types are needed:
    for each prefix: 
    .top topology file (dir + prefix + .top)
    .itp (dir + name + .itp)
    .dhdl.xvg for every state (dir + prefix + .state.dhdl.xvg
    .edr for every state (dir + prefix + .state.edr)
    GROMACS binaries
    
    If optimizing then the following additional file types are needed:
    gaffnonbonded.itp (dir_tops + gaffnonbonded.itp)
    gaforcefield (dir_inputs + gaforcefield.itp)
    for every prefix:
    .top (dir_inputs + prefix + .top)
    .mdp for every state (dir + prefix +.state.mdp)
    .gro for every state (dir + prefix +.state.gro)
    .trr for every state (dir + prefix +.state.trr)
    
    
    """
    n_states = 16 # This should be automated later
    missing_files = 0
    # Check location of GROMACS binaries, assume that if the path exists then all binaries are present. 
    if not os.path.exists(gro_loc):
        print 'GROMACS BINARY LOCATION IS INCORRECT, CHECK -g FLAG\n'
        missing_files += 1
    # Check for simulation files that are required
    for prefix in prefixes:
        res_num = 0

        #Determine simulation type (one molecule in sea of tip4pew water, or pure fluid, or one molecule of water in a sea of another molecule)
        sim_type = dict()
        type = prefix.split('_')[-1]
        sim_type[prefix] = type

        if type == 'pure':
            name = prefix.split('_')[0]
            dir = dir_pure + name + '/'
            dir_inputs = dir_inputs_pure + name + '/'
        elif type == 'idh2o':
            name = prefix.split('_')[0]
            dir = dir_idh2o + name + '/'
            dir_inputs = dir_inputs_idh2o + name + '/'
        elif type == 't4900ewcoul':
            ### This name command is not very resilient, find a way to make it more generic
            name = prefix.split('_')[0][:3]
            dir = dir_t4900ewcoul + name + '/'
            dir_inputs = dir_inputs_t4900ewcoul + name + '/'
        
        file_name = dir + prefix + '.top'
        if not os.path.exists(file_name):
            print 'MISSING TOPOLOGY FILE FOR %(prefix)s \n' % vars()
            print 'MAKE SURE THAT %(file_name)s EXISTS AND IS IN THE CORRECT LOCATION \n' % vars()
            missing_files += 1

        file_name = dir + name + '.itp'
        if not os.path.exists(file_name):
            print 'MISSING TOPOLOGY FILE FOR %(prefix)s \n' % vars()
            print 'MAKE SURE THAT %(file_name)s EXISTS AND IS IN THE CORRECT LOCATION \n' % vars()
            missing_files += 1

        for state in range(n_states):

            file_name = dir + prefix + '.' + str(state) + '.edr'
            if not os.path.exists(file_name):
                print 'MISSING ENERGY FILE FOR %(prefix)s.%(state)s \n' % vars()
                print 'MAKE SURE THAT %(file_name)s EXISTS AND IS IN THE CORRECT LOCATION \n' % vars()
                missing_files += 1
   
            file_name = dir + prefix + '.' + str(state) + energy_suffix
            if not os.path.exists(file_name):
                print 'MISSING DHDL FILE FOR %(prefix)s.%(state)s \n' % vars()
                print 'MAKE SURE THAT %(file_name)s EXISTS AND IS IN THE CORRECT LOCATION \n' % vars()
                missing_files += 1

        if optimize:

            for state in range(n_states):
                file_name = dir + prefix + '.' + str(state) + '.mdp'
                if not os.path.exists(file_name):
                    print 'MISSING SIMULATION INPUT FILE FOR %(prefix)s.%(state)s \n' % vars()
                    print 'MAKE SURE THAT %(file_name)s EXISTS AND IS IN THE CORRECT LOCATION \n' % vars()
                    missing_files += 1

                file_name = dir + prefix + '.' + str(state) + '.gro'
                if not os.path.exists(file_name):
                    print 'MISSING COORDINATE FILE FOR %(prefix)s.%(state)s \n' % vars()
                    print 'MAKE SURE THAT %(file_name)s EXISTS AND IS IN THE CORRECT LOCATION \n' % vars()
                    missing_files += 1

                file_name = dir + prefix + '.' + str(state) + '.trr'
                if not os.path.exists(file_name):
                    print 'MISSING TRAJECTORY FILE FOR %(prefix)s.%(state)s \n' % vars()
                    print 'MAKE SURE THAT %(file_name)s EXISTS AND IS IN THE CORRECT LOCATION \n' % vars()
                    missing_files += 1

    file_name = dir_tops + 'gaffnonbonded.itp'
    if not os.path.exists(file_name):
        print 'MISSING FORCE FIELD TOPOLOGY FILE'
        print 'MAKE SURE THAT %(file_name)s EXISTS AND IS IN THE CORRECT LOCATION \n' % vars()
        missing_files += 1

    file_name = dir_inputs + 'gaforcefield.itp'
    if not os.path.exists(file_name):
        print 'MISSING FORCE FIELD TOPOLOGY FILE'
        print 'MAKE SURE THAT %(file_name)s EXISTS AND IS IN THE CORRECT LOCATION \n' % vars()
        missing_files += 1

    if missing_files != 0:
        sys.exit()
    return()


#=========================================================================
# INITIALIZE LOG FILE
#=========================================================================
def initialize_log_file(cwd, logfile):
    """
    """
    logfile = cwd + logfile
    log = open(logfile, 'w')
    log.write(datetime.datetime.now().ctime())
    log.write('\n=========================================================================\n')
    log.close()

    return()

#===============================================================================
# INITIAL PARAMETERS
#===============================================================================
def initial_parameters():
    """ Sets initial parameters for optimizer.py.
    Inputs: None
    Outputs: paramkey, initparam, numparam
    """
    # Initial parameter space includes charge scaling (first parameter) and sigma and epsilon for each atom type)
    paramkey = ('charge scaling', '004 sigma', '004 epsilon', '015 sigma', '015 epsilon', '026 sigma', '026 epsilon', '053 sigma', '053 epsilon') #where numbers refer to the GAFF molecule type numbering, availsable as reference in gaffnonbonded.itp.
    initparam = np.zeros((1,9))
    initparam[0] = [1, 3.39967e-01, 4.57730e-01, 3.06647e-01, 8.80314e-01, 2.64953e-01, 6.56888e-02, 2.47135e-01, 6.56888e-02 ]

    numparam = np.shape(initparam)[1]
    return(paramkey, initparam, numparam)

#=========================================================================
# GENERATE MATRIX OF PERTURBED PARAMETERS
#=========================================================================
def perturb_params(optimize, paramkey, initparam, numparam, current_params):
    """ Perturb parameters for optimizer.py
    Inputs: optimize (TRUE or FALSE), initparam, numparam, current_params
    Outputs:initperturb, perturbparam, combos
    """
    initperturb = np.ones(numparam)
    if paramkey[0]=='charge scaling':
        initperturb[0] = initperturb[0] * 1e-2
        initperturb[1:] = initperturb[1:] * 1e-5
    else:
        initperturb = np.ones(numparam) * 1e-5
    if not(optimize):
        perturbparam = initparam
        maxcount = 0
    else:
        perturbparam = np.zeros((numparam, 1+2*numparam))
        plusperturb = current_params * np.ones((numparam, numparam))
        minusperturb = current_params * np.ones((numparam, numparam))
        for param in range(numparam):
            plusperturb[param,param] += initperturb[param]
            minusperturb[param,param] -= initperturb[param]
            if minusperturb[param,param] < 0:
                minusperturb[param,param] = 0
            if plusperturb[param,param] < 0:
                 plusperturb[param,param] = 0
        if isinstance(current_params, list) or (np.shape(current_params) == (1,9)):
            perturbparam = np.concatenate((current_params, plusperturb, minusperturb), axis = 0)
        else:
            perturbparam = np.concatenate((current_params.reshape((1,len(current_params))), plusperturb, minusperturb), axis = 0)
    combos = np.shape(perturbparam)[0]
    return(initperturb, perturbparam, combos)

#=========================================================================
# READ TOPOLOGY FILE, DETERMINE TOTAL MASS OF SYSTEM AND SIMULATION TYPE
#=========================================================================
def read_top(prefixes, dir_pure, dir_t4900ewcoul, dir_idh2o, mol_mass):
    """
    Read a set of topoly files, determine the total mass of several systems and the simulation type
    Inputs: prefixes, dir_pure, dir_t4900ewcoul, mol_mass
    Outputs: tmass, sim_type, molnum, 
    """
    tmass = dict()
    sim_type = dict()
    molnum = dict()
    for prefix in prefixes:
        res_num = 0

        #Determine simulation type (one molecule in sea of tip4pew water, or pure fluid, or one molecule of water in a sea of another molecule)
        type = prefix.split('_')[-1]
        sim_type[prefix] = type

        if type == 'pure':
            name = prefix.split('_')[0]
            dir = dir_pure + name + '/'
            tmass[prefix] = mol_mass[name]
        elif type == 'idh2o':
            name = prefix.split('_')[0]
            dir = dir_idh2o + name + '/'
            tmass[prefix] = mol_mass[name]
        elif type == 't4900ewcoul':
            ### This name command is not very resilient, find a way to make it more generic
            name = prefix.split('_')[0][:3]
            dir = dir_t4900ewcoul + name + '/'
            ### may be able to delete this
            tmass[prefix] = mol_mass['SOL']

        topfile = dir + prefix + '.top'
        topin = open(topfile, 'r')
        lines = filter(None, (line.rstrip() for line in topin))
        topin.close()

        past_molecules = False
        for line in lines:
            if line == '[ molecules ]':
                past_molecules = True
                continue
            if past_molecules:
                if line[0] != ';':
                    res_num = res_num + int(line.split()[1])

            molnum[prefix] = res_num
    return(tmass, sim_type, molnum)


#=========================================================================
# GENERATE DELTA P
#=========================================================================
def gen_delta_p(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, v_call, g_energy, EXPkappa, verbose):
    """
    Generates appropriate deltaP for reweighting for a list of systems (given by prefixes)
    inputs: prefixes, sim_type, dir_pure, dir_t4900ewcoul, v_call, g_energy, EXPkappa
    outputs: deltaPdict, max_time
    """
    deltaPdict = dict()
    max_time = dict()
    for prefix in prefixes:

        #Determine simulation type (one molecule in sea of tip4pew water, or pure fluid)
        type = sim_type[prefix]

        if type == 'pure':
            dir = dir_pure + prefix.split('_')[0] + '/'
        if type == 'idh2o':
            dir = dir_idh2o + prefix.split('_')[0] + '/'
        elif type == 't4900ewcoul':
            ### This dir command is not very resilient, find a way to make it more generic
            dir = dir_t4900ewcoul + prefix.split('_')[0][:3] + '/'

        volumefile = '%(prefix)s_volume.0.xvg' % vars()
        if not (os.path.exists(volumefile)):
            if verbose:
                print 'creating %(volumefile)s with GROMACS g_energy command' % vars()
            edr = '%(dir)s%(prefix)s.0.edr' % vars()
            if verbose:
                gen_energy_command = 'echo %(v_call)s | %(g_energy)s -f %(edr)s -o %(volumefile)s -dp' % vars()
            else:
                gen_energy_command = 'echo %(v_call)s | %(g_energy)s -f %(edr)s -o %(volumefile)s -dp &>/dev/null' % vars()
            #os.system(gen_energy_command)
            subprocess.call(gen_energy_command, shell=True)
            while not (os.path.exists(volumefile)):
                time.sleep(2)
        v = np.genfromtxt(volumefile , skip_header = 19)
        mtime = int(v[-1,0])
        volume = v[:,1]
        deltaV = np.std(volume)
        meanV = np.mean(volume)
        if (prefix.split('_')[-1] == 't4900ewcoul') or (prefix.split('_')[-1] == 'idh2o'):
            alpha = EXPkappa['water']
        else:
            alpha = EXPkappa[prefix]
        deltaP = deltaV/(alpha * meanV) # units are GPa
        deltaP = deltaP * 100000000 # units are pascals
        deltaPdict[prefix] = deltaP
        max_time[prefix] = mtime
    return (deltaPdict, max_time)

#=========================================================================
# GENERATE DICTIONARY OF MOLECULE TYPES
#=========================================================================
def gen_mol_types(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o):
    """
    Generate a dictionary of molecule types in each system for a  set of systems.
    Inputs: prefixes, sim_type, dir_pure, dir_t4900ewcoul, 
    Outputs:mol_type_dict
    """
    mol_type_dict = dict()
    for prefix in prefixes:
        mol_type_list = list()
        type = sim_type[prefix]
        if type == 'pure':
            mol_name = prefix.split('_')[0]
            dir = dir_pure + mol_name + '/'
        elif type == 'idh2o':
            mol_name = prefix.split('_')[0]
            dir = dir_idh2o + mol_name + '/'
        elif type == 't4900ewcoul':
            ### This dir command is not very resilient, find a way to make it more generic
            mol_name = prefix.split('_')[0][:3]
            dir = dir_t4900ewcoul + mol_name + '/'

        itpfile = dir + mol_name + '.itp'
        itpin = open(itpfile, 'r')
        lines = itpin.readlines()
        itpin.close()
        record = 0
        for line in lines:
            elements = line.split()
            if line == '[ atoms ]\n':
                record = 1
            if record == 1:
                if (len(elements)>1) and (elements[1].split('_')[0] == 'gaff'):
                    number = elements[1].split('_')[1]
                    skip = 0
                    for recorded in mol_type_list:
                        if number == recorded:
                            skip = 1
                    if skip == 0:
                        mol_type_list.append(number)
            if line == '[ bonds ]\n':
                break
        mol_type_dict[prefix] = mol_type_list
    return (mol_type_dict)
    
#=========================================================================
# READ ALL SNAPSHOT DATA (META DATA)
#=========================================================================
def gen_meta_data(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, energy_suffix, nsamples, verbose):
    """
    Generate meta data for reading dhdl files
    inputs: prefix, sim_type, dir_pure, dir_t4900ewcoul, 
    outputs:
    """
    # Assuming that all simulation sets have the same meta data
    prefix = prefixes[0]
    type = sim_type[prefix]

    if type == 'pure':
        dir = dir_pure + prefix.split('_')[0] + '/'
    if type == 'idh2o':
        dir = dir_idh2o + prefix.split('_')[0] + '/'
    elif type == 't4900ewcoul':
        ### This dir command is not very resilient, find a way to make it more generic
        dir = dir_t4900ewcoul + prefix.split('_')[0][:3] + '/'
    # Generate a list of all .xvg datafiles, sort them numerically if necessary and count them
    filenames = getoutput('ls %(dir)s%(prefix)s*%(energy_suffix)s' % vars()).split()
    def sortbynum(item):
        vals = item.split('.')
        for v in reversed(vals):
            if v.isdigit():
                return int(v)
        print "Warning: No digits found in filename can't sort", item
    filenames.sort(key=sortbynum)
    n_files = len(filenames)

    nsnapshots = np.zeros(n_files, int) # nsnapshots[nf] is the number of snapshots from file nf, no larger than the number of lines.
    # Temporarily read file into memory
    for nf in range(n_files):
        infile = open(filenames[nf], 'r')
        lines = infile.readlines()
        infile.close()

        #Determine the maxnumber of snapshots from quickly parsing the file and ignoring the head lines.
        nsnapshots[nf] = 0
        for line in lines:
            if ((line[0] == '#') or (line[0] == '@')):
                continue
            nsnapshots[nf] += 1
        if nsnapshots[nf] == 0:
            pdb.set_trace()

    # Determine the maximum number of snapshots from any file, assume that this is the max for any state.
    maxn = max(nsnapshots)
    if nsamples:
        maxn = nsamples

    # Read in prefixes from first file
    filename = filenames[0]
    if verbose:
        print 'Reading metadata from %(filename)s...' % vars()
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    # Initialize variables to be set
    n_components = 0
    n_states = 0
    bPV = False
    bExpanded  = False
    bEnergy = False

    maxlambda = 1000 # arbitrary number higher than the maximum number of lambdas
    for line in lines:
        # Split each line into elements
        elements = line.split()

        # This section auto parses the header to count the number of dhdl components, and the number of states at which energies are calculated.
        # This may have to be modified for different file input formats

        if ((line[0] == '#') or (line[0] == '@')):
            if (line[0] == '@'):
                # It's and xvg legend entry -- load in the information
                if (line[2] == 's') and (line[3] != 'u'):
                    # It's an xvg entry, and potentially a lambda component, note it
                    if re.search('Thermodynamic state', line):
                    # This is an expanded ensemble calculation, since thermodynamic state is being tracked
                        bExpanded = True
                        print 'Thermodynamic State is listed as a variable.'
                        print 'This is an expanded ensemble simulation.'
                    elif re.search('Energy', line):
                        # Keeping track of the total energy
                        bEnergy = True
                    elif re.search('-lambda', line):
                        # This is a listing of the lambdas
                        n_components += 1
                        lv = np.zeros([maxlambda, n_components], float)
                    elif re.search("\\\\xD\\\\f", line):
                        lambdaone = re.sub('[(),"]', '', elements[6])
                        lambdatwo = re.sub('[(),"]', '', elements[7])
                        lambdas = [lambdaone, lambdatwo]
                        for i in range(n_components):
                            lv[n_states,i] = float(lambdas[i])
                        n_states += 1
                    elif re.search('pV', line):
                        # Tracking pV term, constant pressure simulations
                        bPV = True;
        else:
            # Finished with the metadata, exit.
            break
    if verbose:
        print 'Completed reading metadata from %(filename)s...' % vars()
    return(nsnapshots, n_files, n_components, n_states, maxn, bPV, bExpanded, bEnergy, lv)

#=========================================================================
# SET SYSTEM VARIABLES
#=========================================================================
def set_system_variables(prefix, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, beta_pure, beta_t4900ewcoul, beta_idh2o, temp_t4900ewcoul, temp_pure, temp_idh2o, deltaPdict, delT, kB):
    # May have to modify this to get temperature from the mdp files
    if sim_type[prefix] == 'pure':
        dir = dir_pure + prefix.split('_')[0] + '/'
        dir_inputs = dir_inputs_pure + prefix.split('_')[0] + '/'
        temp = temp_pure
        beta = beta_pure
    elif sim_type[prefix] == 'idh2o':
        dir = dir_idh2o + prefix.split('_')[0] + '/'
        dir_inputs = dir_inputs_idh2o + prefix.split('_')[0] + '/'
        temp = temp_idh2o
        beta = beta_idh2o
    elif sim_type[prefix] == 't4900ewcoul':
        ### This dir command is not very resilient, find a way to make it more generic
        dir = dir_t4900ewcoul + prefix.split('_')[0][:3] + '/'
        dir_inputs = dir_inputs_t4900ewcoul + prefix.split('_')[0][:3] + '/'
        temp = temp_t4900ewcoul
        beta = beta_t4900ewcoul

    delta_beta = beta - 1/(kB * (temp+delT))
    delP = deltaPdict[prefix]
    return(dir, dir_inputs, temp, beta, delta_beta, delP)
    
#=========================================================================
# EDIT TOPOLOGY FILES WITH PERTURBED PARAMETERS
#=========================================================================
def edit_top(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, cwd, paramkey, initparam, perturbparam, dir_tops, combos, inner_loop):
    """ Edits gaffnonbonded.itp with new sigma and epsilon values, saves to a new file
    Inputs: parameters, current_params, dir_tops, newitp
    Outputs: None (writes a file)
    """
    for combo in range(combos):
        paramset = perturbparam[combo]
        tempparamkey = paramkey
        if combo == 0 and (paramset == initparam).all():
            continue
        else:
            if tempparamkey[0] == 'charge scaling':
                for prefix in prefixes:
                    if sim_type[prefix] == 'pure':
                        origitp = dir_pure + prefix.split('_')[0] + '/' + prefix.split('_')[0] + '.itp'
                        newitp = 'new_' + prefix + str(combo) + '.itp'
                    elif sim_type[prefix] == 'idh2o':
                        origitp = dir_idh2o + prefix.split('_')[0] + '/' + prefix.split('_')[0] + '.itp'
                        newitp = 'new_' + prefix + str(combo) + '.itp'
                    elif sim_type[prefix] == 't4900ewcoul':
                        ### This dir command is not very resilient, find a way to make it more generic
                        origitp = dir_t4900ewcoul + prefix.split('_')[0][:3] + '/' + prefix.split('_')[0][:3] + '.itp'
                        newitp = 'new_' + prefix + str(combo) + '.itp'
                    itpin = open(origitp, 'r')
                    lines = itpin.readlines()
                    itpin.close()
                    itpout = open(newitp, 'w')
                    section = False
                    for line in lines:
                        elements = line.split()
                        if line == '[ atoms ]\n':
                            section = True
                        if line == '[ bonds ]\n':
                            section = False
                        if section:
                            if len(elements)==8 and elements[0] != ';':
                                newcharge = float(elements[6]) * paramset[0]
                                elements[6]= '%.5f' %newcharge
                                a = elements[0]
                                b = elements[1]
                                c = elements[2]
                                d = elements[3]
                                e = elements[4]
                                f = elements[5]
                                g = elements[6]
                                h = elements[7]
                                newline = a.rjust(6) + b.rjust(13) + c.rjust(7) + d.rjust(7) + e.rjust(7) + f.rjust(7) + g.rjust(11) + h.rjust(11) + '\n'
                                itpout.write(newline)
                            else:
                                itpout.write(line)
                        elif section == False:
                            itpout.write(line)
                    itpout.close()
                paramset = np.delete(paramset, 0)
                tempparamkey = np.delete(tempparamkey, 0)

            origitp = dir_tops + 'gaffnonbonded.itp'
            itpin = open(origitp, 'r')
            lines = itpin.readlines()
            itpin.close()
            if inner_loop==True:
                newitp = '%(cwd)sgaffnonbonded.%(combo)s.itp' % vars()
            else:
                newitp = '%(cwd)sgaffnonbonded.opt.itp' % vars()
            itpout = open(newitp, 'w')

            for line in lines:

                elements = line.split()
                if elements[0].split('_')[0] == 'gaff':
                    mnum = elements[0].split('_')[-1]

                    for number in range(len(tempparamkey)):
                        parameter = tempparamkey[number]
                        pnum = parameter.split()[0]
                        pvar = parameter.split()[-1]
                        if mnum == pnum:
                            if pvar == 'sigma':
                                sigma = paramset[number]
                                elements[5] = '%.5e' %sigma
                            if pvar == 'epsilon':
                                epsilon = paramset[number]
                                elements[6] = '%.5e' %epsilon

                    a = elements[0]
                    b = elements[1]
                    c = elements[2]
                    d = elements[3]
                    e = elements[4]
                    f = elements[5]
                    g = elements[6]
                    newline = a.rjust(9) + b.rjust(5) + c.rjust(12) + d.rjust(8) + e.rjust(3) + f.rjust(14) + g.rjust(13) + '\n'
                    itpout.write(newline)
                else:
                    itpout.write(line)
            itpout.close()
    return()

#=========================================================================
# RERUN WITH PERTURBED PARAMETERS
#=========================================================================
def rerun_perturbed(prefixes, sim_type, cwd, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, combos, paramkey, perturbparam, initparam, mol_type_dict, cluster, optimize, pbs_lines, rr_states, mdrun, grompp, max_time, nsamples, verbose, pbsfilename):
    """
    rerun simulations at new parameters (use edit_top to edit topology with new parameters)
    Inputs: prefixes, sim_type, cwd, dir_pure, dir_t4900ewcoul, dir_inputs_pure, dir_inputs_t4900ewcoul, new_gaffnonbonded, combos, paramkey, current_params, mol_type_dict, cluster (TRUE or FALSE), optimize (TRUE    or FALSE), pbs_lines, rr_states, mdrun, grompp
    Outputs:combo_dict
    """
    rerunjobs = 0
    combo_dict = dict()
    for prefix in prefixes:
        # Set variables that depend on simulation type
        if sim_type[prefix] == 'pure':
            dir = dir_pure + prefix.split('_')[0] + '/'
            dir_inputs = dir_inputs_pure + prefix.split('_')[0] + '/'
        elif sim_type[prefix] == 'idh2o':
            dir = dir_idh2o + prefix.split('_')[0] + '/'
            dir_inputs = dir_inputs_idh2o + prefix.split('_')[0] + '/'
        elif sim_type[prefix] == 't4900ewcoul':
            ### This dir command is not very resilient, find a way to make it more generic
            dir = dir_t4900ewcoul + prefix.split('_')[0][:3] + '/'
            dir_inputs = dir_inputs_t4900ewcoul + prefix.split('_')[0][:3] + '/'
        #Determine
        combo_list = list()
        combo_list.append(0)
        for combo in range(combos):
            parameter = paramkey[(combo-1)%np.shape(current_params)[1]]

            if combo == 0 and (perturbparam[0,:] == initparam).all():
                # Switch to original molecule.itp files and gaffnonbonded.itp file is not necessary as long as you are starting out with the correct files referenced.  May want to add that in to the code somewhere.
                continue
            elif combo == 0:
                continue
            else:
                if paramkey[0] == 'charge scaling':
                    newchargeitp = cwd + 'new_' + prefix + str(combo) + '.itp'

                    # Rewrite molecule.top to call the correct molecule itp file
                    origtop = dir_inputs + prefix + '.top'
                    topin = open(origtop, 'r')
                    lines = topin.readlines()
                    topin.close()
                    newtop = cwd + prefix + '.top'
                    topout = open(newtop, 'w')
                    writenewline = False
                    for line in lines:
                        if writenewline:
                            newline = '#include "%(newchargeitp)s"\n' % vars()
                            topout.write(newline)
                            writenewline = False
                            continue
                        if line =='; Include molecule topology:\n':
                            writenewline = True
                        else:
                            writenewline = False
                        topout.write(line)
                    topout.close()
                    subprocess.call('mv %(newtop)s %(origtop)s' % vars(), shell = True)
                """        
                else:
                    if paramkey[0] == 'charge scaling':
                        origtop = dir_inputs + prefix + '.top'
                        subprocess.call('mv originaltop.txt %(origtop)s' % vars(), shell = True)
                """
                new_gaffnonbonded = '%(cwd)sgaffnonbonded.%(combo)s.itp' % vars()
                # Rewrite gaforcefield to call the correct gaffnonbonded.itp
                origitp = dir_inputs + 'gaforcefield.itp'
                itpin = open(origitp, 'r')
                lines = itpin.readlines()
                itpin.close()
                newitp = cwd + 'gaforcefield.itp'
                itpout = open(newitp, 'w')
                for line in lines:
                    if len(line.split()) > 0 and line.split()[0] == '#include':
                        newline = '#include "%(new_gaffnonbonded)s"\n' % vars()
                        itpout.write(newline)
                    else:
                        itpout.write(line)
                itpout.close()
                subprocess.call('mv %(newitp)s %(origitp)s' % vars(), shell = True)

                rerun = False
                if parameter == 'charge scaling':
                    rerun = True
                    combo_list.append(combo)
                # Which molecule's parameter is changing?
                for number in mol_type_dict[prefix]:
                    if number == parameter.split()[0]:
                        rerun = True
                        combo_list.append(combo)

                if rerun:
                    # rerun simulation with new parameters!
                    for state in rr_states:
                        # Rewrite mdp file to make sure that trr files are not being written
                        origmdp = '%(dir)s%(prefix)s.%(state)s.mdp' % vars()
                        mdpin = open(origmdp, 'r')
                        lines = mdpin.readlines()
                        mdpin.close()
                        newmdp = cwd + 'temp_mdp.mdp'
                        mdpout = open(newmdp, 'w')
                        for line in lines:
                            elements = line.split()
                            if len(elements)>0 and elements[0] == 'nstxout':
                                newline = 'nstxout			 = 0\n'
                                mdpout.write(newline)
                            else:
                                mdpout.write(line)
                        mdpout.close()
                        subprocess.call('mv %(newmdp)s %(origmdp)s' % vars(), shell = True)

                        dhdlout ='%(cwd)s%(prefix)s.%(state)s.%(combo)s.dhdl.xvg' % vars()
                        edrout = '%(cwd)s%(prefix)s.%(state)s.%(combo)s.edr' % vars()
                        
                        if (os.path.exists(dhdlout)) and (os.path.exists(edrout)):
                            rm_command ='rm %(edrout)s %(dhdlout)s' % vars()
                            os.system(rm_command)
                        if verbose:
                            print 'Rerunning %(prefix)s state %(state)s at parameter set %(combo)s' % vars()
                            grompp_command = '%(grompp)s -f %(origmdp)s -c %(dir)s%(prefix)s.%(state)s.gro -p %(dir_inputs)s%(prefix)s.top -o %(cwd)s%(prefix)s.%(state)s.%(combo)s.tpr -maxwarn 2' % vars()
                        else:
                            grompp_command = '%(grompp)s -f %(origmdp)s -c %(dir)s%(prefix)s.%(state)s.gro -p %(dir_inputs)s%(prefix)s.top -o %(cwd)s%(prefix)s.%(state)s.%(combo)s.tpr -maxwarn 2 &>/dev/null' % vars()
                        mdrun_command = '%(mdrun)s -ntmpi 1 -pin off -deffnm %(cwd)s%(prefix)s.%(state)s.%(combo)s -e %(edrout)s -dhdl %(dhdlout)s -rerun %(dir)s%(prefix)s.%(state)s.trr' % vars()
                        subprocess.call(grompp_command, shell = True)
                        rm_command = 'rm mdout.mdp %(prefix)s.%(state)s.%(combo)s.log %(prefix)s.%(state)s.%(combo)s.tpr %(prefix)s.%(state)s.%(combo)s.trr \#*\#' % vars()
                        rm_command = 'rm mdout.mdp \#*\#' % vars()
                        subprocess.call('rm \#*\#', shell = True)
                        if cluster:
                            pbs_lines[18] = mdrun_command
                            pbs_lines[22] = rm_command
                            #pbsfilename = 'rerun_lbfgsb.sh'
                            jobfile = open(pbsfilename, 'w')
                            for line in pbs_lines:
                                jobfile.write(line)
                            jobfile.close()
                            os.system('qsub %(pbsfilename)s' % vars())
                            rerunjobs += 1
                        else:
                            subprocess.call(mdrun_command, shell = True)
                            os.system(rm_command)
                
        combo_dict[prefix] = combo_list
    # IF RUNNING ON THE CLUSTER, CHECK IF THE FILES ARE FINISHED RUNNING!
    if cluster and optimize:
        finished = getoutput('qstat -u bsz9ur').split()
        while rerunjobs > 0:
            jobs = len(finished[27:])/11
            rerunjobs = 0
            for job in range(jobs):
                status = finished[36+job*11]
                if status == 'S':
                    pdb.set_trace()
                    kill_command ='qdel ' + finished[16+job*11][:-4]
                    os.system(kill_command)
                if finished[30+job*11] == pbsfilename[0:10]:
                    rerunjobs += 1
            time.sleep(100)
            finished = getoutput('qstat -u bsz9ur').split()
        time.sleep(2)
        rm_command = 'rm %(pbsfilename)s*' % vars()
        os.system(rm_command)
    failed_reruns = 1
    while failed_reruns != 0:
        failed_reruns = 0
        # Check to make sure that all files were completely written (bigtmp has an intermittant read write error)
        for prefix in prefixes:
            # Set variables that depend on simulation type
            if sim_type[prefix] == 'pure':
                dir = dir_pure + prefix.split('_')[0] + '/'
                dir_inputs = dir_inputs_pure + prefix.split('_')[0] + '/'
            if sim_type[prefix] == 'idh2o':
                dir = dir_idh2o + prefix.split('_')[0] + '/'
                dir_inputs = dir_inputs_idh2o + prefix.split('_')[0] + '/'
            elif sim_type[prefix] == 't4900ewcoul':
                ### This dir command is not very resilient, find a way to make it more generic
                dir = dir_t4900ewcoul + prefix.split('_')[0][:3] + '/'
                dir_inputs = dir_inputs_t4900ewcoul + prefix.split('_')[0][:3] + '/'

            for combo in combo_dict[prefix]:
                if combo == 0:
                    continue
                for state in rr_states:
                    dhdlout ='%(cwd)s%(prefix)s.%(state)s.%(combo)s.dhdl.xvg' % vars()
                    edrout = '%(cwd)s%(prefix)s.%(state)s.%(combo)s.edr' % vars()

                    rerunfailed = False
                    if not(os.path.exists(dhdlout)) or not(os.path.exists(edrout)):
                        rerunfailed = True
                    if rerunfailed == False:
                        dhdl_endtime = (getoutput('tail -n 1 %(dhdlout)s' % vars())).split()
                        if len(dhdl_endtime) > 0:
                            dhdl_endtime = float(dhdl_endtime[0])
                        else:
                            rerunfailed = True
                        if dhdl_endtime != float(max_time[prefix]) and nsamples == False:
                           rerunfailed = True
                    if rerunfailed:
                        pdb.set_trace()
                        failed_reruns += 1
                        if paramkey[0] == 'charge scaling':
                            newchargeitp = cwd + 'new_' + prefix + str(combo) + '.itp'
                        # Rewrite molecule.top to call the correct molecule itp file
                        origtop = dir_inputs + prefix + '.top'
                        topin = open(origtop, 'r')
                        lines = topin.readlines()
                        topin.close()
                        newtop = cwd + prefix + '.top'
                        topout = open(newtop, 'w')
                        writenewline = False
                        for line in lines:
                            if writenewline:
                                newline = '#include "%(newchargeitp)s"\n' % vars()
                                topout.write(newline)
                                writenewline = False
                                continue
                            if line =='; Include molecule topology:\n':
                                writenewline = True
                            else:
                                writenewline = False
                            topout.write(line)
                        topout.close()
                        subprocess.call('mv %(newtop)s %(origtop)s' % vars(), shell = True)

                        # Rewrite gaforcefield to call the correct gaffnonbonded.itp
                        origitp = dir_inputs + 'gaforcefield.itp'
                        itpin = open(origitp, 'r')
                        lines = itpin.readlines()
                        itpin.close()
                        newitp = cwd + 'gaforcefield.itp'
                        itpout = open(newitp, 'w')
                        for line in lines:
                            if len(line.split()) > 0 and line.split()[0] == '#include':
                                newline = '#include "%(cwd)sgaffnonbonded.%(combo)s.itp"\n' % vars()
                                itpout.write(newline)
                            else:
                                itpout.write(line)
                        itpout.close()
                        subprocess.call('mv %(newitp)s %(origitp)s' % vars(), shell = True)
                        if os.path.exists(dhdlout):
                            rm_command ='rm %(dhdlout)s' % vars()
                            os.system(rm_command)
                        if os.path.exists(edrout):
                            rm_command = 'rm %(edrout)s' % vars()
                            os.system(rm_command)
                        print 'ATTEMPTING AGAIN: Rerunning %(prefix)s state %(state)s at parameter set %(combo)s' % vars()
                        if verbose:
                            grompp_command = '%(grompp)s -f %(dir)s%(prefix)s.%(state)s.mdp -c %(dir)s%(prefix)s.%(state)s.gro -p %(dir_inputs)s%(prefix)s.top -o %(cwd)s%(prefix)s.%(state)s.%(combo)s.tpr -maxwarn 2' % vars()
                        else:
                            grompp_command = '%(grompp)s -f %(dir)s%(prefix)s.%(state)s.mdp -c %(dir)s%(prefix)s.%(state)s.gro -p %(dir_inputs)s%(prefix)s.top -o %(cwd)s%(prefix)s.%(state)s.%(combo)s.tpr -maxwarn 2 &>/dev/null' % vars()
                        mdrun_command = '%(mdrun)s -ntmpi 1 -pin off -deffnm %(cwd)s%(prefix)s.%(state)s.%(combo)s -e %(edrout)s -dhdl %(dhdlout)s -rerun %(dir)s%(prefix)s.%(state)s.trr' % vars()
                        subprocess.call(grompp_command, shell = True)
                        rm_command = 'rm %(prefix)s.%(state)s.%(combo)s.log %(prefix)s.%(state)s.%(combo)s.tpr %(prefix)s.%(state)s.%(combo)s.trr \#*\#' % vars()
                        subprocess.call('rm \#*\#', shell = True)
                        if cluster:
                            pbs_lines[18] = mdrun_command
                            pbs_lines[22] = rm_command
                            jobfile = open(pbsfilename, 'w')
                            for line in pbs_lines:
                                jobfile.write(line)
                            jobfile.close()
                            os.system('qsub %(pbsfilename)s' % vars())
                            rerunjobs += 1
                        else:
                            #subprocess.call(mdrun_command, shell = True)
                            os.system(rm_command)
        if cluster and optimize:
            # IF RUNNING ON THE CLUSTER, CHECK IF THE FILES ARE FINISHED RUNNING!
            finished = getoutput('qstat -u bsz9ur').split()
            print 'checking if they are finished'
            while rerunjobs > 0:
                rerunjobs = 0
                jobs = len(finished[27:])/11
                rerunjobs = 0
                for job in range(jobs):
                    status = finished[36+job*11]
                    if status == 'S':
                        kill_command ='qdel ' + finished[16+job*11][:-4]
                        os.system(kill_command)
                    if finished[30+job*11] == pbsfilename[0:10]:
                        rerunjobs += 1
                time.sleep(10)
                finished = getoutput('qstat -u bsz9ur').split()
            time.sleep(2)
            rm_command = 'rm rerun*'
            os.system(rm_command)
    return(combo_dict)

#=========================================================================
# RERUN WITH OPTIMIZED PARAMS
#=========================================================================
def rerun_opt_param(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, paramkey, n_states, cwd, grompp, mdrun, cluster, optimize, pbs_lines, verbose, max_time, pbsfilename):
    """
    Rerun simulation output at new parameters (optimized, unsampled parameters)
    Inputs:
    Outputs:
    """
    # Rewrite gaforcefield.itp to call the correct gaffnonbonded
    rerunjobs = 0
    for prefix in prefixes:

        # Set variables that depend on simulation type
        if sim_type[prefix] == 'pure':
            dir = dir_pure + prefix.split('_')[0] + '/'
            dir_inputs = dir_inputs_pure + prefix.split('_')[0] + '/'
        if sim_type[prefix] == 'idh2o':
            dir = dir_idh2o + prefix.split('_')[0] + '/'
            dir_inputs = dir_inputs_idh2o + prefix.split('_')[0] + '/'
        elif sim_type[prefix] == 't4900ewcoul':
            ### This dir command is not very resilient, find a way to make it more generic
            dir = dir_t4900ewcoul + prefix.split('_')[0][:3] + '/'
            dir_inputs = dir_inputs_t4900ewcoul + prefix.split('_')[0][:3] + '/'
    
        if paramkey[0] == 'charge_scaling':
            newchargeitp = cwd  + 'new_' + prefix + '0.itp'
            # Rewrite molecule.top to call the correct molecule itp file
            origtop = dir_inputs + prefix + '.top'
            topin = open(origtop, 'r')
            lines = topin.readlines()
            topin.close()
            newtop = cwd + prefix + '.top'
            topout = open(newtop, 'w')
            writenewline = False
            for line in lines:
                if writenewline:
                    newline = '#include "%(newchargeitp)s"\n' % vars()
                    topout.write(newline)
                    writenewline = False
                    continue
                if line =='; Include molecule topology:\n':
                    writenewline = True
                else:
                    writenewline = False
                    topout.write(line)
                topout.close()
                subprocess.call('mv %(newtop)s %(origtop)s' % vars(), shell = True)

        origitp = dir_inputs + 'gaforcefield.itp'
        itpin = open(origitp, 'r')
        lines = itpin.readlines()
        itpin.close()
        newitp = cwd + 'gaforcefield.itp'
        itpout = open(newitp, 'w')
        for line in lines:
            if len(line.split()) > 0 and line.split()[0] == '#include':
                newline = '#include "%(cwd)sgaffnonbonded.0.itp"\n' % vars()
                itpout.write(newline)
            else:
                itpout.write(line)
        itpout.close()
        subprocess.call('mv %(newitp)s %(origitp)s' % vars(), shell = True)

        for state in range(n_states):
            # rerun simulation with new parameters!
            # Rewrite mdp file to make sure that trr files are not being written
            origmdp = '%(dir)s%(prefix)s.%(state)s.mdp' % vars()
            mdpin = open(origmdp, 'r')
            lines = mdpin.readlines()
            mdpin.close()
            newmdp = cwd + 'temp_mdp.mdp'
            mdpout = open(newmdp, 'w')
            for line in lines:
                elements = line.split()
                if len(elements)>0 and elements[0] == 'nstxout':
                    newline = 'nstxout			 = 200\n'
                    mdpout.write(newline)
                else:
                    mdpout.write(line)
            mdpout.close()
            subprocess.call('mv %(newmdp)s %(origmdp)s' % vars(), shell = True)

            dhdlout = '%(prefix)s.%(state)s.p.dhdl.xvg' % vars()
            edrout = '%(prefix)s.%(state)s.p.edr' % vars()
            tprfile = '%(cwd)s/%(prefix)s.%(state)s.p.tpr' % vars()
            trrfile = '%(prefix)s.%(state)s.p.trr' % vars()
            nlogfile = '%(prefix)s.%(state)s.p.log' % vars()
            if (os.path.exists(dhdlout)) and (os.path.exists(edrout)):
                rm_command = 'rm %(edrout)s %(dhdlout)s' % vars()
                os.system(rm_command)
            #pdb.set_trace()
            if verbose:
                print 'Rerunning %(prefix)s state %(state)s at new parameters' % vars()
                grompp_command = '%(grompp)s -f %(origmdp)s -c %(dir)s%(prefix)s.%(state)s.gro -p %(dir_inputs)s%(prefix)s.top -o %(tprfile)s -maxwarn 2' % vars()
            else:
                grompp_command = '%(grompp)s -f %(origmdp)s -c %(dir)s%(prefix)s.%(state)s.gro -p %(dir_inputs)s%(prefix)s.top -o %(cwd)s%(prefix)s.%(state)s.p.tpr -maxwarn 2 &>/dev/null' % vars()
            mdrun_command = '%(mdrun)s -ntmpi 1 -pin off -s %(tprfile)s -o %(trrfile)s -g %(nlogfile)s -e %(edrout)s -dhdl %(dhdlout)s -rerun %(dir)s%(prefix)s.%(state)s.trr' % vars()
            subprocess.call(grompp_command, shell = True)
            rm_command = 'rm mdout.mdp %(prefix)s.%(state)s.p.log %(prefix)s.%(state)s.p.tpr %(prefix)s.%(state)s.p.trr \#*\#' % vars()
            subprocess.call('rm \#*\#', shell = True)
            if cluster:
                pbs_lines[18] = mdrun_command
                pbs_lines[22] = rm_command
                #pbsfilename = 'rerun_lbfgs.sh'
                jobfile = open(pbsfilename, 'w')
                for line in pbs_lines:
                    jobfile.write(line)
                jobfile.close()
                os.system('qsub %(pbsfilename)s' % vars())
                rerunjobs += 1
            else:
                subprocess.call(mdrun_command, shell = True)
                os.system(rm_command)
    # IF RUNNING ON THE CLUSTER, CHECK IF THE FILES ARE FINISHED RUNNING!
    if cluster and optimize:
        finished = getoutput('qstat -u bsz9ur').split()
        while rerunjobs > 0:
            jobs = len(finished[27:])/11
            rerunjobs = 0
            for job in range(jobs):
                status = finished[36+job*11]
                if status == 'S':
                    pdb.set_trace()
                    kill_command ='qdel ' + finished[16+job*11][:-4]
                    os.system(kill_command)
                if finished[30+job*11] == pbsfilename[0:10]:
                    rerunjobs += 1
            if rerunjobs == 0:
                break
            time.sleep(10)
            finished = getoutput('qstat -u bsz9ur').split()
        rm_command = 'rm %(pbsfilename)s*' % vars()
        os.system(rm_command)

    failed_reruns = 1
    while failed_reruns != 0:
        failed_reruns = 0
        # Check to make sure that all files were completely written (bigtmp has an intermittant read write error)
        for prefix in prefixes:
            # Set variables that depend on simulation type
            if sim_type[prefix] == 'pure':
                dir = dir_pure + prefix.split('_')[0] + '/'
                dir_inputs = dir_inputs_pure + prefix.split('_')[0] + '/'
            elif sim_type[prefix] == 't4900ewcoul':
                ### This dir command is not very resilient, find a way to make it more generic
                dir = dir_t4900ewcoul + prefix.split('_')[0][:3] + '/'
                dir_inputs = dir_inputs_t4900ewcoul + prefix.split('_')[0][:3] + '/'
            if sim_type[prefix] == 'idh2o':
                dir = dir_idh2o + prefix.split('_')[0] + '/'
                dir_inputs = dir_inputs_idh2o + prefix.split('_')[0] + '/'

                for state in range(n_states):
                    dhdlout ='%(cwd)s%(prefix)s.%(state)s.p.dhdl.xvg' % vars()
                    edrout = '%(cwd)s%(prefix)s.%(state)s.p.edr' % vars()

                    rerunfailed = False

                    if not(os.path.exists(dhdlout)) or not(os.path.exists(edrout)):
                        rerunfailed = True
                    if rerunfailed == False:
                        dhdl_endtime = (getoutput('tail -n 1 %(dhdlout)s' % vars())).split()
                        if len(dhdl_endtime) > 0:
                            dhdl_endtime = float(dhdl_endtime[0])
                        else:
                            rerunfailed = True
                        if dhdl_endtime != float(max_time[prefix]) and nsamples == False:
                           rerunfailed = True

                    if rerunfailed:
                        failed_reruns += 1
                        if paramkey[0] == 'charge scaling':
                            newchargeitp = cwd + 'new_' + prefix + 'p.itp'
                            # Rewrite molecule.top to call the correct molecule.itp file
                            origtop = dir_inputs + prefix + '.top'
                            topin = open(origtop, 'r')
                            lines = topin.readlines()
                            topin.close()
                            newtop = cwd + prefix + '.top'
                            topout = open(newtop, 'w')
                            writenewline = False
                            for line in lines:
                                if writenewline:
                                    newline = '#include "%(newchargeitp)s"\n' % vars()
                                    topout.write(newline)
                                    writenewline = False
                                    continue
                                if line =='; Include molecule topology:\n':
                                    writenewline = True
                                else:
                                    writenewline = False
                                    topout.write(line)
                            topout.close()
                            subprocess.call('mv %(newtop)s %(origtop)s' % vars(), shell = True)

                        # Rewrite gaforcefield to call the correct gaffnonbonded.itp
                        origitp = dir_inputs + 'gaforcefield.itp'
                        itpin = open(origitp, 'r')
                        lines = itpin.readlines()
                        itpin.close()
                        newitp = cwd + 'gaforcefield.itp'
                        itpout = open(newitp, 'w')
                        for line in lines:
                            if len(line.split()) > 0 and line.split()[0] == '#include':
                                newline = '#include "%(cwd)sgaffnonbonded.0.itp"\n' % vars()
                                itpout.write(newline)
                            else:
                                itpout.write(line)
                        itpout.close()
                        subprocess.call('mv %(newitp)s %(origitp)s' % vars(), shell = True)
                        if os.path.exists(dhdlout):
                            rm_command ='rm %(dhdlout)s' % vars()
                            os.system(rm_command)
                        if os.path.exists(edrout):
                            rm_command = 'rm %(edrout)s' % vars()
                            os.system(rm_command)
                        if verbose:
                            print 'ATTEMPTING AGAIN: Rerunning %(prefix)s state %(state)s at optimized parameters' % vars()
                        grompp_command = '%(grompp)s -f %(dir)s%(prefix)s.%(state)s.mdp -c %(dir)s%(prefix)s.%(state)s.gro -p %(dir_inputs)s%(prefix)s.top -o %(cwd)s%(prefix)s.%(state)s.p.tpr -maxwarn 2 &>/dev/null' % vars()
                        mdrun_command = '%(mdrun)s -ntmpi 1 -pin off -deffnm %(cwd)s%(prefix)s.%(state)s.p -e %(edrout)s -dhdl %(dhdlout)s -rerun %(dir)s%(prefix)s.%(state)s.trr' % vars()
                        subprocess.call(grompp_command, shell = True)
                        rm_command = '#rm %(prefix)s.%(state)s.p.log %(prefix)s.%(state)s.p.tpr %(prefix)s.%(state)s.p.trr \#*\#' % vars()
                        subprocess.call('rm \#*\#', shell = True)
                        if cluster:
                            pbs_lines[18] = mdrun_command
                            pbs_lines[22] = rm_command
                            jobfile = open(pbsfilename, 'w')
                            for line in pbs_lines:
                                jobfile.write(line)
                            jobfile.close()
                            os.system('qsub %(pbsfilename)s' % vars())
                            rerunjobs += 1
                        else:
                            subprocess.call(mdrun_command, shell = True)
                            os.system(rm_command)
        # IF RUNNING ON THE CLUSTER, CHECK IF THE FILES ARE FINISHED RUNNING!
        if cluster and optimize:
            finished = getoutput('qstat -u bsz9ur').split()
            while rerunjobs > 0:
                rerunjobs = 0
                jobs = len(finished[27:])/11
                rerunjobs = 0
                for job in range(jobs):
                    status = finished[36+job*11]
                    if status == 'S':
                        pdb.set_trace()
                        kill_command ='qdel ' + finished[16+job*11][:-4]
                        os.system(kill_command)
                    if finished[30+job*11] == pbsfilename[0:10]:
                        rerunjobs += 1
                time.sleep(10)
                finished = getoutput('qstat -u bsz9ur').split()
            rm_command = 'rm rerun*'
            os.system(rm_command)

    return()

#=========================================================================
# GENERATE SUB_U_KLT
#=========================================================================
def gen_u_klt(prefix, state, max_time, dhdlfile, edrfile, energyfile, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, K, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, nsamples, verbose):
    """
    Read in data from dhdl and edr files and generate a K x K x n_samples u_klt matrix
    """
    something_is_wrong = 0
    
    pe = np.zeros(K, np.float64) # temporary storage for energies
    sub_u_klt = np.zeros([1, K, maxn], np.float64)
    sub_dhdlt = np.zeros([1, n_components, maxn], float)
    sub_tp_pe =np.zeros([1, maxn], float)
    sub_tp_pv =np.zeros([1, maxn], float)
    
    # Get kinetic energy from edr file
    if bEnergy:
        # This might not be needed
        #if verbose:
        #    gen_energy_command = 'echo %(ke_call)s | %(g_energy)s -f %(edrfile)s -o %(energyfile)s -dp' % vars()
        #else:
        #    gen_energy_command = 'echo %(ke_call)s | %(g_energy)s -f %(edrfile)s -o %(energyfile)s -dp &>/dev/null' % vars()
        #subprocess.call(gen_energy_command, shell = True)
        #waiting = 0
        #while not(os.path.exists(energyfile)):
        if verbose:
            print 'creating %(energyfile)s with GROMACS g_energy' % vars()
            gen_energy_command = 'echo %(ke_call)s | %(g_energy)s -f %(edrfile)s -o %(energyfile)s -dp' % vars()
        else:
            gen_energy_command = 'echo %(ke_call)s | %(g_energy)s -f %(edrfile)s -o %(energyfile)s -dp &>/dev/null' % vars()
        #waiting +=1
        subprocess.call(gen_energy_command, shell = True)
        waiting = 0
        while not(os.path.exists(energyfile)):
            import time
            time.sleep(2)
            waiting += 1
            if waiting > 10:
                subprocess.call(gen_energy_command, shell = True)
                waiting = 0
                while not(os.path.exists(energyfile)):
                    import time
                    time.sleep(2)
                    waiting += 1
                    if waiting > 10:
                        pdb.set_trace()
        energy_maxtime = float(getoutput('tail -n 1 %(energyfile)s' % vars()).split()[0])
        if energy_maxtime != max_time[prefix] and nsamples == False:
            subprocess.call(gen_energy_command, shell = True)
        ke = np.genfromtxt(energyfile, skip_header = 19)
        ke = ke[:,1]
        if nsamples != False:
            ke = ke[:nsamples]
        if len(ke) < maxn:
            while len(ke)< maxn:
                ke = np.genfromtxt(energyfile, skip_header = 19)
                ke = ke[:,1]
                if nsamples != False:
                    ke = ke[:nsamples]
        energy_count = 0

    # Read dhdl file
    if verbose:
        print 'Reading %(dhdlfile)s ...' % vars()
    infile = open(dhdlfile, 'r')
    lines = infile.readlines()
    infile.close()
    while len(lines) == 0:
        infile = open(dhdlfile, 'r')
        lines = infile.readlines()
        infile.close()
        if len(lines)==0:
            something_is_wrong = 1
            pdb.set_trace()

    bEquil = False
    tot_nsnapshots = 0
    nsnapshots[state] = 0
    
    for line in lines:
        if ((line[0] != '#') and (line[0] != '@')):
            if not bEquil:
                if bExpanded:
                    if tot_nsnapshots > nequil:
                        tot_nsnapshots = 0 #start over again
                        bEquil = True
                    else:
                        if nsnapshots[state] > nequil:
                            nsnapshots[state] = 0
                            bEquil = True
            elements = line.split()

            # Record the time of the sample (get rid of the first element)
            time = float(elements.pop(0))

            if bExpanded:
                state = int(elements.pop(0))
            else:
                state = state

            if bEnergy:
                energy = float(elements.pop(0)) - ke[energy_count]
                energy_count += 1

                sub_tp_pe[0, nsnapshots[state]] = energy
            else:
                energy = 0
            for nl in range(n_components):
                sub_dhdlt[0, nl, nsnapshots[state]] = float(elements.pop(0))
            if something_is_wrong == 1:
                pdb.set_trace()
            for l in range(K):
                pe[l] = float(elements.pop(0)) + energy
            if bPV:
                pv = float(elements.pop(0))
                sub_tp_pv[0, nsnapshots[state]] = pv

            # Compute and store reduced potential energy at each state
            for l in range(K):
                sub_u_klt[0, l, nsnapshots[state]] = beta * (pe[l] + pv)

            if bExpanded:
                u_t[tot_nsnapshots] = sub_u_klt[state, state, nsnapshots[state]]
            tot_nsnapshots += 1
            nsnapshots[state] += 1

            if nsnapshots[state] == maxn:
                break
    if nsnapshots[state] != maxn:
        pdb.set_trace()
    return(sub_u_klt, sub_dhdlt, sub_tp_pe, sub_tp_pv)

#=========================================================================
# GENERATE U_KLT FOR INITIAL PARAMETERS (SAMPLED STATES)
#=========================================================================
def gen_u_klt_initial(prefix, n_files, dir, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, K, max_time, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, nsamples, verbose):
    """
    generate u_klt for a single set of sampled states
    Inputs:
    Outputs: u_klt, dhdlt,tp_pe0, tp_pv0, tp_pe15, tp_pv15
    """
    for state in range(n_files):
        dhdlfile = '%(dir)s%(prefix)s.%(state)s.dhdl.xvg' % vars()
        edrfile = '%(dir)s%(prefix)s.%(state)s.edr' % vars()
        energyfile = '%(prefix)s_energy.%(state)s.xvg' % vars()
        if os.path.exists(energyfile):
            os.system('rm %(energyfile)s' % vars())
        (sub_u_klt, sub_dhdlt, sub_tp_pe, sub_tp_pv) = gen_u_klt(prefix, state, max_time, dhdlfile, edrfile, energyfile, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, K, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, nsamples, verbose)
        u_klt[state, :K, :] = sub_u_klt
        dhdlt[state, :K, :] = sub_dhdlt
        if state == 0:
            tp_pe0[0, :] = sub_tp_pe
            tp_pv0[0, :] = sub_tp_pv
        if state == n_files-1:
            tp_pe15[0, :] = sub_tp_pe
            tp_pv15[0, :] = sub_tp_pv
    return(u_klt, dhdlt,tp_pe0, tp_pv0, tp_pe15, tp_pv15)

#=========================================================================
# ADD TO U_KLT FOR OPTIMIZED PARAMETERS (UNSAMPLED STATES)
#=========================================================================
def add_u_klt_optimized(prefix, n_states, initparam, current_params, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, max_time, K, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, cwd, nsamples, verbose):
    """
    """
    shift = 1
    for state in range(n_states):
        dhdlfile = '%(cwd)s%(prefix)s.%(state)s.p.dhdl.xvg' % vars()
        edrfile = '%(cwd)s%(prefix)s.%(state)s.p.edr' % vars()
        energyfile = '%(prefix)s_energy.%(state)s.p.xvg' % vars()
        if os.path.exists(energyfile):
            os.system('rm %(energyfile)s' % vars())
        (sub_u_klt, sub_dhdlt, sub_tp_pe, sub_tp_pv) = gen_u_klt(prefix, state, max_time, dhdlfile, edrfile, energyfile, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, K, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, nsamples, verbose)
        u_klt[state,shift*K:(shift+1)*K, :] = sub_u_klt[0, :, :]
        if state == 0:
            tp_pe0[shift, :] = sub_tp_pe
            tp_pv0[shift, :] = sub_tp_pv
        if state == 15:
            tp_pe15[shift, :] = sub_tp_pe
            tp_pv15[shift, :] = sub_tp_pv
    return(u_klt, tp_pe0, tp_pv0, tp_pe15, tp_pv15)

#=========================================================================
# ADD TO U_KLT FOR PERTURBED PARAMETERS (UNSAMPLED STATES)
#=========================================================================
def add_u_klt_perturbed(prefix, combo_dict, rr_states, n_states, initparam, current_params, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, K, max_time, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, cwd, nsamples, verbose):
    """
    """
    if isinstance(current_params, list):
        if (current_params != initparam):
            shift = 1
        else:
            shift = 0
    else:
        if (current_params != initparam).any():
            shift = 1
        else:
            shift = 0

    counter = 0
    for combo in combo_dict[prefix]:
        if combo == 0:
            continue
        counter += 1
        for state in rr_states:
            dhdlfile = '%(cwd)s%(prefix)s.%(state)s.%(combo)s.dhdl.xvg' % vars()
            edrfile = '%(cwd)s%(prefix)s.%(state)s.%(combo)s.edr' % vars()
            energyfile = '%(prefix)s_energy.%(state)s.%(combo)s.xvg' % vars()
            if os.path.exists(energyfile):
                os.system('rm %(energyfile)s' % vars())
            (sub_u_klt, sub_dhdlt, sub_tp_pe, sub_tp_pv) = gen_u_klt(prefix, state, max_time, dhdlfile, edrfile, energyfile, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, K, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, nsamples, verbose)
            u_klt[state,(counter+shift)*K:(counter+shift+1)*K, :] = sub_u_klt[0, :, :]
            if state == 0:
                tp_pe0[counter + shift, :] = sub_tp_pe
                tp_pv0[counter + shift, :] = sub_tp_pv
            if state == 15:
                tp_pe15[counter + shift, :] = sub_tp_pe
                tp_pv15[counter + shift, :] = sub_tp_pv
    return(u_klt, tp_pe0, tp_pv0, tp_pe15, tp_pv15)

#=========================================================================
# SUBSAMPLE DATA TO OBTAIN UNCORRELATED SAMPLES
#=========================================================================
def subsample(prefix, nsnapshots, u_klt, shift, combo_dict, K, n_components, bExpanded, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, verbose):
    """
    """
    maxn = np.max(nsnapshots)
    u_kln = np.zeros(np.shape(u_klt), np.float64) # u_kln[k,m,n] is the reduced potential energy of uncorrelated sample index n from state k evaluated at state m
    N_k = np.zeros((len(combo_dict[prefix]) + shift) * K, int) # N_k[k] is the number of uncorrelated samples from state k

    dhdl = np.zeros([K, n_components, maxn], float) #dhdl is the value for dhdl for each component in the file at each time
    
    if verbose:
        print 'Now computing correlation times'

    # Generate a dummy set to calculate indices
    if bExpanded: # Expanded ensemble uses a different correlation time
        # Use the expanded ensemble energy
        g = timeseries.statisticalInefficiency(u_t[0:tot_nsnapshots])
        g = g * np.ones(K, np.float64)
    else:
        g = np.zeros(K,float) # Autocorrelation times for the data
        for k in range(K):
            # Determine indices of uncorrelated samples from potential autocorrelation analysis at state k.
            # Alternatively, the energy differences could also be used.  Here we use the total dhdl.
            dhdl_sum = np.sum(dhdlt[k, : , 0:nsnapshots[k]], axis = 0)
            g[k] = timeseries.statisticalInefficiency(dhdl_sum)

    # Determine the indices to use
    for k in range(K):
        indices = np.array(timeseries.subsampleCorrelatedData(u_klt[k, k, 0:nsnapshots[k]], g =g[k], verbose = verbose)) # Indices of uncorrelated samples
        N = len(indices) # Number of uncorrelated samples

        for n in range(n_components):
            if not(np.shape(indices)[0]>0):
                pdb.set_trace()
            if isinstance(n, int) and isinstance(k, int) and isinstance(indices[0], int):
                pass
            else:
                pass
        dhdl[k,n,0:N] = dhdlt[k,n,indices]
        for l in range((K)):
            u_kln[k,l,0:N] = u_klt[k,l,indices]
            for combo in range(np.shape(u_klt)[1]/K):
                u_kln[k,l + (combo) * K ,0:N] = u_klt[k,l + (combo) * K,indices]

        if k == 0:
            tp_poten0 = np.zeros(N)
            tp_poten0 = tp_pe0[:,indices]
            tp_pV0 = np.zeros(N)
            tp_pV0 = tp_pv0[:,indices]
            N0 = N
        if k == 15:
            tp_poten15 = np.zeros(N)
            tp_poten15 = tp_pe15[:,indices]
            tp_pV15 = np.zeros(N)
            tp_pV15 = tp_pv15[:,indices]
            N15 = N

        N_k[k] = N
    if verbose:
         print '\nCorrelation times:'
         print g
         print ''
         print 'number of uncorrelated samples:'
         print N_k
         print ''
    return(u_kln, tp_poten0, tp_poten15, tp_pV0, tp_pV15, N_k, N0, N15)

#=========================================================================
# REWEIGHT AT T AND P
#=========================================================================
def reweight_T_P(current_params, initparam, z, prefix, N0, N15, beta, delta_beta, press, delP, tp_poten0, tp_pV0, tp_poten15, tp_pV15, combo_dict):
    """
    shift = 0 if count = 0, shift = 1 if count != 0
    """
    if isinstance(current_params, list):
        if (current_params != initparam):
            shift = 1
        else:
            shift = 0
    else:
        if (current_params != initparam).any():
            shift = 1
        else:
            shift = 0
    #tp_u_kln0 is the u_kln matrix for reweighted T and P for state 0, tp_u_kln15 is reweighted T and P for state 15
    # Changed last dimension in tp_u_kln matrices from z to N0 and N15
    tp_u_kln0 = np.zeros((5*(len(combo_dict[prefix]) + shift),5*(len(combo_dict[prefix]) + shift),N0))
    tp_N_k0 = np.zeros(5*(len(combo_dict[prefix]) + shift))
    tp_N_k0[0] = N0
    tp_u_kln15 = np.zeros((5*(len(combo_dict[prefix]) + shift),5*(len(combo_dict[prefix]) + shift),N15))
    tp_N_k15 = np.zeros(5*(len(combo_dict[prefix]) + shift))
    tp_N_k15[0] = N15

    if shift == 1:
        tp_u_kln0[0, 0,:N0] = beta * (tp_poten0[0] + tp_pV0[0])
        tp_u_kln0[0, 1,:N0] = beta * (tp_poten0[0] + tp_pV0[0] * ((press + delP) / press))
        tp_u_kln0[0, 2,:N0] = beta * (tp_poten0[0] + tp_pV0[0] * ((press - delP) / press))
        tp_u_kln0[0, 3,:N0] = (1+delta_beta) * beta * (tp_poten0[0] + tp_pV0[0])
        tp_u_kln0[0, 4,:N0] = (1-delta_beta) * beta * (tp_poten0[0] + tp_pV0[0])

        tp_u_kln15[0, 0,:N15] = beta * (tp_poten15[0] + tp_pV15[0])
        tp_u_kln15[0, 1,:N15] = beta * (tp_poten15[0] + tp_pV15[0] * ((press + delP) / press))
        tp_u_kln15[0, 2,:N15] = beta * (tp_poten15[0] + tp_pV15[0] * ((press - delP) / press))
        tp_u_kln15[0, 3,:N15] = (1+delta_beta) * beta * (tp_poten15[0] + tp_pV15[0])
        tp_u_kln15[0, 4,:N15] = (1-delta_beta) * beta * (tp_poten15[0] + tp_pV15[0])

    for combo in range(len(combo_dict[prefix])):

        tp_u_kln0[0,(combo + shift) * 5 + 0,:N0] = beta * (tp_poten0[combo + shift] + tp_pV0[combo + shift])
        tp_u_kln0[0,(combo + shift) * 5 + 1,:N0] = beta * (tp_poten0[combo + shift] + tp_pV0[combo + shift] * ((press + delP) / press))
        tp_u_kln0[0,(combo + shift) * 5 + 2,:N0] = beta * (tp_poten0[combo + shift] + tp_pV0[combo + shift] * ((press - delP) / press))
        tp_u_kln0[0,(combo + shift) * 5 + 3,:N0] = (1+delta_beta) * beta * (tp_poten0[combo + shift] + tp_pV0[combo + shift])
        tp_u_kln0[0,(combo + shift) * 5 + 4,:N0] = (1-delta_beta) * beta * (tp_poten0[combo + shift] + tp_pV0[combo + shift])

        tp_u_kln15[0,(combo + shift) * 5 + 0,:N15] = beta * (tp_poten15[combo + shift] + tp_pV15[combo + shift])
        tp_u_kln15[0,(combo + shift) * 5 + 1,:N15] = beta * (tp_poten15[combo + shift] + tp_pV15[combo + shift] * ((press + delP) / press))
        tp_u_kln15[0,(combo + shift) * 5 + 2,:N15] = beta * (tp_poten15[combo + shift] + tp_pV15[combo + shift] * ((press - delP) / press))
        tp_u_kln15[0,(combo + shift) * 5 + 3,:N15] = (1+delta_beta) * beta * (tp_poten15[combo + shift] + tp_pV15[combo + shift])
        tp_u_kln15[0,(combo + shift) * 5 + 4,:N15] = (1-delta_beta) * beta * (tp_poten15[combo + shift] + tp_pV15[combo + shift])    
    return(tp_N_k0, tp_N_k15, tp_u_kln0, tp_u_kln15)

#=========================================================================
# MBAR 
#=========================================================================
def compute_fe(fedict, dfedict, prefix, u_kln, N_k, verbose, relative_tolerance, tp_u_kln0, tp_N_k0, tp_u_kln15, tp_N_k15):
    """
    """
    # Initialize MBAR (computing free energy estimates, which may take a while)
    if verbose:
        print "Computing free energy differences..."
    sub_fedict = dict()
    sub_dfedict = dict()
    MBAR = pymbar.MBAR(u_kln, N_k, verbose = verbose, method = 'adaptive', relative_tolerance = relative_tolerance, initialize = 'BAR') # Use fast Newton-Raphson solver
    # Got matrix of dimensionless free energy differences and uncertainty estimate. 
    if verbose:
        print 'Computing covariance matrix...'

    (Deltaf_ij, dDeltaf_ij, Theta_ij) = MBAR.getFreeEnergyDifferences(return_theta = True)

    MBAR = pymbar.MBAR(tp_u_kln0, tp_N_k0, verbose = verbose, method = 'adaptive', relative_tolerance = relative_tolerance, initialize = 'BAR') # Use fast Newton-Raphson solver
    # Got matrix of dimensionless free energy differences and uncertainty estimate.
    if verbose:
        print 'Computing covariance matrix...'

    (tp_Deltaf_ij0, tp_dDeltaf_ij0, tp_Theta_ij0) = MBAR.getFreeEnergyDifferences(return_theta = True)

    MBAR = pymbar.MBAR(tp_u_kln15, tp_N_k15, verbose = verbose, method = 'adaptive', relative_tolerance = relative_tolerance, initialize = 'BAR') # Use fast Newton-Raphson solver
        # Got matrix of dimensionless free energy differences and uncertainty estimate.
    if verbose:
        print 'Computing covariance matrix...'

    (tp_Deltaf_ij15, tp_dDeltaf_ij15, tp_Theta_ij15) = MBAR.getFreeEnergyDifferences(return_theta = True)
    sub_fedict[0] = Deltaf_ij
    sub_dfedict[0] = dDeltaf_ij
    sub_fedict[1] = tp_Deltaf_ij0
    sub_dfedict[1] = tp_dDeltaf_ij0
    sub_fedict[2] = tp_Deltaf_ij15
    sub_dfedict[2] = tp_dDeltaf_ij15

    fedict[prefix] = sub_fedict
    dfedict[prefix] = sub_dfedict
    return(fedict, dfedict)

#=========================================================================
# DEFINE PROPERTY CALCULATIONS
#=========================================================================
# Compute Density
def compDENS(prefix, combo, shift, beta, molnum, fedict, dfedict, deltaPdict, tmass):
    """
    Compute the average and standard deviation of the density of the system.
    Units are g/cm3
    Pure fluid simulation
    rho = mass/volume = mass/ (dG/dP)
    estimate dG/dP with the central difference formula
    """
    dens = molnum[prefix] * (tmass[prefix]/((1000*fedict[prefix][1][2+(combo+shift)*5,1+(combo+shift)*5]/beta)/(2*deltaPdict[prefix])))*(1./(100.**3)) #g/cm3
    ddens = np.abs(molnum[prefix] * (tmass[prefix]*2*deltaPdict[prefix]*beta/(1000*(100**3)*(fedict[prefix][1][2+(combo+shift)*5,1+(combo+shift)*5]**2)))*dfedict[prefix][1][2+(combo+shift)*5,1+(combo+shift)*5])

    return(dens,ddens)

# Compute Pure Fluid Heat Capacity
def compHCAP(prefix, combo, shift, kB, beta, delta_beta, temp, delT, molnum, fedict, dfedict):
    """
    CHECK DERIVATION
    Compute the Pure Fluid Constant Pressure Heat Capacity of the system.
    Units are J/ (mol K)
    Pure fluid simulation
    """

    delta_beta = beta - 1/(kB * (temp+delT))

    hCap = (-kB * beta**2 * (fedict[prefix][1][0+(combo+shift)*5,3+(combo+shift)*5] + fedict[prefix][1][0+(combo+shift)*5,4+(combo+shift)*5])/(molnum[prefix]*(0.5*delta_beta)**2))*1000
    dhCap = np.sqrt(np.abs(((-1000*beta**2)/(molnum[prefix]*delta_beta**2))*(dfedict[prefix][1][0+(combo+shift)*5,3+(combo+shift)*5]**2 + dfedict[prefix][1][0+(combo+shift)*5,4+(combo+shift)*5]**2)))

    return(hCap, dhCap)

# Compute Pure Fluid Isothermal Compressibility
def compKappa(prefix, combo, shift, fedict, dfedict, deltaPdict):
    """
    Compute the Isothermal Compressibility of the pure fluid (Kappa)
    Units are 1/GPa
    Pure fluid simulation
    """
    kappa = 2*(((fedict[prefix][1][2+(combo+shift)*5,0+(combo+shift)*5] + fedict[prefix][1][1+(combo+shift)*5,0+(combo+shift)*5])/(deltaPdict[prefix]))/(fedict[prefix][1][2+(combo+shift)*5,1+(combo+shift)*5]))*(10**9)

    dkappa = np.abs(np.sqrt((((10**9)/(deltaPdict[prefix]*fedict[prefix][1][2+(combo+shift)*5,1+(combo+shift)*5]))**2)*(dfedict[prefix][1][1+(combo+shift)*5,0+(combo+shift)*5]**2 + dfedict[prefix][1][2+(combo+shift)*5,0+(combo+shift)*5]**2-(((fedict[prefix][1][1+(combo+shift)*5,0+(combo+shift)*5]+fedict[prefix][1][2+(combo+shift)*5,0+(combo+shift)*5])/fedict[prefix][1][2+(combo+shift)*5,1+(combo+shift)*5])**2)*dfedict[prefix][1][2+(combo+shift)*5,1+(combo+shift)*5]**2)))

    return(kappa, dkappa)

# Calculate Infinite Dilution Chemical Potential
def compIDCHEMPOT(prefix, combo, shift, K, beta, fedict, dfedict):
    """
    Compute the Infinite Dilution Chemical Potential of the system.
    Units are kJ/mol
    Single molecule solvation (NPT) simulation (single molecule in a box of waters).
    """
    idChemPot = fedict[prefix][0][K-1+(combo+shift)*K,0+(combo+shift)*K]/beta
    didChemPot = dfedict[prefix][0][K-1+(combo+shift)*K,0+(combo+shift)*K]/beta

    return(idChemPot, didChemPot)

# Calculate Heat of Vaporization
def compHvap(prefix, combo, shift, beta, kB, temp, delT, K, fedict, dfedict):
    """
    Compute the Heat of Vaporization of the system.
    Units are kJ/mol
    """
    delta_beta = beta - 1/(kB * (temp+delT))
    Hvap = fedict[prefix][1][3+(combo+shift)*5,4+(combo+shift)*5]/(2*delta_beta) - fedict[prefix][2][3+(combo+shift)*5,4+(combo+shift)*5]/(2*delta_beta) + fedict[prefix][0][0+(combo+shift)*K,K-1+(combo+shift)*K]/beta
    dHvap = (((1/(2*delta_beta))**2)*((dfedict[prefix][1][3+(combo+shift)*5,4+(combo+shift)*5])**2 + (dfedict[prefix][2][3+(combo+shift)*5,4+(combo+shift)*5])**2)+dfedict[prefix][0][0+(combo+shift)*K,K-1+(combo+shift)*K]/beta)**0.5

    return(Hvap, dHvap)

# Calculate Pure Chemical Potential
def compPCHEMPOT(prefix, combo, shift, K, beta, fedict, dfedict):
    """
    Compute the Pure Chemical Potential of the system.
    Units are kJ/mol
    Single molecule solvation (NPT) simulation (single molecule is a pure fluid box)
    """
    pChemPot = fedict[prefix][0][K-1+(combo+shift)*K,0+(combo+shift)*K]/beta
    dpChemPot = dfedict[prefix][0][K-1+(combo+shift)*K,0+(combo+shift)*K]/beta

    return(pChemPot, dpChemPot)

# Calculate Infinite Dilution Water Chemical Potential
def compIDH2OCHEMPOT(prefix, combo, shift, K, beta, fedict, dfedict):
    """
    Compute the Pure Chemical Potential of the system.
    Units are kJ/mol
    Single molecule solvation (NPT) simulation (single molecule is a pure fluid box)
    """
    idh2oChemPot = fedict[prefix][0][K-1+(combo+shift)*K,0+(combo+shift)*K]/beta
    didh2oChemPot = dfedict[prefix][0][K-1+(combo+shift)*K,0+(combo+shift)*K]/beta

    return(idh2oChemPot, didh2oChemPot)

#=========================================================================
# GENERATE PROPERTIES
#=========================================================================
def gen_properties(prefixes, current_params, initparam, combos, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, beta_pure, beta_t4900ewcoul, beta_idh2o, temp_t4900ewcoul, temp_pure, temp_idh2o, deltaPdict, delT, kB, combo_dict, K, molnum, tmass, fedict, dfedict, convertE):
    """
    """
    # DICTIONARIES TO STORE CALCULATED PROPERTY VALUES IN
    idChemPotdict = dict()
    didChemPotdict = dict()
    idChemPotErrdict = dict()

    idh2oChemPotdict = dict()
    didh2oChemPotdict = dict()
    idh2oChemPotErrdict = dict()

    pChemPotdict = dict()
    dpChemPotdict = dict()
    pChemPotErrdict = dict()

    densdict = dict()
    ddensdict = dict()
    densErrdict = dict()

    hCapdict = dict()
    dhCapdict = dict()
    hCapErrdict = dict()

    kappadict = dict()
    dkappadict = dict()
    kappaErrdict = dict()

    hVapdict = dict()
    dhVapdict = dict()
    hVapErrdict = dict()

    if isinstance(current_params, list):
        if (current_params != initparam):
            shift = 1
        else:
            shift = 0
    else:
        if (current_params != initparam).any():
            shift = 1
        else:
            shift = 0

    returns = list()

    didChemPotdict_initial = dict()
    didh2oChemPotdict_initial = dict()
    dpChemPotdict_initial = dict()
    ddensdict_initial = dict()
    dhCapdict_initial = dict()
    dkappadict_initial = dict()
    dhVapdict_initial = dict()
        
    # Generate properties for all simulations
    for combo in range(combos):
        sub_idChemPotdict = dict()
        sub_didChemPotdict = dict()
        sub_idChemPotErrdict = dict()

        sub_idh2oChemPotdict = dict()
        sub_didh2oChemPotdict = dict()
        sub_idh2oChemPotErrdict = dict()

        sub_pChemPotdict = dict()
        sub_dpChemPotdict = dict()
        sub_pChemPotErrdict = dict()

        sub_densdict = dict()
        sub_ddensdict = dict()
        sub_densErrdict = dict()

        sub_hCapdict = dict()
        sub_dhCapdict = dict()
        sub_hCapErrdict = dict()

        sub_kappadict = dict()
        sub_dkappadict = dict()
        sub_kappaErrdict = dict()

        sub_hVapdict = dict()
        sub_dhVapdict = dict()
        sub_hVapErrdict = dict()

        (EXPdens, EXPdensERR, EXPhCap, EXPhCapERR, EXPkappa, EXPkappaERR, EXPidChemPot, EXPidChemPotERR, EXPhVap, EXPhVapERR, EXPpChemPot, EXPpChemPotERR, EXPidh2oChemPot, EXPidh2oChemPotERR) = experimental_data(convertE, beta_idh2o)

        for prefix in prefixes:
            (dir, dir_inputs, temp, beta, delta_beta, delP) = set_system_variables(prefix, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, beta_pure, beta_t4900ewcoul, beta_idh2o, temp_t4900ewcoul, temp_pure, temp_idh2o, deltaPdict, delT, kB)
            perturbed = False
            counter = -1
            if sim_type[prefix] == 't4900ewcoul':
                (idChemPot, didChemPot) = compIDCHEMPOT(prefix, 0, 0, K, beta, fedict, dfedict)
                didChemPotdict_initial[prefix] = didChemPot
                for number in combo_dict[prefix]:
                    counter += 1
                    if combo == number:
                        perturbed = True
                        comboB = counter

                if not perturbed:
                    comboB = 0
                (idChemPot, didChemPot) = compIDCHEMPOT(prefix, comboB, shift, K, beta, fedict, dfedict)
                idChemPotErr = np.abs((idChemPot - EXPidChemPot[prefix])/EXPidChemPot[prefix])

                sub_idChemPotdict[prefix] = idChemPot
                sub_didChemPotdict[prefix] = didChemPot
                sub_idChemPotErrdict[prefix] = idChemPotErr

            if sim_type[prefix] == 'idh2o':
                (idh2oChemPot, didh2oChemPot) = compIDH2OCHEMPOT(prefix, 0, 0, K, beta, fedict, dfedict)
                didh2oChemPotdict_initial[prefix] = didh2oChemPot
                for number in combo_dict[prefix]:
                    counter += 1
                    if combo == number:
                        perturbed = True
                        comboB = counter

                if not perturbed:
                    comboB = 0
                (idh2oChemPot, didh2oChemPot) = compIDH2OCHEMPOT(prefix, comboB, shift, K, beta, fedict, dfedict)
                idh2oChemPotErr = np.abs((idh2oChemPot - EXPidh2oChemPot[prefix])/EXPidh2oChemPot[prefix])

                sub_idh2oChemPotdict[prefix] = idh2oChemPot
                sub_didh2oChemPotdict[prefix] = didh2oChemPot
                sub_idh2oChemPotErrdict[prefix] = idh2oChemPotErr    

            if sim_type[prefix] == 'pure':
                (dens, ddens) = compDENS(prefix, 0, 0, beta, molnum, fedict, dfedict, deltaPdict, tmass)
                ddensdict_initial[prefix] = ddens
                (hCap, dhCap) = compHCAP(prefix, 0, 0, kB, beta, delta_beta, temp, delT, molnum, fedict, dfedict)
                dhCapdict_initial[prefix] = dhCap
                (kappa, dkappa) = compKappa(prefix, 0, 0, fedict, dfedict, deltaPdict)
                dkappadict_initial[prefix] = dkappa
                (hVap, dhVap) = compHvap(prefix, 0, 0, beta, kB, temp, delT, K, fedict, dfedict)
                dhVapdict_initial[prefix] = dhVap
                (pChemPot, dpChemPot) = compPCHEMPOT(prefix, 0, 0, K, beta, fedict, dfedict)
                dpChemPotdict_initial[prefix] = dpChemPot

                for number in combo_dict[prefix]:
                    counter += 1
                    if combo == number:
                        perturbed = True
                        comboB = counter

                if not perturbed:
                    comboB = 0
                (dens, ddens) = compDENS(prefix, comboB, shift, beta, molnum, fedict, dfedict, deltaPdict, tmass)
                densErr = np.abs((dens - EXPdens[prefix])/EXPdens[prefix])
                sub_densdict[prefix] = dens
                sub_ddensdict[prefix] = ddens
                sub_densErrdict[prefix] = densErr

                (hCap, dhCap) = compHCAP(prefix, comboB, shift, kB, beta, delta_beta, temp, delT, molnum, fedict, dfedict)
                hCapErr = np.abs((hCap - EXPhCap[prefix])/EXPhCap[prefix])
                sub_hCapdict[prefix] = hCap
                sub_dhCapdict[prefix] = dhCap
                sub_hCapErrdict[prefix] = hCapErr

                (kappa, dkappa) = compKappa(prefix, comboB, shift, fedict, dfedict, deltaPdict)
                kappaErr = np.abs((kappa - EXPkappa[prefix])/EXPkappa[prefix])
                sub_kappadict[prefix] = kappa
                sub_dkappadict[prefix] = dkappa
                sub_kappaErrdict[prefix] = kappaErr

                (hVap, dhVap) = compHvap(prefix, comboB, shift, beta, kB, temp, delT, K, fedict, dfedict)
                hVapErr = np.abs((hVap - EXPhVap[prefix])/EXPhVap[prefix])
                sub_hVapdict[prefix] = hVap
                sub_dhVapdict[prefix] = dhVap
                sub_hVapErrdict[prefix] = hVapErr

                (pChemPot, dpChemPot) = compPCHEMPOT(prefix, comboB, shift, K, beta, fedict, dfedict)
                pChemPotErr = np.abs((pChemPot - EXPpChemPot[prefix])/EXPpChemPot[prefix])
                sub_pChemPotdict[prefix] = pChemPot
                sub_dpChemPotdict[prefix] = dpChemPot
                sub_pChemPotErrdict[prefix] = pChemPotErr

        if len(sub_idChemPotdict)>0:
            idChemPotdict[combo] = sub_idChemPotdict
            didChemPotdict[combo] = sub_didChemPotdict
            idChemPotErrdict[combo] = sub_idChemPotErrdict
        if len(sub_idh2oChemPotdict)>0:
            idh2oChemPotdict[combo] = sub_idh2oChemPotdict
            didh2oChemPotdict[combo] = sub_didh2oChemPotdict
            idh2oChemPotErrdict[combo] = sub_idh2oChemPotErrdict
        if len(sub_pChemPotdict)>0:
            pChemPotdict[combo] = sub_pChemPotdict
            dpChemPotdict[combo] = sub_dpChemPotdict
            pChemPotErrdict[combo] = sub_pChemPotErrdict
        if len(sub_densdict)>0:
            densdict[combo] = sub_densdict
            ddensdict[combo] = sub_ddensdict
            densErrdict[combo] = sub_densErrdict
        if len(sub_hCapdict)>0:
            hCapdict[combo] = sub_hCapdict
            dhCapdict[combo] = sub_dhCapdict
            hCapErrdict[combo] = sub_hCapErrdict
        if len(sub_kappadict)>0:
            kappadict[combo] = sub_kappadict
            dkappadict[combo] = sub_dkappadict
            kappaErrdict[combo] = sub_kappaErrdict
        if len(sub_hVapdict)>0:
            hVapdict[combo] = sub_hVapdict
            dhVapdict[combo] = sub_dhVapdict
            hVapErrdict[combo] = sub_hVapErrdict

    return(idChemPotdict, didChemPotdict, idChemPotErrdict, idh2oChemPotdict, didh2oChemPotdict, idh2oChemPotErrdict, pChemPotdict, dpChemPotdict, pChemPotErrdict, densdict, ddensdict, densErrdict, hCapdict, dhCapdict, hCapErrdict, kappadict, dkappadict, kappaErrdict, hVapdict, dhVapdict, hVapErrdict, didChemPotdict_initial, didh2oChemPotdict_initial, dpChemPotdict_initial, ddensdict_initial, dhCapdict_initial, dkappadict_initial, dhVapdict_initial)

#=========================================================================
# CALCULATE OBJECTIVE FUNCTION
#=========================================================================
def calcOBJ(prefixes, sim_type, combo, idChemPotErrdict, idChemPotdict, didChemPotdict_initial, EXPidChemPot, EXPidChemPotERR, idh2oChemPotErrdict, idh2oChemPotdict, didh2oChemPotdict_initial, EXPidh2oChemPot, EXPidh2oChemPotERR, densErrdict, densdict, ddensdict_initial, EXPdens, EXPdensERR, hCapErrdict, hCapdict, dhCapdict_initial, EXPhCap, EXPhCapERR, hVapErrdict, hVapdict, dhVapdict_initial, EXPhVap, EXPhVapERR, kappaErrdict, kappadict, dkappadict_initial, EXPkappa, EXPkappaERR, pChemPotErrdict, pChemPotdict, dpChemPotdict_initial, EXPpChemPot, EXPpChemPotERR, weight):
    errorsum = 0
    for prefix in prefixes:

        if sim_type[prefix] == 't4900ewcoul':
            errorsum = errorsum + (weight['idChemPot']/((didChemPotdict_initial[prefix]/EXPidChemPot[prefix])**2 + (idChemPotdict[combo][prefix]/EXPidChemPot[prefix])**2 * (EXPidChemPotERR[prefix]/EXPidChemPot[prefix])**2))*idChemPotErrdict[combo][prefix]

        if sim_type[prefix] == 'idh2o':
            errorsum = errorsum + (weight['idh2oChemPot']/((didh2oChemPotdict_initial[prefix]/EXPidh2oChemPot[prefix])**2 + (idh2oChemPotdict[combo][prefix]/EXPidh2oChemPot[prefix])**2 * (EXPidh2oChemPotERR[prefix]/EXPidh2oChemPot[prefix])**2))*idh2oChemPotErrdict[combo][prefix]

        if sim_type[prefix] == 'pure':
            errorsum = errorsum + (weight['dens']/((ddensdict_initial[prefix]/EXPdens[prefix])**2 + (densdict[combo][prefix]/EXPdens[prefix])**2 * (EXPdensERR[prefix]/EXPdens[prefix])**2))*densErrdict[combo][prefix]
            errorsum = errorsum + (weight['hCap']/((dhCapdict_initial[prefix]/EXPhCap[prefix])**2 + (hCapdict[combo][prefix]/EXPhCap[prefix])**2 * (EXPhCapERR[prefix]/EXPhCap[prefix])**2))*hCapErrdict[combo][prefix]
            errorsum = errorsum + (weight['hVap']/((dhVapdict_initial[prefix]/EXPhVap[prefix])**2 + (hVapdict[combo][prefix]/EXPhVap[prefix])**2 * (EXPhVapERR[prefix]/EXPhVap[prefix])**2))*hVapErrdict[combo][prefix]
            errorsum = errorsum + (weight['kappa']/((dkappadict_initial[prefix]/EXPkappa[prefix])**2 + (kappadict[combo][prefix]/EXPkappa[prefix])**2 * (EXPkappaERR[prefix]/EXPkappa[prefix])**2))*kappaErrdict[combo][prefix]
            errorsum = errorsum + (weight['pChemPot']/((dpChemPotdict_initial[prefix]/EXPpChemPot[prefix])**2 + (pChemPotdict[combo][prefix]/EXPpChemPot[prefix])**2 * (EXPpChemPotERR[prefix]/EXPpChemPot[prefix])**2))*pChemPotErrdict[combo][prefix]

    return(errorsum)

#=========================================================================
# WRITE LOG FILE
#=========================================================================
def write_log_file(logfile, figurefile, datafile, combos, shift, perturbparam, convertE, idChemPotdict, didChemPotdict, idChemPotErrdict, idh2oChemPotdict, didh2oChemPotdict, idh2oChemPotErrdict, pChemPotdict, dpChemPotdict, pChemPotErrdict, densdict, ddensdict, densErrdict, hCapdict, dhCapdict, hCapErrdict, kappadict, dkappadict, kappaErrdict, hVapdict, dhVapdict, hVapErrdict, objdict, beta_idh2o, errorsum_log, idChemPot_log, didChemPot_log, idh2oChemPot_log, didh2oChemPot_log, pChemPot_log, dpChemPot_log, dens_log, ddens_log, hCap_log, dhCap_log, kappa_log, dkappa_log, hVap_log, dhVap_log, directory):
    """
    """
    (EXPdens, EXPdensERR, EXPhCap, EXPhCapERR, EXPkappa, EXPkappaERR, EXPidChemPot, EXPidChemPotERR, EXPhVap, EXPhVapERR, EXPpChemPot, EXPpChemPotERR, EXPidh2oChemPot, EXPidh2oChemPotERR) = experimental_data(convertE, beta_idh2o)

    # Open log file to write properties.
    log = open(logfile, 'a')
    log.write(datetime.datetime.now().ctime())
    log.write('\n')
    for combo in range(combos):
        pkstring = ''
        pstring = ''
        if combo == 0:
            if shift == 0:
                log.write('PARAMETER SET: 0\n')
            else:
                log.write('NEW PARAMETER SET\n')
        else:
            log.write('PERTURBATION: %(combo)s\n' % vars())
        for value in range(len(paramkey)):
            pkstring += paramkey[value].rjust(13)
            if len(np.shape(perturbparam)) == 1:
                paramvalue = '%.6f'%perturbparam[value]
            else:
                paramvalue = '%.6f'%perturbparam[combo][value]
            pstring += paramvalue.rjust(13)
        pkstring +='\n'
        pstring += '\n'
        log.write(pkstring)
        log.write(pstring)
        errorsum = objdict[combo]
        line = '\nERROR SUM = %(errorsum)s' % vars()
        log.write(line)
        if len(idChemPotdict)>0:
            log.write('\n')
            log.write('Infinite Dilution Chemical Potential (kJ/mol):\n')
            line = 'SYSTEM'.ljust(20) + 'SIM'.rjust(10) + '+/-'.rjust(10) + 'EXP'.rjust(10) + '%ERR'.rjust(10) + '\n'
            log.write(line)
            for prefix in idChemPotdict[combo]:
                a = '%.5f'%idChemPotdict[combo][prefix]
                b = '%.5f'%didChemPotdict[combo][prefix]
                c = '%.5f'%EXPidChemPot[prefix]
                d = '%.5f'%idChemPotErrdict[combo][prefix]
                line = prefix.ljust(20) + a.rjust(10) + b.rjust(10) + c.rjust(10) + d.rjust(10) + '\n'
                log.write(line)

        if len(idh2oChemPotdict)>0:
            log.write('\n')
            log.write('Infinite Dilution Water Chemical Potential (kJ/mol):\n')
            line = 'SYSTEM'.ljust(20) + 'SIM'.rjust(10) + '+/-'.rjust(10) + 'EXP'.rjust(10) + '%ERR'.rjust(10) + '\n'
            log.write(line)
            for prefix in idh2oChemPotdict[combo]:
                a = '%.5f'%idh2oChemPotdict[combo][prefix]
                b = '%.5f'%didh2oChemPotdict[combo][prefix]
                c = '%.5f'%EXPidh2oChemPot[prefix]
                d = '%.5f'%idh2oChemPotErrdict[combo][prefix]
                line = prefix.ljust(20) + a.rjust(10) + b.rjust(10) + c.rjust(10) + d.rjust(10) + '\n'
                log.write(line)

        if len(pChemPotdict)>0:
            log.write('\n')
            log.write('Pure Chemical Potential (kJ/mol):\n')
            line = 'SYSTEM'.ljust(20) + 'SIM'.rjust(10) + '+/-'.rjust(10) + 'EXP'.rjust(10) + '%ERR'.rjust(10) + '\n'
            log.write(line)
            for prefix in pChemPotdict[combo]:
                a = '%.5f'%pChemPotdict[combo][prefix]
                b = '%.5f'%dpChemPotdict[combo][prefix]
                c = '%.5f'%EXPpChemPot[prefix]
                d = '%.5f'%pChemPotErrdict[combo][prefix]
                line = prefix.ljust(20) + a.rjust(10) + b.rjust(10) + c.rjust(10) + d.rjust(10) + '\n'
                log.write(line)

        if len(densdict)>0:
            log.write('\n')
            log.write('Density (g/cm3):\n')
            line = 'SYSTEM'.ljust(20) + 'SIM'.rjust(10) + '+/-'.rjust(10) + 'EXP'.rjust(10) + '%ERR'.rjust(10) + '\n'
            log.write(line)
            for prefix in densdict[combo]:
                a = '%.5f'%densdict[combo][prefix]
                b = '%.5f'%ddensdict[combo][prefix]
                c = '%.5f'%EXPdens[prefix]
                d = '%.5f'%densErrdict[combo][prefix]
                line = prefix.ljust(20) + a.rjust(10) + b.rjust(10) + c.rjust(10) + d.rjust(10) + '\n'
                log.write(line)

        if len(hCapdict)>0:
            log.write('\n')
            log.write('Heat Capacity (J/mol K):\n')
            line = 'SYSTEM'.ljust(20) + 'SIM'.rjust(10) + '+/-'.rjust(10) + 'EXP'.rjust(10) + '%ERR'.rjust(10) + '\n'
            log.write(line)
            for prefix in hCapdict[combo]:
                a = '%.5f'%hCapdict[combo][prefix]
                b = '%.5f'%dhCapdict[combo][prefix]
                c = '%.5f'%EXPhCap[prefix]
                d = '%.5f'%hCapErrdict[combo][prefix]
                line = prefix.ljust(20) + a.rjust(10) + b.rjust(10) + c.rjust(10) + d.rjust(10) + '\n'
                log.write(line)

        if len(hVapdict)>0:
            log.write('\n')
            log.write('Heat of Vaporization (kJ/mol):\n')
            line = 'SYSTEM'.ljust(20) + 'SIM'.rjust(10) + '+/-'.rjust(10) + 'EXP'.rjust(10) + '%ERR'.rjust(10) + '\n'
            log.write(line)
            for prefix in hVapdict[combo]:
                a = '%.5f'%hVapdict[combo][prefix]
                b = '%.5f'%dhVapdict[combo][prefix]
                c = '%.5f'%EXPhVap[prefix]
                d = '%.5f'%hVapErrdict[combo][prefix]
                line = prefix.ljust(20) + a.rjust(10) + b.rjust(10) + c.rjust(10) + d.rjust(10) + '\n'
                log.write(line)

        if len(kappadict)>0:
            log.write('\n')
            log.write('Isothermal Compressibility (1/GPa):\n')
            line = 'SYSTEM'.ljust(20) + 'SIM'.rjust(10) + '+/-'.rjust(10) + 'EXP'.rjust(10) + '%ERR'.rjust(10) + '\n'
            log.write(line)
            for prefix in kappadict[combo]:
                a = '%.5f'%kappadict[combo][prefix]
                b = '%.5f'%dkappadict[combo][prefix]
                c = '%.5f'%EXPkappa[prefix]
                d = '%.5f'%kappaErrdict[combo][prefix]
                line = prefix.ljust(20) + a.rjust(10) + b.rjust(10) + c.rjust(10) + d.rjust(10) + '\n'
                log.write(line)
        log.write('\n\n')
    log.close()

    # Open datafile to write data for figures
    data = open(datafile, 'w')
    data.write(datetime.datetime.now().ctime())
    data.write('\n')
    data.write('\n')

    # Only record data for combo = 0 in these files (data from initial paramters and optimized parameters, not perturbed parameters) 
    combo = 0
    #Append new data:
    errorsum_log.append(objdict[combo])
    data.write('Objective Function\n')
    for item in errorsum_log:
        data.write('%s        ' % item)
    data.write('\n')
    
    if len(idChemPotdict)>0:
        data.write('\n')
        data.write('Infinite Dilution Chemical Potentials (kJ/mol)\n')
        for prefix in idChemPotdict[combo]:
            expvalue = EXPidChemPot[prefix]
            experror = EXPidChemPotERR[prefix]
            line = '%(prefix)s     experimental value: %(expvalue)s +- %(experror)s'% vars()
            data.write(line)
            data.write('\n')
            idChemPot_log[prefix].append(idChemPotdict[combo][prefix])
            for item in idChemPot_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')
            didChemPot_log[prefix].append(didChemPotdict[combo][prefix])
            for item in didChemPot_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')

    if len(idh2oChemPotdict)>0:
        data.write('\n')
        data.write('Inifinite Dilution Water Chemical Potentials (kJ/mol)\n')
        for prefix in idh2oChemPotdict[combo]:
            expvalue = EXPidh2oChemPot[prefix]
            experror = EXPidh2oChemPotERR[prefix]
            line = '%(prefix)s     experimental value: %(expvalue)s +- %(experror)s'% vars()
            data.write(line)
            data.write('\n')
            idh2oChemPot_log[prefix].append(idh2oChemPotdict[combo][prefix])
            for item in idh2oChemPot_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')
            didh2oChemPot_log[prefix].append(didh2oChemPotdict[combo][prefix])
            for item in didh2oChemPot_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')

    if len(pChemPotdict)>0:
        data.write('\n')
        data.write('Pure Chemical Potentials (kJ/mol)\n')
        for prefix in pChemPotdict[combo]:
            expvalue = EXPpChemPot[prefix]
            experror = EXPpChemPotERR[prefix]
            line = '%(prefix)s     experimental value: %(expvalue)s +- %(experror)s'% vars()
            data.write(line)
            data.write('\n')
            pChemPot_log[prefix].append(pChemPotdict[combo][prefix])
            for item in pChemPot_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')
            dpChemPot_log[prefix].append(dpChemPotdict[combo][prefix])
            for item in dpChemPot_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')

    if len(densdict)>0:
        data.write('\n')
        data.write('Densities (g/cm3)\n')
        for prefix in densdict[combo]:
            expvalue = EXPdens[prefix]
            experror = EXPdensERR[prefix]
            line = '%(prefix)s     experimental value: %(expvalue)s +- %(experror)s'% vars()
            data.write(line)
            data.write('\n')
            dens_log[prefix].append(densdict[combo][prefix])
            for item in dens_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')
            ddens_log[prefix].append(ddensdict[combo][prefix])
            for item in ddens_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')

    if len(hCapdict)>0:
        data.write('\n')
        data.write('Heat Capacitites(J/mol K)\n')
        for prefix in hCapdict[combo]:
            expvalue = EXPhCap[prefix]
            experror = EXPhCapERR[prefix]
            line = '%(prefix)s     experimental value: %(expvalue)s +- %(experror)s'% vars()
            data.write(line)
            data.write('\n')
            hCap_log[prefix].append(hCapdict[combo][prefix])
            for item in hCap_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')
            dhCap_log[prefix].append(dhCapdict[combo][prefix])
            for item in dhCap_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')

    if len(hVapdict)>0:
        data.write('\n')
        data.write('Enthalpy of Vaporization (kJ/mol)\n')
        for prefix in hVapdict[combo]:
            expvalue = EXPhVap[prefix]
            experror = EXPhVapERR[prefix]
            line = '%(prefix)s     experimental value: %(expvalue)s +- %(experror)s'% vars()
            data.write(line)
            data.write('\n')
            hVap_log[prefix].append(hVapdict[combo][prefix])
            for item in hVap_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')
            dhVap_log[prefix].append(dhVapdict[combo][prefix])
            for item in dhVap_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')

    if len(kappadict)>0:
        data.write('\n')
        data.write('Isothermal Compressibilities (1/GPa)\n')
        for prefix in kappadict[combo]:
            expvalue = EXPkappa[prefix]
            experror = EXPkappaERR[prefix]
            line = '%(prefix)s     experimental value: %(expvalue)s +- %(experror)s'% vars()
            data.write(line)
            data.write('\n')
            kappa_log[prefix].append(kappadict[combo][prefix])
            for item in kappa_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')
            dkappa_log[prefix].append(dkappadict[combo][prefix])
            for item in dkappa_log[prefix]:
                data.write('%s        ' % item)
            data.write('\n')
    """
    figures = open(figurefile, 'w')
    figures.write(datetime.datetime.now().ctime())
    figures.write('\n')
    
    copy_command = 'cp %(logfile)s %(directory)s/%(logfile)s' % vars()
    subprocess.call(copy_command, shell = True) 
    
    copy_command = 'cp %(datafile)s %(directory)s/%(datafile)s' % vars()
    subprocess.call(copy_command, shell = True)
 
    copy_command = 'cp %(figurefile)s %(directory)s/%(figurefile)s' % vars()
    subprocess.call(copy_command, shell = True) 
    """
    return()

#=========================================================================
# MAIN: GENERATE SINGLE POINT CALCULATION OF OBJECTIVE FUNCTION
#=========================================================================
def main_single_point(current_params):

    """
    """
    #pdb.set_trace()
    print datetime.datetime.now().ctime()
    # PARSE INPUTS
    (cluster, optimize, verbose, temp_t4900ewcoul, temp_pure, temp_idh2o, delT, press, nequil, logfile, figurefile, datafile, pbsfilename, single, gro_loc, job_script, orestart, directory) = parse_inputs()
    # SET PARAMETERS
    (prefixes, nsamples, cwd, dir_idh2o, dir_t4900ewcoul, dir_pure, dir_inputs_idh2o, dir_inputs_t4900ewcoul, dir_inputs_pure, dir_tops, g_energy, grompp, mdrun, ke_call, pe_call, v_call, energy_suffix, weight, pbs_lines, maxcount, errortol, error, rr_states, not_rr_states, n_systems) = set_parameters(single, gro_loc)
    (kB, NA, relative_tolerance, convertP, convertE, mol_mass, beta_t4900ewcoul, beta_pure, beta_idh2o) = set_constants()
    # EXPERIMENTAL DATA
    (EXPdens, EXPdensERR, EXPhCap, EXPhCapERR, EXPkappa, EXPkappaERR, EXPidChemPot, EXPidChemPotERR, EXPhVap, EXPhVapERR, EXPpChemPot, EXPpChemPotERR, EXPidh2oChemPot, EXPidh2oChemPotERR) = experimental_data(convertE, beta_idh2o)
    # INITIAL PARAMETERS
    (paramkey, initparam, numparam) = initial_parameters()
    # READ TOPOLOGIES
    (tmass, sim_type, molnum) = read_top(prefixes, dir_pure, dir_t4900ewcoul, dir_idh2o, mol_mass)
    # GENERATE DELTA P
    (deltaPdict, max_time) = gen_delta_p(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, v_call, g_energy, EXPkappa, verbose)
    # GENERATE META DATA
    (nsnapshots, n_files, n_components, n_states, maxn, bPV, bExpanded, bEnergy, lv) = gen_meta_data(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, energy_suffix, nsamples, verbose)
    # GENERATE DICTIONARIES OF MOLECULE TYPES INCLUDED IN EACH SYSTEM
    (mol_type_dict) = gen_mol_types(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o)

    combos = 1        
    if (current_params == initparam).all():
        shift = 0
    else:
        shift = 1

    # GENERATE U_KLT FOR CURRENT PARAMETERS
    if shift == 0:
        pass
    else:
        # RERUN AT OPTIMIZED PARAMETERS
        # EDIT TOPOLOGIES
        perturbparam = np.ones((1,len(current_params)))*current_params
        inner_loop = True
        edit_top(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, cwd, paramkey, initparam, perturbparam, dir_tops, combos, inner_loop)
        #RERUN AT OPTIMIZED PARAMETERS
        rerun_opt_param(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, paramkey, n_states, cwd, grompp, mdrun, cluster, optimize, pbs_lines, verbose, max_time, pbsfilename)
    
    # INITIALIZE DICTIONARIES FOR STORING FREE ENERGIES
    fedict = dict()
    dfedict = dict()
    combo_dict = dict()

    for prefix in prefixes:
        combo_dict[prefix] = [0]
        (dir, dir_inputs, temp, beta, delta_beta, delP) = set_system_variables(prefix, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, beta_pure, beta_t4900ewcoul, beta_idh2o, temp_t4900ewcoul, temp_pure, temp_idh2o, deltaPdict, delT, kB)

        # INITIALIZE MATRICES FOR STORING DATA
        K = n_states
        length = K * (1 + shift)
        lv = lv[0:K,0:n_components] # lambda values
        dhdlt = np.zeros([K,n_components, maxn], float) # dhdlt[k,n,t] is the derivative of energy component n with respect to state k of snapshot t
        u_klt = np.zeros([length, length, maxn], np.float64) # u_klt[k,m,t] is the reduced potential energy of snapshot t of state k evaluated at state m
        tp_pe0 = np.zeros((len(combo_dict[prefix]) + shift, maxn))
        tp_pv0 = np.zeros((len(combo_dict[prefix]) + shift, maxn))
        tp_pe15 = np.zeros((len(combo_dict[prefix]) + shift, maxn))
        tp_pv15 = np.zeros((len(combo_dict[prefix]) + shift, maxn))
        nsnapshots = np.zeros(length, int) #nsnapshots[k] is the number of states from file k

        if bExpanded:
            u_t = np.zeros([maxn], np.float64) # u_t[k,m,t] is the reduced potential energy as a function of time

        # GENERATE U_KLT FOR INITIAL PARAMETERS
        (u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15) = gen_u_klt_initial(prefix, n_states, dir, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, K, max_time, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, nsamples, verbose)

        if shift == 1:
            # ADD TO U_KLT
            (u_klt, tp_pe0, tp_pv0, tp_pe15, tp_pv15) = add_u_klt_optimized(prefix, n_states, initparam, current_params, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, max_time, K, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, cwd, nsamples, verbose)

        # SUBSAMPLE DATA
        (u_kln, tp_poten0, tp_poten15, tp_pV0, tp_pV15, N_k, N0, N15) = subsample(prefix, nsnapshots, u_klt, shift, combo_dict, K, n_components, bExpanded, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, verbose)
        z = np.shape(u_kln)[2]

        # REWEIGHT AT +/- T AND +/- P
        (tp_N_k0, tp_N_k15, tp_u_kln0, tp_u_kln15) = reweight_T_P(current_params, initparam, z, prefix, N0, N15, beta, delta_beta, press, delP, tp_poten0, tp_pV0, tp_poten15, tp_pV15, combo_dict)

        # CALCULATE WITH MBAR
        (fedict, dfedict) = compute_fe(fedict, dfedict, prefix, u_kln, N_k, verbose, relative_tolerance, tp_u_kln0, tp_N_k0, tp_u_kln15, tp_N_k15)
    # CALCULATE AND STORE PROPERTIES FOR ALL SYSTEMS
    (idChemPotdict, didChemPotdict, idChemPotErrdict, idh2oChemPotdict, didh2oChemPotdict, idh2oChemPotErrdict, pChemPotdict, dpChemPotdict, pChemPotErrdict, densdict, ddensdict, densErrdict, hCapdict, dhCapdict, hCapErrdict, kappadict, dkappadict, kappaErrdict, hVapdict, dhVapdict, hVapErrdict, didChemPotdict_initial, didh2oChemPotdict_initial, dpChemPotdict_initial, ddensdict_initial, dhCapdict_initial, dkappadict_initial, dhVapdict_initial)= gen_properties(prefixes, current_params, initparam, combos, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, beta_pure, beta_t4900ewcoul, beta_idh2o, temp_t4900ewcoul, temp_pure, temp_idh2o, deltaPdict, delT, kB, combo_dict, K, molnum, tmass, fedict, dfedict, convertE)

    # CALCULATE OBJECTIVE FUNCTION
    objdict = dict()
    for combo in range(combos):
        objdict[combo] = calcOBJ(prefixes, sim_type, combo, idChemPotErrdict, idChemPotdict, didChemPotdict_initial, EXPidChemPot, EXPidChemPotERR, idh2oChemPotErrdict, idh2oChemPotdict, didh2oChemPotdict_initial, EXPidh2oChemPot, EXPidh2oChemPotERR, densErrdict, densdict, ddensdict_initial, EXPdens, EXPdensERR, hCapErrdict, hCapdict, dhCapdict_initial, EXPhCap, EXPhCapERR, hVapErrdict, hVapdict, dhVapdict_initial, EXPhVap, EXPhVapERR, kappaErrdict, kappadict, dkappadict_initial, EXPkappa, EXPkappaERR, pChemPotErrdict, pChemPotdict, dpChemPotdict_initial, EXPpChemPot, EXPpChemPotERR, weight)
    # WRITE LOG FILE
    if optimize:
        perturbparam = current_params
        write_log_file(logfile, figurefile, datafile, combos, shift, perturbparam, convertE, idChemPotdict, didChemPotdict, idChemPotErrdict, idh2oChemPotdict, didh2oChemPotdict, idh2oChemPotErrdict, pChemPotdict, dpChemPotdict, pChemPotErrdict, densdict, ddensdict, densErrdict, hCapdict, dhCapdict, hCapErrdict, kappadict, dkappadict, kappaErrdict, hVapdict, dhVapdict, hVapErrdict, objdict, beta_idh2o, errorsum_log, idChemPot_log, didChemPot_log, idh2oChemPot_log, didh2oChemPot_log, pChemPot_log, dpChemPot_log, dens_log, ddens_log, hCap_log, dhCap_log, kappa_log, dkappa_log, hVap_log, dhVap_log, directory)

    else:
        perturbparam = current_params
        write_log_file(logfile, figurefile, datafile, combos, shift, perturbparam, convertE, idChemPotdict, didChemPotdict, idChemPotErrdict, idh2oChemPotdict, didh2oChemPotdict, idh2oChemPotErrdict, pChemPotdict, dpChemPotdict, pChemPotErrdict, densdict, ddensdict, densErrdict, hCapdict, dhCapdict, hCapErrdict, kappadict, dkappadict, kappaErrdict, hVapdict, dhVapdict, hVapErrdict, objdict, beta_idh2o, errorsum_log, idChemPot_log, didChemPot_log, idh2oChemPot_log, didh2oChemPot_log, pChemPot_log, dpChemPot_log, dens_log, ddens_log, hCap_log, dhCap_log, kappa_log, dkappa_log, hVap_log, dhVap_log, directory)


    return(objdict[0])
    
        
#=========================================================================
# MAIN: GENERATE DERIVATIVES OF OBJECTIVE FUNCTION
#=========================================================================
def main_derivatives(current_params):

    """
    """
    #pdb.set_trace()
    # PARSE INPUTS
    (cluster, optimize, verbose, temp_t4900ewcoul, temp_pure, temp_idh2o, delT, press, nequil, logfile, figurefile, datafile, pbsfilename, single, gro_loc, job_script, orestart, directory) = parse_inputs()
    # SET PARAMETERS
    (prefixes, nsamples, cwd, dir_idh2o, dir_t4900ewcoul, dir_pure, dir_inputs_idh2o, dir_inputs_t4900ewcoul, dir_inputs_pure, dir_tops, g_energy, grompp, mdrun, ke_call, pe_call, v_call, energy_suffix, weight, pbs_lines, maxcount, errortol, error, rr_states, not_rr_states, n_systems) = set_parameters(single, gro_loc)
    (kB, NA, relative_tolerance, convertP, convertE, mol_mass, beta_t4900ewcoul, beta_pure, beta_idh2o)= set_constants()
    # EXPERIMENTAL DATA
    (EXPdens, EXPdensERR, EXPhCap, EXPhCapERR, EXPkappa, EXPkappaERR, EXPidChemPot, EXPidChemPotERR, EXPhVap, EXPhVapERR, EXPpChemPot, EXPpChemPotERR, EXPidh2oChemPot, EXPidh2oChemPotERR) = experimental_data(convertE, beta_idh2o)
    # INITIAL PARAMETERS
    (paramkey, initparam, numparam) = initial_parameters()
    # READ TOPOLOGIES
    (tmass, sim_type, molnum) = read_top(prefixes, dir_pure, dir_t4900ewcoul, dir_idh2o, mol_mass)
    # GENERATE DELTA P
    (deltaPdict, max_time) = gen_delta_p(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, v_call, g_energy, EXPkappa, verbose)
    # GENERATE META DATA
    (nsnapshots, n_files, n_components, n_states, maxn, bPV, bExpanded, bEnergy, lv) = gen_meta_data(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, energy_suffix, nsamples, verbose)
    # GENERATE DICTIONARIES OF MOLECULE TYPES INCLUDED IN EACH SYSTEM
    (mol_type_dict) = gen_mol_types(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o)
    # PERTURB PARAMETERS
    (initperturb, perturbparam, combos) = perturb_params(optimize, paramkey, initparam, numparam, current_params)
    # EDIT TOPOLOGIES
    inner_loop = True
    edit_top(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, cwd, paramkey, initparam, perturbparam, dir_tops, combos, inner_loop)
    # RERUN AT PERTURBED PARAMETERS
    (combo_dict) = rerun_perturbed(prefixes, sim_type, cwd, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, combos, paramkey, perturbparam, initparam, mol_type_dict, cluster, optimize, pbs_lines, rr_states, mdrun, grompp, max_time, nsamples, verbose, pbsfilename)

    # CALCULATE SHIFT
    if (current_params == initparam).all():
        shift = 0
    else:
        shift = 1

    # INITIALIZE DICTIONARIES FOR STORING FREE ENERGIES
    fedict = dict()
    dfedict = dict()
    for prefix in prefixes:
        # Set variables that depend on simulation type
        if sim_type[prefix] == 'pure':
            dir = dir_pure + prefix.split('_')[0] + '/'
            dir_inputs = dir_inputs_pure + prefix.split('_')[0] + '/'
            temp = temp_pure
            beta = beta_pure
        elif sim_type[prefix] == 't4900ewcoul':
            ### This dir command is not very resilient, find a way to generalize it
            dir = dir_t4900ewcoul + prefix.split('_')[0][:3] + '/'
            dir_inputs = dir_inputs_t4900ewcoul + prefix.split('_')[0][:3] + '/'
            temp = temp_t4900ewcoul
            beta = beta_t4900ewcoul
        elif sim_type[prefix] == 'idh2o':
            ### This dir command is not very resilient, find a way to generalize it
            dir = dir_idh2o + prefix.split('_')[0] + '/'
            dir_inputs = dir_inputs_idh2o + prefix.split('_')[0] + '/'
            temp = temp_idh2o
            beta = beta_idh2o
        delta_beta = beta - 1/(kB * (temp+delT))
        delP = deltaPdict[prefix]

        # INITIALIZE MATRICES FOR STORING DATA
        K = n_states
        length = K * (len(combo_dict[prefix]) + shift)
        lv = lv[0:K, 0:n_components]
        dhdlt = np.zeros([K, n_components, maxn], float) # dhdlt[k,n,t] is the derivative of energy component n with respect to state k of snapshot t
        u_klt = np.zeros([length, length,maxn], np.float64) # u_klt[k,m,t] is the reduced potential energy of snapshot t of state k evaluated at state m
        tp_pe0 = np.zeros((len(combo_dict[prefix]) + shift, maxn))
        tp_pv0 = np.zeros((len(combo_dict[prefix]) + shift, maxn))
        tp_pe15 = np.zeros((len(combo_dict[prefix]) + shift, maxn))
        tp_pv15 = np.zeros((len(combo_dict[prefix]) + shift, maxn))
        nsnapshots = np.zeros(length, int) #nsnapshots[k] is the number of states from file k

        if bExpanded:
            u_t = np.zeros([maxn], np.float64) # u_t[k,m,t] is the reduced potential energy as a function of time

        # GENERATE U_KLT FOR INITIAL PARAMETERS
        (u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15) = gen_u_klt_initial(prefix, n_states, dir, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, K, max_time, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, nsamples, verbose)

        # GENERATE U_KLT FOR CURRENT PARAMETERS
        if shift == 0:
            pass
        else:
            (u_klt, tp_pe0, tp_pv0, tp_pe15, tp_pv15) = add_u_klt_optimized(prefix, n_states, initparam, current_params, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, max_time, K, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, cwd, nsamples, verbose)

        # GENERATE U_KLT FOR PERTURBED PARAMETERS
        (u_klt, tp_pe0, tp_pv0, tp_pe15, tp_pv15) = add_u_klt_perturbed(prefix, combo_dict, rr_states, n_states, initparam, current_params, u_klt, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, K, max_time, n_components, maxn, g_energy, ke_call, nsnapshots, bEnergy, bExpanded, bPV, beta, cwd, nsamples, verbose)

        # SUBSAMPLE DATA
        (u_kln, tp_poten0, tp_poten15, tp_pV0, tp_pV15, N_k, N0, N15) = subsample(prefix, nsnapshots, u_klt, shift, combo_dict, K, n_components, bExpanded, dhdlt, tp_pe0, tp_pv0, tp_pe15, tp_pv15, verbose)
        z = np.shape(u_kln)[2]
 
        # REWEIGHT AT +/- T AND +/- P
        (tp_N_k0, tp_N_k15, tp_u_kln0, tp_u_kln15) = reweight_T_P(current_params, initparam, z, prefix, N0, N15, beta, delta_beta, press, delP, tp_poten0, tp_pV0, tp_poten15, tp_pV15, combo_dict)

        # CALCULATE WITH MBAR
        (fedict, dfedict) = compute_fe(fedict, dfedict, prefix, u_kln, N_k, verbose, relative_tolerance, tp_u_kln0, tp_N_k0, tp_u_kln15, tp_N_k15)
    # CALCULATE AND STORE PROPERTIES FOR ALL SYSTEMS
    (idChemPotdict, didChemPotdict, idChemPotErrdict, idh2oChemPotdict, didh2oChemPotdict, idh2oChemPotErrdict, pChemPotdict, dpChemPotdict, pChemPotErrdict, densdict, ddensdict, densErrdict, hCapdict, dhCapdict, hCapErrdict, kappadict, dkappadict, kappaErrdict, hVapdict, dhVapdict, hVapErrdict, didChemPotdict_initial, didh2oChemPotdict_initial, dpChemPotdict_initial, ddensdict_initial, dhCapdict_initial, dkappadict_initial, dhVapdict_initial)= gen_properties(prefixes, current_params, initparam, combos, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, beta_pure, beta_t4900ewcoul, beta_idh2o, temp_t4900ewcoul, temp_pure, temp_idh2o, deltaPdict, delT, kB, combo_dict, K, molnum, tmass, fedict, dfedict, convertE)

    # CALCULATE OBJECTIVE FUNCTION
    objdict = dict()
    for combo in range(combos):
        objdict[combo]= calcOBJ(prefixes, sim_type, combo, idChemPotErrdict, idChemPotdict, didChemPotdict_initial, EXPidChemPot, EXPidChemPotERR, idh2oChemPotErrdict, idh2oChemPotdict, didh2oChemPotdict_initial, EXPidh2oChemPot, EXPidh2oChemPotERR, densErrdict, densdict, ddensdict_initial, EXPdens, EXPdensERR, hCapErrdict, hCapdict, dhCapdict_initial, EXPhCap, EXPhCapERR, hVapErrdict, hVapdict, dhVapdict_initial, EXPhVap, EXPhVapERR, kappaErrdict, kappadict, dkappadict_initial, EXPkappa, EXPkappaERR, pChemPotErrdict, pChemPotdict, dpChemPotdict_initial, EXPpChemPot, EXPpChemPotERR, weight)


    # WRITE LOG FILE
    #write_log_file(logfile, figurefile, datafile, combos, shift, perturbparam, convertE, idChemPotdict, didChemPotdict, idChemPotErrdict, idh2oChemPotdict, didh2oChemPotdict, idh2oChemPotErrdict, pChemPotdict, dpChemPotdict, pChemPotErrdict, densdict, ddensdict, densErrdict, hCapdict, dhCapdict, hCapErrdict, kappadict, dkappadict, kappaErrdict, hVapdict, dhVapdict, hVapErrdict, objdict, beta_idh2o, errorsum_log, idChemPot_log, didChemPot_log, idh2oChemPot_log, didh2oChemPot_log, pChemPot_log, dpChemPot_log, dens_log, ddens_log, hCap_log, dhCap_log, kappa_log, dkappa_log, hVap_log, dhVap_log, directory)

    # CALCULATE DERIVATIVES
    derivative = np.zeros(numparam, np.float64)
    for num in range(numparam):
        derivative[num] = (objdict[num + 1] - objdict[num + 1 + numparam])/(2*initperturb[num])

    derivative = derivative
    return(derivative)

#=========================================================================
# RESIMULATE AT OPTIMIZED PARAMETERS
#=========================================================================
def resim_opt(prefixes, dir_pure, dir_inputs_pure, dir_idh2o, dir_inputs_idh2o, dir_t4900ewcoul, dir_inputs_t4900ewcoul, cwd, n_states, gro_loc, verbose, job_script, cluster, optimize, sim_type):
    qsub_command = 'qsub'
    resim_jobs = 0
    for prefix in prefixes:
        # Set variables that depend on simulation type
        if sim_type[prefix] == 'pure':
            dir = dir_pure + prefix.split('_')[0] + '/'
            dir_inputs = dir_inputs_pure + prefix.split('_')[0] + '/'
        if sim_type[prefix] == 'idh2o':
            dir = dir_idh2o + prefix.split('_')[0] + '/'
            dir_inputs = dir_inputs_idh2o + prefix.split('_')[0] + '/'
        elif sim_type[prefix] == 't4900ewcoul':
            ### This dir command is not very resilient, find a way to make it more generic
            dir = dir_t4900ewcoul + prefix.split('_')[0][:3] + '/'
            dir_inputs = dir_inputs_t4900ewcoul + prefix.split('_')[0][:3] + '/'

        origitp = dir_inputs + 'gaforcefield.itp'
        itpin = open(origitp, 'r')
        lines = itpin.readlines()
        itpin.close()
        newitp = cwd + 'gaforcefield.itp'
        itpout = open(newitp, 'w')
        for line in lines:
            if len(line.split()) > 0 and line.split()[0] == '#include':
                newline = '#include "%(cwd)sgaffnonbonded.o.itp"'
                itpout.write(newline)
            else:
                itpout.write(line)
            itpout.close()
            subprocess.call('mv %(newitp)s %(origitp)s' % vars(), shell = True)

            # Run grompp over all states
            for state in range(n_states):
                os.system('%(gro_loc)s/grompp_d -f %(dir)s%(prefix)s.%(state)s.mdp -c %(dir)s%(prefix)s.gro -p %(dir)s%(prefix)s.top -o %(prefix)s.%(state)s.tpr -maxwarn 2' % locals())
            if verbose:
                print 'finished running grompp the %(n_states)s of %(prefix)s' % vars()

            #Create dictionary with parameter values for modifying generic jobfile
            parameterdict = {'LOCATION' : gro_loc, 'NAME' : prefix}
            
            for state in range(n_states):
                #Create shell script to submit PBS job (single simulation at lambda=state)
                #Add parameter dictionary member for current lambda state
                parameterdict['NUM'] = state
                
                #Open generic job file for reading
                genericjobfile = open(job_script, 'r')
                
                #Create a new job file for reading
                jobfilename = 'pbsjob%s.sh' % (state)
                jobfile = open(jobfilename, 'w')
                
                #Iterate through the lines of the generic job file
                for line in genericjobfile:
                    #Substitute keys for values in the template line and output result
                    jobfile.write(string.Template(line).substitute(parameterdict))
                genericjobfile.close()
                jobfile.close()

            #Call shell script to sumbit PBS job (single simulation at lambda = x)

            if cluster:
                os.system('%s %s' % (qsub_command, jobfilename))
            else:
                #st = os.state(jobfilename)
                #os.chmod(jobfilename, st.st_mode \ stat.S_IXUSR)
                #jobfilename = './' + jobfilename
                #os.system(jobfilename)
                pass
            # Sleep a little while to make sure that PBS has time to copy the script, and doesn't overload qsub
            os.system('sleep 1')
            resim_jobs += 1
    #IF RUNNING ON THE CLUSTER, CHECK IF THE FILES ARE FINISHED RUNNING!
    if cluster and optimize:
        finished = getoutput('qstat -u bsz9ur').split()
        ### CHANGE so that user name is an input (bsz9ur is an input and can be easily changed)
        while resim_jobs > 0:
            resim_jobs = 0
            jobs = len(finished[27:])/11
            resim_jobs = 0
            for job in range(jobs):
                status = finished[36+job*11]
                if status == 'S':
                    pdb.set_trace()
                    kill_command ='qdel ' + finished[16+job*11][:-4]
                    os.system(kill_command)
                if finished[30+job*11] == jobfilename[0:10]:
                    resim_jobs += 1
            time.sleep(600)
            finished = getoutput('qstat -u bsz9ur').split()
        #rm_command = 'rm rerun*'
        #os.system(rm_command)
    
    #Check to make sure that all files were completely written (bigtmp has an intermittant read write error)
    ### WRITE THIS PART LATER

    return()

#=========================================================================
#MAIN: RESIMULATE AT OPTIMIZED PARAMETERS
#=========================================================================
def main_resim(prefixes, dir_pure, dir_inputs_pure, dir_t4900ewcoul, dir_inputs_t4900ewcoul, dir_idh2o, dir_inputs_idh2o, cwd, paramkey, initparam, current_params, dir_tops, gro_loc, verbose, job_script, cluster, optimize):
    pdb.set_trace()
    sim_type = dict()
    for prefix in prefixes:
        type = prefix.split('_')[-1]
        sim_type[prefix] = type
    perturbparam = current_params
    combos = 1
    inner_loop = False
    n_states = 16 # fix this eventually
    edit_top(prefixes, sim_type, dir_pure, dir_t4900ewcoul, dir_idh2o, cwd, paramkey, initparam, perturbparam, dir_tops, combos, inner_loop)
    resim_opt(prefixes, dir_pure, dir_inputs_pure, dir_idh2o, dir_inputs_idh2o, dir_t4900ewcoul, dir_inputs_t4900ewcoul, cwd, n_states, gro_loc, verbose, job_script, cluster, optimize, sim_type)
    return()

#=========================================================================
# RESTART
#=========================================================================
def restart_optimization(logfile, datafile, initparam, prefixes):
    #Find the parameters that yield the lowest objective function
    lines_in = open(logfile, 'r')
    lines = lines_in.readlines()
    lines_in.close()
    counter = 0
    current_params = initparam
    for line in lines:
        elements = line.split()
        if len(elements)>3 and elements[0] == 'ERROR':
            subcounter = 0
            parameters = lines[counter-2].split()
            for parameter in parameters:
                current_params[0,subcounter] = float(parameter)
                subcounter += 1
        counter += 1
    log = open(logfile, 'a')
    log.write('\n RESTARTING OPTIMIZATION\n')
    log.close()

    # Get previous optimization data from datafile
    data_in = open(datafile, 'r')
    lines = data_in.readlines()
    data_in.close()
    counter = 0

    record_idChemPot = False
    record_pChemPot = False
    record_dens = False
    record_hCap = False
    record_hVap = False
    record_kappa = False
    record_idh2oChemPot = False

    for line in lines:
        elements = line.split()
        if line == 'Objective Function\n':
            values = lines[counter+1].split()
            for value in values:
                errorsum_log.append(value)
        if line == 'Infinite Dilution Chemical Potentials (kJ/mol)\n':
            record_idChemPot = True
            record_pChemPot = False
            record_dens = False
            record_hCap = False
            record_hVap = False
            record_kappa = False
            record_idh2oChemPot = False
        if line == 'Pure Chemical Potentials (kJ/mol)\n':
            record_idChemPot = False
            record_pChemPot = True
            record_dens = False
            record_hCap = False
            record_hVap = False
            record_kappa = False
            record_idh2oChemPot = False
        if line == 'Densities (g/cm3)\n':
            record_idChemPot = False
            record_pChemPot = False
            record_dens = True
            record_hCap = False
            record_hVap = False
            record_kappa = False
            record_idh2oChemPot = False
        if line == 'Heat Capacities (J/mol K)\n':
            record_idChemPot = False
            record_pChemPot = False
            record_dens = False
            record_hCap = True
            record_hVap = False
            record_kappa = False
            record_idh2oChemPot = False
        if line == 'Enthalpy of Vaporization (kJ/mol)\n':
            record_idChemPot = False
            record_pChemPot = False
            record_dens = False
            record_hCap = False
            record_hVap = True
            record_kappa = False
            record_idh2oChemPot = False
        if line == 'Isothermal Compressibilities (1/GPa)\n':
            record_idChemPot = False
            record_pChemPot = False
            record_dens = False
            record_hCap = False
            record_hVap = False
            record_kappa = True
            record_idh2oChemPot = False
        if line == 'Inifinite Dilution Water Chemical Potentials (kJ/mol)\n':
            record_idChemPot = False
            record_pChemPot = False
            record_dens = False
            record_hCap = False
            record_hVap = False
            record_kappa = False
            record_idh2oChemPot = True
        
        if record_idChemPot and len(elements)>0:
            for prefix in prefixes:
                if elements[0] == prefix:
                    idChemPot_log[prefix] = []
                    values = lines[counter+1].split()
                    for value in values:
                        idChemPot_log[prefix].append(value)
                    didChemPot_log[prefix] = []
                    values = lines[counter+2].split()
                    for value in values:
                        didChemPot_log[prefix].append(value)
        if record_idh2oChemPot and len(elements)>0:
            for prefix in prefixes:
                if elements[0] == prefix:
                    idh2oChemPot_log[prefix] = []
                    values = lines[counter+1].split()
                    for value in values:
                        idh2oChemPot_log[prefix].append(value)
                    didh2oChemPot_log[prefix] = []
                    values = lines[counter+2].split()
                    for value in values:
                        didh2oChemPot_log[prefix].append(value)
        if record_pChemPot and len(elements)>0:
            for prefix in prefixes:
                if elements[0] == prefix:
                    pChemPot_log[prefix] = []
                    values = lines[counter+1].split()
                    for value in values:
                        pChemPot_log[prefix].append(value)
                    dpChemPot_log[prefix] = []
                    values = lines[counter+2].split()
                    for value in values:
                        dpChemPot_log[prefix].append(value)
        if record_dens and len(elements)>0:
            for prefix in prefixes:
                if elements[0] == prefix:
                    dens_log[prefix] = []
                    values = lines[counter+1].split()
                    for value in values:
                        dens_log[prefix].append(value)
                    ddens_log[prefix] = []
                    values = lines[counter+2].split()
                    for value in values:
                        ddens_log[prefix].append(value)
        if record_hCap and len(elements)>0:
            for prefix in prefixes:
                if elements[0] == prefix:
                    hCap_log[prefix] = []
                    values = lines[counter+1].split()
                    for value in values:
                        hCap_log[prefix].append(value)
                    dhCap_log[prefix] = []
                    values = lines[counter+2].split()
                    for value in values:
                        dhCap_log[prefix].append(value)
        if record_hVap and len(elements)>0:
            for prefix in prefixes:
                if elements[0] == prefix:
                    hVap_log[prefix] = []
                    values = lines[counter+1].split()
                    for value in values:
                        hVap_log[prefix].append(value)
                    dhVap_log[prefix] = []
                    values = lines[counter+2].split()
                    for value in values:
                        dhVap_log[prefix].append(value)
        if record_kappa and len(elements)>0:
            for prefix in prefixes:
                if elements[0] == prefix:
                    kappa_log[prefix] = []
                    values = lines[counter+1].split()
                    for value in values:
                        kappa_log[prefix].append(value)
                    dkappa_log[prefix] = []
                    values = lines[counter+2].split()
                    for value in values:
                        dkappa_log[prefix].append(value)
        counter += 1
    return(current_params)
    

#========================================================================
# RESTART AT PARAMETERS WITH LOWEST OBJECTIVE FUNCTION
#========================================================================
def restart_lowest(logfile, initparam):
    lines_in = open(logfile, 'r')
    lines = lines_in.readlines()
    lines_in.close()
    counter = 0
    old_error = 100000000 #arbitrary number larger than the value of the objective function
    current_params = initparam
    for line in lines:
        elements = line.split()
        if len(elements)>3 and elements[0] == 'ERROR':
            subcounter = 0
            parameters = lines[counter-2].split()
            current_error = float(elements[-1])
            if current_error<old_error:
                old_error = current_error
                for parameter in parameters:
                    current_params[0,subcounter]=float(parameter)
                    subcounter += 1
        counter += 1
            
    log = open(logfile, 'a')
    log.write('\n RESTARTING OPTIMIZATION\n')
    log.close()
    return(current_params)
    

#========================================================================
# GENERATE NEW PARAMETERS
#========================================================================
# Can't use Newton's method because can't calculate mixed second derivatives
def gradient_descent(old_params,derivative):
    #pgrad = np.zeros(numparam, np.float64)
    #for x1 in range(numparam):
    #    pgrad[x1]=(objdict[x1+1]-objdict[x1+1+numparam])/(2*initperturb[x1])
    # Calculate a step size\
    boundaries = list()
    for param in old_params[0]:
        boundaries.append(((0.5)*param,1.5*param))
    param_step = 1/np.max(np.abs(derivative))**2
    current_params = old_params-param_step*derivative
    for param in range(len(current_params[0])):
        if current_params[0][param]>np.max(boundaries[param]):
            current_params[0][param] = np.max(boundaries[param])
        if current_params[0][param]<np.min(boundaries[param]):
            current_params[0][param] = np.min(boundaries[param])
    return (current_params)

#=========================================================================
# MAIN
#=========================================================================
if __name__ == "__main__":
    # Method can be 'derivative' or 'single'
    # Edit this stuff to be more resilient
    # INITIAL PARAMETERS
    (paramkey, initparam, numparam) = initial_parameters()
    # PARSE INPUTS
    (cluster, optimize, verbose, temp_t4900ewcoul, temp_pure, temp_idh2o, delT, press, nequil, logfile, figurefile, datafile, pbsfilename, single, gro_loc, job_script, orestart, directory) = parse_inputs()
    # SET PARAMETERS
    (prefixes, nsamples, cwd, dir_idh2o, dir_t4900ewcoul, dir_pure, dir_inputs_idh2o, dir_inputs_t4900ewcoul, dir_inputs_pure, dir_tops, g_energy, grompp, mdrun, ke_call, pe_call, v_call, energy_suffix, weight, pbs_lines, maxcount, errortol, error, rr_states, not_rr_states, n_systems) = set_parameters(single, gro_loc)
    # CHECK THAT FILES EXIST
    file_check(prefixes, dir_pure, dir_t4900ewcoul, dir_idh2o, dir_inputs_pure, dir_inputs_t4900ewcoul, dir_inputs_idh2o, dir_tops, optimize)

    errorsum_log = []

    idChemPot_log = dict()
    didChemPot_log = dict()

    idh2oChemPot_log = dict()
    didh2oChemPot_log = dict()

    pChemPot_log = dict()
    dpChemPot_log = dict()

    dens_log = dict()
    ddens_log = dict()

    hCap_log = dict()
    dhCap_log = dict()

    kappa_log = dict()
    dkappa_log = dict()

    hVap_log = dict()
    dhVap_log = dict()
    
    if (orestart) == False:
        # INITIALIZE LOG FILE
        initialize_log_file(cwd, logfile)
        current_params = initparam
    else:
        (current_params) = restart_optimization(logfile, datafile, initparam, prefixes)

    for prefix in prefixes:
        idChemPot_log[prefix] = list()
        didChemPot_log[prefix] = list()

        idh2oChemPot_log[prefix] = list()
        didh2oChemPot_log[prefix] = list()

        pChemPot_log[prefix] = list()
        dpChemPot_log[prefix] = list()

        dens_log[prefix] = list()
        ddens_log[prefix] = list()

        hCap_log[prefix] = list()
        dhCap_log[prefix] = list()

        kappa_log[prefix] = list()
        dkappa_log[prefix] = list()

        hVap_log[prefix] = list()
        dhVap_log[prefix] = list()

        method = 'L-BFGS-B'


    if optimize:
       import scipy.optimize as so
       boundaries = list()
       for param in initparam[0]:
           boundaries.append(((0.5)*param,1.5*param))
       so.fmin_tnc(main_single_point, current_params, fprime=main_derivatives, bounds=boundaries, disp=2)
       #pdb.set_trace()
       #main_resim(prefixes, dir_pure, dir_inputs_pure, dir_t4900ewcoul, dir_inputs_t4900ewcoul, dir_idh2o, dir_inputs_idh2o, cwd, paramkey, initparam, current_params, dir_tops, gro_loc, verbose, job_script, cluster, optimize)
        
    else:
        current_params = initparam
        main_single_point(current_params)
        #main_derivatives(current_params)

