"""
This file contains scripts for calculating J Couplings from backbone dihedrals.

References:

Self-consistent 3J coupling analysis for the joint calibration of Karplus coefficients and evaluation of torsion angles

Limits on variations in protein backbone dynamics from precise measurements of scalar couplings

Structure and dynamics of the homologous series of alanine peptides: a joint molecular dynamics/NMR study
"""
import numpy as np
import pandas as pd


def J3_HN_HA_schwalbe(phi):
    """
    RMS = 0.39
    Personal RMS on ubiquitin = 0.254
    Originally from Hu and Bax
    """
    phi = phi*np.pi/180.
    phi0 = -60*np.pi/180.
    A = 7.09
    B = -1.42
    C = 1.55
    return A*np.cos(phi + phi0)**2. + B*np.cos(phi + phi0) + C


def J3_HN_HA_ruterjans(phi):
    """RMS = 0.25
    """
    phi = phi*np.pi/180.
    phi0 = -60*np.pi/180.
    A = 7.90
    B = -1.05
    C = 0.65
    return A*np.cos(phi + phi0)**2. + B*np.cos(phi + phi0) + C


def J3_HN_HA_bax(phi):
    """RMS = 0.36
    """
    phi = phi*np.pi/180.
    phi0 = -60*np.pi/180.
    A = 8.4
    B = -1.36
    C = 0.33
    return A*np.cos(phi + phi0)**2. + B*np.cos(phi + phi0) + C


def J3_HN_HA_karplus(phi):
    """RMS = 0.36
    """
    phi = phi*np.pi/180.
    phi0 = 0.0 * np.pi/180.
    A = 6.4
    B = -1.4
    C = 1.9
    return A*np.cos(phi + phi0)**2. + B*np.cos(phi + phi0) + C


J3_HN_HA = J3_HN_HA_bax

uncertainties = pd.Series({
"J3_HN_HA":0.36, "J3_HN_Cprime":0.30, "J3_HN_CB":0.22, "J1_N_CA":0.53, "J2_N_CA":0.48,"J3_HA_Cprime":0.44
})
