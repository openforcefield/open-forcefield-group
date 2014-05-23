#!/usr/bin/env python
# -*- coding: utf-8 -*- 

#================================================#
#|                                              |#
#| Minimum-energy two dimensional torsion scans |#
#|         Lee-Ping Wang, October 2013          |#
#|                                              |#
#| This script launches multiple dihedral scans |#
#| and manages resources using the Work Queue   |#
#| library.                                     |#
#|                                              |#
#| One key aspect of this script is that it has |#
#| an internal queue which tries to finish      |#
#| individual scans while maximizing the use of |#
#| resources.                                   |#
#|                                              |#
#================================================#

import os, sys
import numpy as np
import itertools
from copy import deepcopy
from collections import defaultdict
from molecule import *
from nifty import _exec
from nifty import *

if __name__ == "__main__":
    # Create the Work Queue.  We will only have one instance because the
    # script will try to maximize the efficiency of resource utilization.
    wq_port = 7323
    work_queue.set_debug_flag('all')
    wq = work_queue.WorkQueue(port=wq_port, exclusive=False, shutdown=False)
    wq.tasks_failed = 0 # Counter for tasks that fail at the application level
    wq.specify_keepalive_interval(8640000)
    wq.specify_name('dihedral')
    print('Work Queue listening on %d' % (wq.port))

# Directory containing Python scripts and parameters.
based = os.path.split(os.path.abspath(__file__))[0]

# TINKER .key file template for running minimizations.
keyfile="""
parameters {prm}
digits 10
maxiter 10000

restrain-torsion {a} {b} {c} {d} {k1} {dih1}
restrain-torsion {e} {f} {g} {h} {k2} {dih2}
"""

# Q-Chem input file template for dual basis RI-MP2.
# with the 6-31G* basis set.
rimp2_qcfile="""
$molecule
{chg} {mult}
$end

$rem
jobtype              opt
exchange             hf
correlation          rimp2
basis                6-31g*
basis2               r64g
aux_basis            rimp2-vdz
dual_basis_energy    true
n_frozen_core        fc
mem_static           800
geom_opt_max_cycles  150
$end
"""

rimp2_qcfile_1="""
$molecule
{chg} {mult}
$end

$rem
jobtype              opt
exchange             hf
correlation          rimp2
basis                6-31+g*
aux_basis            aux_gen
n_frozen_core        fc
mem_static           800
geom_opt_max_cycles  150
$end

$aux_basis
H 0
rimp2-cc-pVDZ
****
C 0
rimp2-aug-cc-pVDZ
****
N 0
rimp2-aug-cc-pVDZ
****
O 0
rimp2-aug-cc-pVDZ
****
P 0
rimp2-aug-cc-pVDZ
****
S 0
rimp2-aug-cc-pVDZ
$end
"""

# Q-Chem input file template for B97-D
# with the 6-31G* basis set.
b97d_qcfile="""
$molecule
{chg} {mult}
$end

$rem
jobtype              opt
exchange             b97-d
basis                6-31{plus}g*
geom_opt_max_cycles  150
$end
"""

# Q-Chem input file template for Hartree-Fock
# with the 3-21G basis set.
hf_qcfile="""
$molecule
{chg} {mult}
$end

$rem
jobtype              opt
exchange             hf
basis                3-21g
geom_opt_max_cycles  150
$end
"""

chi12_qcfile="""

$opt
constraint
tors {a} {b} {c} {d} {phi}
tors {b} {c} {d} {e} {psi}
tors {b} {c} {f} {g} {chi1}
tors {c} {f} {g} {h} {chi2}
endconstraint
$end
"""

chi1_qcfile="""

$opt
constraint
tors {a} {b} {c} {d} {phi}
tors {b} {c} {d} {e} {psi}
tors {b} {c} {f} {g} {chi1}
endconstraint
$end
"""

dih12_qcfile="""

$opt
constraint
tors {a} {b} {c} {d} {dih1}
tors {e} {f} {g} {h} {dih2}
endconstraint
$end
"""

# Dictionary that maps methods to input files / programs.
MDict = {"AMOEBA"    : ("tinker", "amoebapro13.prm"),
         "AMBER99SB" : ("tinker", "amber99sb.prm"),
         "RIMP2"     : ("qchem", rimp2_qcfile),
         "B97D"      : ("qchem", b97d_qcfile),
         "HF"        : ("qchem", hf_qcfile),
         "None"      : (None, None)}

AngList = ["phipsi", "phichi1", "psichi1", "chi1chi2"]

# Listing of amino acid charges.
AAChg = {'LYS':1,
         'ARG':1,
         'ASP':-1,
         'GLU':-1,
         'SEP':-2}

# Sorted alphabetically but with glycine first
# No proline!!
AAList = ["GLY", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HID", "HIE", "ILE", "LEU", "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL", "SEP"]

# Sorted by number of atoms
AAList = ["GLY", "ALA", "CYS", "SER", "ASP", "ASN", "THR", "GLU", "VAL", "GLN", "MET", "HID", "HIE", "ILE", "LEU", "LYS", "PHE", "TYR", "ARG", "TRP", "SEP"][::-1]

AAList = ["GLY", "ALA", "CYS", "SER", "ASN", "THR", "VAL", "GLN", "MET", "HID", "HIE", "ILE", "LEU", "LYS", "PHE", "TYR", "ARG", "TRP", "SEP"][::-1]

AAList = ["ASP", "GLU"]

AAList = ["SEP"]

AAList = ["PRO"]

# Amino acids with scannable side chains.
Has_Chi1 = ["CYS", "SER", "ASP", "ASN", "THR", "GLU", "VAL", "GLN", "MET", "HID", "HIE", "ILE", "LEU", "LYS", "PHE", "TYR", "ARG", "TRP", "SEP"]
Has_Chi2 = ["CYS", "SER", "ASP", "ASN", "THR", "GLU", "GLN", "MET", "HID", "HIE", "ILE", "LEU", "LYS", "PHE", "TYR", "ARG", "TRP", "SEP"]

# Grid spacing in degrees.
grd = 15

# Range of values from -165 to +180 (in the case of a 15 degree grid.)
rng = np.array([int(i) for i in np.linspace(-180+grd, 180, 360/grd)])

# Empirical formula
def ef(Formula):
    return ''.join([('%s%i' % (k, Formula.count(k)) if Formula.count(k) > 1 else '%s' % k) for k in sorted(set(Formula))])

# Functions that make sure we stay within (-180, 180.)
def increment(angle, incr):
    new = angle + incr
    if new > 180: new -= 360
    if new <= -180: new += 360
    return new

# This function is used for measuring and we want to be biased towards +180.
def bracket(angle):
    new = angle
    if new > 181: new -= 360
    if new <= -179: new += 360
    return new

def bracket360(angle):
    new = angle
    if new >= 360: new -= 360
    if new < 0: new += 360
    return new

#=======================#
#| Managing topologies |#
#=======================#

def walk_bonds(MyG):
    """
    For each atom, walk along the topology and return a dictionary that looks like: 
    
    {number of bonds away: elements of atoms this many bonds away} 
    
    Useful for identifying atoms.
    """
    GDat = MyG.nodes(data=True)
    GDict = {}
    for i in GDat:
        GDict[i[0]] = i[1]
    
    PairPaths = nx.all_pairs_shortest_path_length(MyG)
    Walks = []
    Max   = 3
    for A in PairPaths:
        Walk = defaultdict(list)
        for B in PairPaths[A]:
            if PairPaths[A][B] > Max: continue
            Walk[PairPaths[A][B]].append(GDict[B]['e'])
        for idx, elist in Walk.items():
            Walk[idx] = ''.join([('%s%i' % (k, elist.count(k)) if elist.count(k) > 1 else '%s' % k) for k in sorted(set(elist))])
        Walks.append(Walk)
    return Walks

def iupac_chain(M, chain=[], chains=[], high=True):

    """ Given a starting chain of atoms in a molecule (i.e. alpha
    carbon and beta carbon), find the highest (or lowest) chain of
    atoms that starts from the beta carbon following atomic number
    rules in the IUPAC standard.

    This is a recursive function.

    Input: 
    
    M = ForceBalance Molecule object
    chain = Starting chain [c_alpha, c_beta]
    high = Set to True to get highest chains, set to False to get
    lowest chains
    
    Input/Output: chains = List of accumulated highest/lowest chains
    that are consistent with IUPAC [[c_alpha, c_beta, c_gamma1],
    [c_alpha, c_beta, c_gamma2]]
    """

    if chains == []: chains = [chain]
    bset = [set(i) for i in M.bonds]
    Branches = []
    BranchE = []
    # Find "branches" that contain atoms not already in the list of ones we have
    for k in range(M.na):
        if set((chain[-1], k)) in bset and k not in chain:
            Branches.append(k)
            BranchE.append(Elements.index(M.elem[k]))
    # Keep only the branches that have the highest/lowest atomic number.
    Gammas = []
    for k in range(len(Branches)):
        if BranchE[k] == (max(BranchE) if high else min(BranchE)):
            Gammas.append(Branches[k])
    if len(Gammas) >= 1:
        chains.remove(chain)
        for k in range(len(Gammas)):
            chains.append(chain+[Gammas[k]])
            iupac_chain(M, chain+[Gammas[k]], chains, high=high)
    # Keep only the chains that are equal to the highest/lowest chain.
    # Get the indices that would sort a list.
    s = [[Elements.index(M.elem[j]) for j in chains[i]] for i in range(len(chains))]
    cidx = sorted(range(len(chains)), key=lambda k: s[k])
    get = [Elements.index(M.elem[j]) for j in chains[cidx[-1 if high else 0]]]
    chains = [c for c in chains if [Elements.index(M.elem[j]) for j in c] == get]
    return chains

def iupac_sequence_rule(M, sn, n, ca, cb):
    """
    For two atoms A->B (e.g. alpha carbon and beta carbon), 
    determine the highest priority atom bonded to B
    for the purposes of torsion angle scanning.
    http://www.chem.qmul.ac.uk/iupac/misc/ppep1.html
    
    sn = Snapshot used in measuring (only for very edge cases)
    """
    Highest = iupac_chain(M, chain=[ca, cb], high=True)
    Lowest = iupac_chain(M, chain=[ca, cb], high=False)
    Gammas = sorted(list(set([i[2] for i in Highest])))
    Omegas = sorted(list(set([i[2] for i in Lowest])))

    bset = [set(i) for i in M.bonds]
    Branches = []
    BranchE = []
    # Find "branches" - i.e. bonded to the beta carbon
    for k in range(M.na):
        if set((cb, k)) in bset and k not in [ca]:
            Branches.append(k)
            BranchE.append(Elements.index(M.elem[k]))

    if len(Branches) > 0:
        # print "Atoms attached to bond %i -> %i ;" % (ca, cb), Branches, [M.elem[k] for k in Branches]
        if len(Gammas) == 1:
            # print "Found highest branch"
            return Gammas[0]
        else:
            print "Identical branches were found among", Gammas
            if len(Branches) == 3:
                if len(Gammas) == 3:
                    # All three branches are identical;
                    # should never be called for amino acids.
                    # Phis = []
                    # print "Three identical branches; invoking IUPAC rule 2.2.3 (smallest measured dihedral)"
                    # for k in range(len(Gammas)):
                    #     Phis.append(bracket(M.measure_dihedrals(n, ca, cb, Gammas[k])[sn]))
                    # return Gammas[np.argmin(np.abs(Phis))]
                    print "Three identical branches"
                    return Gammas
                elif len(Gammas) == 2:
                    # Two identical branches.
                    print "Invoking IUPAC rule 2.2.2 (LEU, VAL; clockwise from lowest branch)",
                    if len(Omegas) > 1:
                        raise Exception("With 2+ identical largest branches, there should only be one smallest branch")
                    Omega = Omegas[0]
                    PhiOmega = M.measure_dihedrals(n, ca, cb, Omega)[sn]
                    PhiGamma1 = M.measure_dihedrals(n, ca, cb, Gammas[0])[sn]
                    PhiGamma2 = M.measure_dihedrals(n, ca, cb, Gammas[1])[sn]
                    print PhiOmega, PhiGamma1, PhiGamma2
                    if bracket360(PhiGamma1 - PhiOmega) < bracket360(PhiGamma2 - PhiOmega):
                        return Gammas[0]
                    else:
                        return Gammas[1]
            if len(Branches) == 2:
                # Phis = []
                # print "Invoking IUPAC rule 2.3.2 (TYR / PHE; branch with smallest measured dihedral)"
                # for k in range(len(Gammas)):
                #     Phis.append(bracket(M.measure_dihedrals(n, ca, cb, Gammas[k])[sn]))
                # return Gammas[np.argmin(np.abs(Phis))]
                print "Two identical branches"
                return Gammas
            return Gammas[0]
    else:
        return

# Fingerprints that identify the "terminal" carbon and nitrogen atoms in a dipeptide.
fprint_ct = {0: 'C', 1: 'CNO', 2: 'CH4', 3: 'C2H'}
fprint_nt = {0: 'N', 1: 'C2H', 2: 'CH3O', 3: 'CHN'}
fprint_ctg = {0: 'C', 1: 'CNO', 2: 'CH4', 3: 'CH2'}
fprint_ntg = {0: 'N', 1: 'C2H', 2: 'CH3O', 3: 'H2N'}
fprint_ctp = {0: 'C', 1: 'CNO', 2: 'C2H3', 3: 'C3H3'}
fprint_cme = {0: 'C', 1: 'CH3', 2: 'NO', 3:'CH'}
fprint_cmep = {0: 'C', 1: 'CH3', 2: 'NO', 3:'C2'}
fprint_nme = {0: 'C', 1: 'H3N', 2: 'CH', 3:'CO'}
    
#==================#
#| Internal queue |#
#==================#

# Task identifier for internal queue.
# Elements of the tuple are:
# Name of the amino acid (plus dihedral angle identifier)
# Work Queue task object
# Scheduling priority (my own thing, not part of WQ)
# Priority is equal to 100*iteration + 20 - amino acid order
# Once the job is set to Work Queue, priority is set to -1.
# Work queue task ID; set to -1 prior to scheduling.
MyTask = namedtuple('MyTask', ['name', 'task', 'prio', 'id'])
IQTasks = []
def schedule(command, input_files, output_files, name="None", tag=None, prio=0, verbose=True):
    """ Create Work Queue task object and store in the internal queue. """
    global WQIDS
    task = work_queue.Task(command)
    for f in input_files:
        task.specify_input_file(f[0],f[1],cache=False)
    for f in output_files:
        task.specify_output_file(f[0],f[1],cache=False)
    task.specify_algorithm(work_queue.WORK_QUEUE_SCHEDULE_RAND)
    if tag == None:
        task.specify_tag(command)
    else:
        task.specify_tag(tag)
    IQTasks.append(MyTask(name, task, prio, -1))

def print_summary():
    """ Print Work Queue summary and internal queue which contains held tasks. """
    QCounts = OrderedDict(list(itertools.chain(*[[(i+'-'+j, 0) for i in AAList] for j in AngList])))
    HCounts = OrderedDict(list(itertools.chain(*[[(i+'-'+j, 0) for i in AAList] for j in AngList])))
    for T in IQTasks:
        if T.name not in QCounts:
            QCounts[T.name] = 0
            HCounts[T.name] = 0
        if T.prio == -1:
            QCounts[T.name] += 1
        else:
            HCounts[T.name] += 1
    print "Jobs sent to Work Queue:",
    print ', '.join(["%s x %i" % (AA, QCounts[AA]) for AA in QCounts.keys() if QCounts[AA] > 0])
    print "Jobs to be queued later:", 
    print ', '.join(["%s x %i" % (AA, HCounts[AA]) for AA in HCounts.keys() if HCounts[AA] > 0])

def manage_wq(wait_intvl=1, print_time=1200):
    """ 
    An adapted version of wq_wait1 from ForceBalance. 
    
    Main functions are:
    1) Dispatch highest-priority batches of jobs to the Work Queue if there are any workers free.
    2) If the Work Queue is full, then keep jobs in the internal queue (reservoir).
    3) Submit a single job to Work Queue at the very beginning.
    4) Perform book-keeping to keep track of which jobs are in the Work Queue and which ones are in the reservoir.
    """

    cooldown_time = 10

    # No tasks to schedule
    if len(IQTasks) == 0: return

    # Determine if any jobs need to be submitted.
    jobs_needed = wq.stats.workers_ready
    if wq.stats.total_tasks_dispatched == 0 and jobs_needed == 0:
        jobs_needed = 1
        
    # Submit all jobs that have the highest priority; this may
    # overfill the queue slightly, and that is intentional.
    # All jobs are kept in the IQTasks dictionary.
    while True:
        highest_prio = sorted(list(t.prio for t in IQTasks))[-1]
        if time.time() - manage_wq.cool < cooldown_time: break # Cooldown time of 10 seconds for submitting batches of jobs.
        if highest_prio < 0: break
        if jobs_needed <= 0: break
        for T in IQTasks:
            # id == -1 indicates that the task hasn't been submitted yet.
            if T.prio == highest_prio and T.id == -1:
                taskid = wq.submit(T.task)
                # print "Submitted task: Name %s Priority %i taskid %i jobs_needed %i" % (T.name, T.prio, taskid, jobs_needed)
                # Replace the namedtuple with a new one that has priority -1 and 
                # appropriate taskid signifying that the job has been submitted.
                # Doing this because namedtuples are immutable -_-
                IQTasks[IQTasks.index(T)] = MyTask(T.name, T.task, -1, taskid)
                jobs_needed -= 1
                if jobs_needed <= 0:
                    manage_wq.cool = time.time()

    # If any workers are connected, this will wait for ten seconds.
    task = wq.wait(wait_intvl)
    if task:
        exectime = task.cmd_execution_time/1000000
        if task.result != 0:
            oldid = task.id
            oldhost = task.hostname
            # If the task has failed, resubmit it and update the book-keeping in IQTasks.
            # The task should be replaced exactly once - otherwise raise an exception.
            taskid = wq.submit(task)
            Replaced = False
            for T in IQTasks:
                if T.id == oldid:
                    if Replaced: raise Exception('Already replaced')
                    IQTasks[IQTasks.index(T)] = MyTask(T.name, task, -1, taskid)
                    Replaced = True
            if not Replaced: raise Exception('Did not replace')
            logger.warning("Command '%s' (task %i) failed on host %s (%i seconds), resubmitted: taskid %i\n" % (task.command, oldid, oldhost, exectime, taskid))
            wq.tasks_failed += 1
        else:
            if exectime > print_time: # Assume that we're only interested in printing jobs that last longer than a minute.
                logger.info("Command '%s' (task %i) finished successfully on host %s (%i seconds)\n" % (task.command, task.id, task.hostname, exectime))
            Deleted = False
            # If the task has succeeded, delete it and update the book-keeping in IQTasks.
            # The task should be deleted exactly once - otherwise raise an exception.
            for T in IQTasks:
                if T.id == task.id:
                    if Deleted: raise Exception('Already deleted')
                    IQTasks.remove(T)
                    Deleted = True
            if not Deleted: raise Exception('Did not delete')
            del task

    # Print some information showing work queue status and number of jobs in reservoir.
    try:
        # Full workers were added with CCTools 4.0.1
        nbusy = wq.stats.workers_busy + wq.stats.workers_full
    except:
        nbusy = wq.stats.workers_busy
    try:
        Complete = wq.stats.total_tasks_complete - wq.tasks_failed
        Total = wq.stats.total_tasks_dispatched - wq.tasks_failed
    except:
        logger.warning("wq object has no tasks_failed attribute, please use createWorkQueue() function.\n")
        Complete = wq.stats.total_tasks_complete
        Total = wq.stats.total_tasks_dispatched
    logger.info("\r%s : %i/%i workers busy; %i/%i jobs complete\r" % \
                    (time.ctime(), nbusy, (wq.stats.total_workers_joined - wq.stats.total_workers_removed), Complete, Total)) 
    # Print a newline if some time has elapsed.
    if time.time() - manage_wq.t0 > 60:
        manage_wq.t0 = time.time()
        logger.info('\n')
        print_summary()
manage_wq.cool = time.time()
manage_wq.t0 = time.time()
    
def identify(M, name=None):
    """
    Given a Molecule object, this function identifies the ACE carbon,
    the backbone nitrogen (N-H), the alpha carbon, the carbonyl carbon
    (C=O), and the NME nitrogen.  Also determines chirality of alpha carbon.
    """
    cme = None
    nme = None
    ct = None
    n1 = None
    ca1 = None
    c1 = None
    nt = None
    h1 = None
    cb1 = None
    walks = walk_bonds(M.topology)
    # First figure out acetyl carbon and N-methyl nitrogen.
    for i in range(M.na):
        if dict(walks[i]) == fprint_ct:
            if ct == None:
                ct = i
            else:
                raise Exception('Found duplicate ACE carbon')
        if dict(walks[i]) == fprint_nt:
            if nt == None:
                nt = i
            else:
                raise Exception('Found duplicate NME nitrogen')
    
    Glycine = False
    Proline = False
    # Glycine special case
    if ct == None and nt == None:
        Glycine = True
        for i in range(M.na):
            if dict(walks[i]) == fprint_ctg:
                if ct == None:
                    ct = i
                else:
                    raise Exception('Found duplicate ACE carbon')
            if dict(walks[i]) == fprint_ntg:
                if nt == None:
                    nt = i
                else:
                    raise Exception('Found duplicate NME nitrogen')
    # Proline special case
    elif ct == None:
        Proline = True
        for i in range(M.na):
            if dict(walks[i]) == fprint_ctp:
                if ct == None:
                    ct = i
                else:
                    raise Exception('Found duplicate ACE carbon')

    # Figure out terminal methyl carbons.
    for i in range(M.na):
        if (not Proline and dict(walks[i]) == fprint_cme) or (Proline and dict(walks[i]) == fprint_cmep):
            if cme == None:
                cme = i
            else:
                raise Exception('Found duplicate ACE methyl')
        if dict(walks[i]) == fprint_nme:
            if nme == None:
                nme = i
            else:
                raise Exception('Found duplicate NME methyl')
    # print cme, nme
    
    # Next figure out C=O and NH.
    for i, j in M.bonds:
        for k in range(M.na):
            if set((k, ct)) == set((i, j)) and M.elem[k] == "N":
                if n1 == None:
                    n1 = k
                else:
                    raise Exception('Found duplicate NH')
            if set((k, nt)) == set((i, j)) and M.elem[k] == "C" and walks[k][1] == 'CNO':
                if c1 == None:
                    c1 = k
                else:
                    raise Exception('Found duplicate C=O')
    # Finally figure out alpha carbon.
    bset = [set(i) for i in M.bonds]
    for k in range(M.na):
        if set((n1, k)) in bset and set((c1, k)) in bset:
            if ca1 == None:
                ca1 = k
            else:
                raise Exception('Found duplicate alpha carbon')
    
    # Figure out alpha hydrogen.
    for k in range(M.na):
        if M.elem[k] == 'H' and set((ca1, k)) in bset:
            if h1 == None:
                h1 = k
            elif Glycine and cb1 == None:
                cb1 = k
            else:
                raise Exception('Found duplicate alpha hydrogen')
    
    # Figure out beta carbon.
    for k in range(M.na):
        if M.elem[k] != 'H' and set((ca1, k)) in bset and k not in [n1, c1, h1]:
            if cb1 == None:
                cb1 = k
            else:
                raise Exception('Found duplicate beta carbon')
    
    # Figure out chirality.
    # Vector going from alpha carbon to C=O.
    caco = M.xyzs[0][c1] - M.xyzs[0][ca1]
    caco /= np.linalg.norm(caco)
    # Vector going from alpha carbon to N-H.
    can = M.xyzs[0][n1] - M.xyzs[0][ca1]
    can /= np.linalg.norm(can)
    # Vector going from alpha carbon to H.
    cah = M.xyzs[0][h1] - M.xyzs[0][ca1]
    cah /= np.linalg.norm(cah)
    # Cross product.
    xcon = np.cross(caco, can)
    xcon /= np.linalg.norm(xcon)
    ang = 180.0*np.arccos(np.dot(xcon, cah))/np.pi
    chiral = "D" if (not Glycine and ang > 90.0) else "L"

    if chiral == "D":
        print "Correcting chirality"
        for i in range(len(M)):
            M.xyzs[i][:,0] *= -1

    gamma1 = None
    delta1 = None
    if name in Has_Chi1:
        gamma1 = iupac_sequence_rule(M, 0, n1, ca1, cb1)
        if name in Has_Chi2:
            delta1 = iupac_sequence_rule(M, 0, ca1, cb1, gamma1)

    return ct, n1, ca1, c1, nt, cb1, h1, cme, nme, gamma1, delta1, chiral

# Figures out the four neighbors of the dihedral point.
def neighbors(dih1, dih2):
    if not (dih1 in rng and dih2 in rng):
        raise
    return [(increment(dih1, -grd), dih2),
            (dih1, increment(dih2, -grd)),
            (increment(dih1, +grd), dih2),
            (dih1, increment(dih2, +grd))]

def unmangle(M1, M2):
    """ 
    Create a mapping that takes M1's atom indices to M2's atom indices based on position.  
    
    If we start with atoms in molecule "PDB", and the new molecule "M" 
    contains re-numbered atoms, then this code works:
    
    M.elem = list(np.array(PDB.elem)[unmangled])
    """
    if len(M1) != 1 or len(M2) != 1:
        print "I only deal with length-1 molecule objects"
        sys.exit()
    unmangler = {}
    for i in range(M1.na):
        for j in range(M2.na):
            if np.linalg.norm(M1.xyzs[0][i] - M2.xyzs[0][j]) < 0.1:
                unmangler[j] = i
    unmangled = [unmangler[i] for i in sorted(unmangler.keys())]
    return unmangled

def read_psi(psiout):
    # Read Psi4 output file for geometries and final CBS energy.
    XMode = 0
    EMode = 0
    xyzs = []
    xyz = []
    elem = []
    cbs = 0
    for line in open(psiout):
        s = line.split()
        if XMode == 1:
            if len(s) == 4 and isfloat(s[1]) and isfloat(s[2]) and isfloat(s[3]):
                e = s[0]
                xyz.append([float(i) for i in s[1:4]])
                if EMode == 1:
                    elem.append(e)
            elif len(xyz) > 0:
                xyzs.append(np.array(xyz))
                xyz = []
                XMode = 0
        if line.strip().startswith("Geometry (in Angstrom)"):
            XMode = 1
            EMode = len(elem) == 0
        if len(s) == 3 and s[0] == 'total' and s[1] == 'CBS':
            cbs = float(s[2])
    if len(xyzs) == 0:
        raise Exception('%s has length zero' % psiout)
    return xyzs, elem, cbs

def read_psi_grad(psiout):
    # Read Psi4 output file for geometries and final CBS energy.
    XMode = 0
    EMode = 0
    GMode = 0
    gxyzs = []
    gxyz = []
    xyzs = []
    xyz = []
    elem = []
    cbs = 0
    for line in open(psiout):
        s = line.split()
        if XMode == 1:
            if len(s) == 4 and isfloat(s[1]) and isfloat(s[2]) and isfloat(s[3]):
                e = s[0]
                xyz.append([float(i) for i in s[1:4]])
                if EMode == 1:
                    elem.append(e)
            elif len(xyz) > 0:
                xyzs.append(np.array(xyz))
                xyz = []
                XMode = 0
        if line.strip().startswith("Geometry (in Angstrom)"):
            XMode = 1
            EMode = len(elem) == 0
        if GMode == 1:
            if len(s) == 4 and isfloat(s[1]) and isfloat(s[2]) and isfloat(s[3]):
                gxyz.append([float(i) for i in s[1:4]])
            elif len(gxyz) > 0:
                gxyzs.append(np.array(gxyz))
                gxyz = []
                GMode = 0
        if line.strip().startswith("-Total Gradient:"):
            GMode = 1
        if len(s) == 3 and s[0] == 'total' and s[1] == 'CBS':
            cbs = float(s[2])
    if len(xyzs) == 0:
        raise Exception('%s has length zero' % psiout)
    return xyzs, elem, gxyzs, cbs

class DihedralGrid(object):
    """
    Class which represents a two-dimensional grid of dihedral angles for a molecule.
    Provides methods for iteratively minimizing the energy on the two-dimensional surface
    until a continuous surface is obtained.
    """
    #------
    # The old constructor, specific to amino acids
    #------
    # def __init__(self, fnm, method, angles, arc=None, root=None, name="NoName"):
    #     # Load in the original coordinates as a PDB.
    #     PDB = Molecule(fnm)
    #     if arc != None: XYZ = Molecule(arc)
    #     else: XYZ = PDB
    #     self.name = name
    #     # Set Work Queue priority based on amino acid name.
    #     if name in AAList:
    #         self.prio = len(AAList) - AAList.index(name)
    #     else:
    #         self.prio = 0
    #     if root == None:
    #         # Root folder containing the PDB file.
    #         self.root = os.path.split(os.path.abspath(fnm))[0]
    #     else:
    #         self.root = root
    #     # Initialize energy and geometry storage.
    #     self.energies = OrderedDict([(i, None) for i in list(itertools.product(rng, repeat=2))])
    #     self.geoms = OrderedDict([(i, None) for i in list(itertools.product(rng, repeat=2))])
    #     # List of geometries to be optimized next.
    #     # Actually a list of lists [[(phi1a, psi1a), (phi1b, psi1b)], (phi0, psi0)]
    #     # indicating the constraint points and the initial geometry.
    #     self.optnext = []
    #     # Record of which points were reoptimized and found to be lower.
    #     self.lower = []
    #     # Method (a few are hard coded in this script)
    #     self.method = method
    #     # Angles to scan over (choose from phi-psi, phi-chi1, psi-chi1, chi1-chi2)
    #     self.angles = angles
    #     self.aname = self.name + '-' + self.angles
    #     # Software and template files for the method
    #     self.engine, self.template = MDict[method]
    #     if self.method == "AMBER99SB":
    #         self.tprm = "amber99sb.prm"
    #     else:
    #         self.tprm = "amoebapro13.prm"
    #     # State of the object (current iteration, and is it running any jobs?)
    #     self.iteration = 0
    #     self.running = 0
    #     self.stordir = []
    #     # Call TINKER to write a TINKER-formatted .xyz file using the
    #     # atom types from the force field.  This is because TINKER
    #     # requires specially formatted files with the atom type on
    #     # each line.
    #     _exec("rm -f *.xyz* *.seq*", print_command=False)
    #     tpdb = os.path.join(self.root, "Scan", self.angles, self.method, "temp.pdb")
    #     txyz = os.path.join(self.root, "Scan", self.angles, self.method, "temp.xyz")
    #     mdnm = os.path.join(self.root, "Scan", self.angles, self.method)
    #     if not os.path.exists(mdnm): os.makedirs(mdnm)
    #     _exec("rm -f %s/temp.*" % os.path.join(self.root, "Scan", self.angles, self.method), print_command=False)
    #     PDB[0].write(tpdb)
    #     o = _exec("pdbxyz %s %s/%s" % (tpdb, based, self.tprm), print_command=False)
    #     # TINKER mangles atom ordering so we go with TINKER atom orders
    #     self.M          = Molecule(txyz, ftype="tinker")
    #     unmangled       = unmangle(PDB[0], self.M)
    #     self.M.elem     = list(np.array(PDB.elem)[unmangled])
    #     self.M.atomname = list(np.array(PDB.atomname)[unmangled])
    #     self.M.resid    = list(np.array(PDB.resid)[unmangled])
    #     self.M.resname  = list(np.array(PDB.resname)[unmangled])
    #     self.M.xyzs     = list([XYZ.xyzs[i][unmangled] for i in range(len(XYZ))])
    #     self.M.comms    = PDB.comms
    #     self.M[0].write(os.path.join(self.root, "Scan", self.angles, self.method, "unmangled.pdb"))
    #     # Determine the charge on the molecule (this works for tripeptides and beyond)
    #     rid = 0
    #     # self.chg = 0
    #     # Don't depend on the amino acid labeling in the PDB file, it can be wrong
    #     self.chg = AAChg.get(self.name, 0)
    #     self.mult = 1
    #     # for i in range(self.M.na):
    #     #     if self.M.resid[i] != rid:
    #     #         rid = self.M.resid[i]
    #     #         if self.M.resname[i] in AAChg:
    #     #             self.chg += AAChg[self.M.resname[i]]
    #     # Determine the important atoms and correct chirality (only works for dipeptides)
    #     self.ct, self.n1, self.ca1, self.c1, self.nt, self.cb1, self.h1, self.cme, self.nme, self.gam1, self.del1, chiral = identify(self.M, self.name)
    #     # Sometimes there are two identical delta-atoms.
    #     if type(self.del1) is list:
    #         self.del1 = self.del1[0]
    #     # Set the molecule object.
    #     self.M = M
    #     phiatoms = (self.ct, self.n1, self.ca1, self.c1)
    #     psiatoms = (self.n1, self.ca1, self.c1, self.nt)
    #     chi1atoms = (self.n1, self.ca1, self.cb1, self.gam1)
    #     chi2atoms = (self.ca1, self.cb1, self.gam1, self.del1)
    #     if self.angles == "phipsi":
    #         self.iA, self.iB, self.iC, self.iD = phiatoms
    #         self.iE, self.iF, self.iG, self.iH = psiatoms
    #     if self.angles == "phichi1":
    #         self.iA, self.iB, self.iC, self.iD = phiatoms
    #         self.iE, self.iF, self.iG, self.iH = chi1atoms
    #     if self.angles == "psichi1":
    #         self.iA, self.iB, self.iC, self.iD = psiatoms
    #         self.iE, self.iF, self.iG, self.iH = chi1atoms
    #     if self.angles == "chi1chi2":
    #         self.iA, self.iB, self.iC, self.iD = chi1atoms
    #         self.iE, self.iF, self.iG, self.iH = chi2atoms
    #     self.valid = all([i != None for i in [self.iA, self.iB, self.iC, self.iD, self.iE, self.iF, self.iG, self.iH]])
    #     # Switch to augmented basis (haDZ) for negatively charged residues.
    #     if self.method == "RIMP2" and self.chg < 0:
    #         self.template = rimp2_qcfile_1

    def __init__(self, M, chg, mult, method, da, db, name):
        # Initialize energy and geometry storage.
        self.energies = OrderedDict([(i, None) for i in list(itertools.product(rng, repeat=2))])
        self.geoms = OrderedDict([(i, None) for i in list(itertools.product(rng, repeat=2))])
        # Set some labels and folders.
        self.name = name
        self.root = os.path.join("Scan", name)
        # Molecule object.
        self.M = M
        # List of geometries to be optimized next.
        # Actually a list of lists [[(phi1a, psi1a), (phi1b, psi1b)], (phi0, psi0)]
        # indicating the constraint points and the initial geometry.
        self.optnext = []
        # Record of which points were reoptimized and found to be lower.
        self.lower = []
        # Method (a few are hard coded in this script)
        self.method = method
        self.engine, self.template = MDict[method]
        if self.engine == "tinker":
            raise RuntimeError("TINKER not supported yet (need to read PDB files, unmangle atoms and whatnot)")
        # State of the object (current iteration, and is it running any jobs?)
        self.iteration = 0
        self.running = 0
        self.stordir = []
        # Set charge and multiplicity
        self.chg = chg
        self.mult = mult
        # Set atoms to be constrained.
        self.iA, self.iB, self.iC, self.iD = da
        self.iE, self.iF, self.iG, self.iH = db
        if self.iB == self.iE and self.iC == self.iF and self.iD == self.iG:
            self.angles = "%i-%i-%i-%i-%i" % (self.iA, self.iB, self.iC, self.iD, self.iH)
        else:
            self.angles = "%i-%i-%i-%i.%i-%i-%i-%i" % (self.iA, self.iB, self.iC, self.iD, self.iE, self.iF, self.iG, self.iH)
        # Specifies the name of the job. (More useful when we have many scans)
        self.aname = self.name + '.' + self.angles
        # Specifies the priority of the job. (More useful when we have many scans)
        self.prio = 0
        # Switch to augmented basis (haDZ) for negatively charged molecules.
        if self.method == "RIMP2" and self.chg < 0:
            self.template = rimp2_qcfile_1

    def run_minimize(self, dih12=None, dih120=None, frame=0):
        """ 
        Run a single minimization job (actually, determine if the job
        has already been run and queue it up if necessary.)
        
        dih12 = Perform the optimization with dihedrals constrained
        to the (phi, psi) values at this grid point.  If not provided,
        use the grid point closest to the (phi, psi) value of the
        initial geometry.

        dih120 = Use the stored geometry at this grid point to
        initialize the optimization.  If not provided, use the initial
        geometry.
        
        This function will queue up the optimization job and return
        the working directory and dihedral angles for the queued job.
        If it detects that the exact optimization has already been
        done, then it will not queue up the job but still return the
        directory (the read_minimize function won't know the difference)
        """
        M0 = self.M[0]
        if dih12 == None:
            # If the restraint dihedrals are not provided,
            # use the closest grid point to the current dihedrals.
            M0.xyzs[0] = self.M.xyzs[frame]
            dih1m = bracket(M0.measure_dihedrals(self.iA, self.iB, self.iC, self.iD)[0])
            dih2m = bracket(M0.measure_dihedrals(self.iE, self.iF, self.iG, self.iH)[0])
            dih1 = rng[np.argmin(np.abs(rng - dih1m))]
            dih2 = rng[np.argmin(np.abs(rng - dih2m))]
            dih12 = (dih1, dih2)
        else:
            dih1, dih2 = dih12
        if dih12 not in self.geoms: # This should never happen...
            raise Exception('dih12 = %s is not on the 15-degree grid' % str(dih12))
        if dih120 != None:
            # Look up initial geometry from dictionary of stored
            # geometries that were saved from the previous iteration.
            M0.xyzs[0] = self.geoms[dih120]
            if M0.xyzs[0] == None:
                raise Exception('dih120 = %s is not a stored geometry' % str(dih120))
        # Actually run the minimization.
        # The directory structure looks like "ALA/Scan/HF/+180/-060/5" which correspond to: 
        # Name of the amino acid
        # Directory containing scans (may name after the angles being scanned)
        # Method used to do the scan
        # The dihedral angles for the constraint
        # A sequential directory number for storing past calculations
        dn = 0
        dnm = os.path.join(self.root, self.angles, self.method, "%+04i" % dih1, "%+04i" % dih2, "%i" % dn)
        while os.path.exists(dnm):
            # Look over output files that already exist.
            # If we discover a finished optimization with the same starting (or ending) 
            # point that we are about to do, then we don't need to actually run the optimization.
            # If we need to run a new optimization, then create a new directory following the sequence.
            if self.engine == "tinker":
                if dnm not in self.stordir and all([os.path.exists(os.path.join(dnm,i)) for i in ["temp.xyz_2", "optimize.log"]]):
                    try:
                        M1 = Molecule(os.path.join(dnm,"temp.xyz"), ftype="tinker", build_topology=False)
                        M2 = Molecule(os.path.join(dnm,"temp.xyz_2"), ftype="tinker", build_topology=False)
                        max1 = abs(M1.xyzs[0] - M0.xyzs[0]).max()
                        max2 = abs(M2.xyzs[0] - M0.xyzs[0]).max()
                        if max1 < 0.01 or max2 < 0.01:
                            E = [float(line.split()[-1]) for line in open(os.path.join(dnm,"optimize.log")) if "Potential at Optimized Geometry:" in line][0]
                            print "\rFound: %12s -> %12s" % (dih120, dih12),
                            self.stordir.append(dnm)
                            return dih12, dnm
                    except: pass
            if self.engine == "qchem":
                if dnm not in self.stordir and all([os.path.exists(os.path.join(dnm,i)) for i in ["energy.txt", "opt.xyz"]]):
                    try:
                        M2 = Molecule(os.path.join(dnm,"opt.xyz"), build_topology=False)
                        max2 = abs(M2.xyzs[0] - M0.xyzs[0]).max()
                        if max2 < 0.01:
                            E = 627.51*float(np.loadtxt(os.path.join(dnm,"energy.txt")))
                            print "\rFound: %12s -> %12s" % (dih120, dih12),
                            self.stordir.append(dnm)
                            return dih12, dnm
                        M1 = Molecule(os.path.join(dnm,"qchem.in"), ftype="qcin", build_topology=False)
                        max1 = abs(M1.xyzs[0] - M0.xyzs[0]).max()
                        if max1 < 0.01:
                            E = 627.51*float(np.loadtxt(os.path.join(dnm,"energy.txt")))
                            print "\rFound: %12s -> %12s" % (dih120, dih12),
                            self.stordir.append(dnm)
                            return dih12, dnm
                    except: pass
            dn += 1
            dnm = os.path.join(self.root, self.angles, self.method, "%+04i" % dih1, "%+04i" % dih2, "%i" % dn)
        os.makedirs(dnm)
        print "\rLaunch: %12s -> %12s" % (dih120, dih12),
        # Submit the job to the Work Queue.
        if self.engine == "tinker":
            M0.write(os.path.join(dnm,'temp.xyz'),ftype="tinker")
            with open(os.path.join(dnm,'temp.key'),'w') as f: print >> f, \
                    keyfile.format(a=self.iA+1, b=self.iB+1, c=self.iC+1, d=self.iD+1, 
                                   e=self.iE+1, f=self.iF+1, g=self.iG+1, h=self.iH+1,
                                   k1=100, k2=100, dih1=dih1, dih2=dih2, prm=self.template)
            schedule("python opt-tinker.py &> optimize.log", verbose=False,
                     input_files=[(os.path.join(dnm,"temp.xyz"),"temp.xyz"),
                                  (os.path.join(dnm,"temp.key"),"temp.key"),
                                  (os.path.join(based,self.tprm),self.tprm),
                                  (os.path.join(based,"opt-tinker.py"),"opt-tinker.py")],
                     output_files=[(os.path.join(dnm,"optimize.log"),"optimize.log"),
                                   (os.path.join(dnm,"temp.xyz_2"),"temp.xyz_2")],
                     name=self.aname, prio=self.prio+self.iteration*100)
        elif self.engine=="qchem":
            with open(os.path.join(dnm,'qtemp.in'),'w') as f: print >> f, \
                    self.template.format(chg=self.chg, mult=self.mult, plus="+" if self.chg < 0 else "")
            with open(os.path.join(dnm,'qtemp.in'),'a') as f: print >> f, \
                    dih12_qcfile.format(a=self.iA+1, b=self.iB+1, c=self.iC+1, d=self.iD+1, 
                                        e=self.iE+1, f=self.iF+1, g=self.iG+1, h=self.iH+1,
                                        dih1=dih1, dih2=dih2)

            Q = deepcopy(M0)
            Q.add_quantum(os.path.join(dnm,'qtemp.in'))
            Q.write(os.path.join(dnm,'qchem.in'), ftype="qcin")
            schedule("python opt-qchem.py", verbose=False, # Single core
                     input_files=[(os.path.join(dnm,"qchem.in"),"qchem.in"),
                                  (os.path.join(based,"opt-qchem.py"),"opt-qchem.py")],
                     output_files=[(os.path.join(dnm,"qchem.out.bz2"),"qchem.out.bz2"),
                                   (os.path.join(dnm,"energy.txt"),"energy.txt"),
                                   (os.path.join(dnm,"opt.xyz"),"opt.xyz")],
                     name=self.aname, tag=dnm, prio=self.prio+self.iteration*100)
        else:
            print "Only the TINKER or Q-Chem engines are supported right now"
            sys.exit()
        return dih12, dnm

    def read_minimize(self, dih12, dnm):
        """ Read the output of a single minimization job. """
        dih1, dih2 = dih12
        worked = True
        # Try to read the output data in the directory.
        # By the time this function is called, all of the jobs ought to be done.
        if self.engine == "tinker":
            try:
                E = [float(line.split()[-1]) for line in open(os.path.join(dnm,"optimize.log")) if "Potential at Optimized Geometry:" in line][0]
                M1 = Molecule(os.path.join(dnm,"temp.xyz_2"), ftype="tinker", build_topology=True)
            except:
                print dih12, "optimization failed"
                worked = False
        elif self.engine == "qchem":
            try:
                E = 627.51*float(np.loadtxt(os.path.join(dnm,"energy.txt")))
                M1 = Molecule(os.path.join(dnm,"opt.xyz"))
            except:
                print dih12, "optimization failed"
                worked = False
        else:
            print "Only the TINKER or Q-Chem engines are supported right now"
            sys.exit()
        # Get the status of the minimization, and store the converged energy and geometry.
        # -1 : The optimization failed.
        #  0 : The optimization converged and it's the first data point we have for this dihedral angle.
        #  1 : The optimization converged to a lower energy than what we had stored previously.
        #  In cases (0) and (1), we launch new optimizations at the neighboring points.
        #  This is how we "iron out" all of the discontinuities.
        #  2 : The optimization converged to an equal or higher energy.
        status = -1
        if worked:
            status = 0
            if M1.bonds != self.M.bonds:
                status = 2
            elif self.energies[dih12] == None:
                self.energies[dih12] = E
                self.optnext.append([neighbors(dih1, dih2), (dih1, dih2)])
                self.geoms[dih12] = M1.xyzs[0]
            elif E + 0.01 < self.energies[dih12]: 
                status = 1
                self.energies[dih12] = E
                if [neighbors(dih1, dih2), (dih1, dih2)] not in self.optnext:
                    self.optnext.append([neighbors(dih1, dih2), (dih1, dih2)])
                self.geoms[dih12] = M1.xyzs[0]
            elif E - 0.01 < self.energies[dih12]:
                status = 2
            else:
                status = 2
        return status

    def print_ascii_image(self):
        # Print an ascii image that shows the progress of the optimization.
        queued = OrderedDict([(i, 0) for i in list(itertools.product(rng, repeat=2))])
        haves = OrderedDict([(i, 0) for i in list(itertools.product(rng, repeat=2))])
        for fys, fy0 in self.optnext:
            for fy in fys:
                if fy[0] - fy0[0] in [-grd, 360-grd, -360-grd]: 
                    if queued[fy] == 0:
                        queued[fy] = 1                            # Left arrow.
                    elif queued[fy] == 3:
                        queued[fy] = 13                           # Down-left.
                    elif queued[fy] == 4:
                        queued[fy] = 14                           # Up-left.
                    else:
                        queued[fy] = 5
                elif fy[0] - fy0[0] in [grd, 360+grd, -360+grd]:  
                    if queued[fy] == 0:
                        queued[fy] = 2                            # Right arrow.
                    elif queued[fy] == 3:
                        queued[fy] = 23                           # Down-left.
                    elif queued[fy] == 4:
                        queued[fy] = 24                           # Up-left.
                    else:
                        queued[fy] = 5
                elif fy[1] - fy0[1] in [-grd, 360-grd, -360-grd]:
                    if queued[fy] == 0:
                        queued[fy] = 3
                    elif queued[fy] == 1:
                        queued[fy] = 13
                    elif queued[fy] == 2:
                        queued[fy] = 23
                    else:
                        queued[fy] = 5
                elif fy[1] - fy0[1] in [grd, 360+grd, -360+grd]: 
                    if queued[fy] == 0:
                        queued[fy] = 4
                    elif queued[fy] == 1:
                        queued[fy] = 14
                    elif queued[fy] == 2:
                        queued[fy] = 24
                    else:
                        queued[fy] = 5
                else:
                    queued[fy] = 5
        for i in self.energies:
            if self.energies[i]:
                haves[i] = 1
        print
        print "--== Ramachandran Plot of Optimization Status for %s ==--" % self.aname
        print "--==  Blue: Optimized, Green: Found Lower, Red: Next  ==--"
        print " ",
        for i in rng[::3]:
            print "%5i" % i,
        print
        for i in rng[::-1]:
            line = '%4i ' % i
            # Unicode characters, woohoo!
            for j in rng:
                if queued[(j,i)] == 1:
                    line += '\x1b[1;41m＜\x1b[0m'
                elif queued[(j,i)] == 2:
                    line += '\x1b[1;41m＞\x1b[0m'
                elif queued[(j,i)] == 3:
                    line += '\x1b[1;41m\\/\x1b[0m'
                elif queued[(j,i)] == 4:
                    line += '\x1b[1;41m/\\\x1b[0m'
                elif queued[(j,i)] == 13:
                    line += '\x1b[41m＼\x1b[0m'
                elif queued[(j,i)] == 14:
                    line += '\x1b[41m／\x1b[0m'
                elif queued[(j,i)] == 23:
                    line += '\x1b[41m／\x1b[0m'
                elif queued[(j,i)] == 24:
                    line += '\x1b[41m＼\x1b[0m'
                elif queued[(j,i)] == 5:
                    line += '\x1b[1;41m><\x1b[0m'
                elif (j,i) in self.lower:
                    line += '\x1b[42m--\x1b[0m'
                elif haves[(j,i)]:
                    line += '\x1b[44m--\x1b[0m'
                else:
                    line += '  '
            print line

    def launch_iteration(self):
        """ Launch an iteration (which looks like a wavefront propagation on the dihedral surface). """
        # We should not start any new calculations if there are already some running.
        if self.running: return 1
        # This list contains the phi-psi values and working directory of each minimization.
        self.optinfo = []
        if self.iteration == 0:
            # Run a constrained minimization for each geometry in the provided file.
            for i in range(len(self.M)):
                fy1, dnm = self.run_minimize(frame=i)
                self.optnext.append(([fy1],fy1))
                self.optinfo.append((fy1, dnm))
            print "Initializing, %i optimizations in this cycle" % (len(self.optnext))
            self.print_ascii_image()
        else:
            # These constrained minimizations should have been scheduled by read_minimize().
            print
            print "Iteration %i, %ix4 optimizations in the queue" % (self.iteration, len(self.optnext))
            self.print_ascii_image()
            if len(self.optnext) == 0:
                print "Scan is finished"
                return 0
            for fys, fy0 in self.optnext:
                for fy in fys:
                    fy1, dnm = self.run_minimize(fy, fy0)
                    self.optinfo.append((fy1, dnm))
        # Clear the list of optimizations that need to be done next, because we just launched them.
        self.optnext = []
        self.running = 1
        print
        return 1

    def read_iteration(self):
        """ Read the results of an iteration. """
        if self.running: 
            # Set the "running" variable to zero if all of the Work Queue tasks are complete.
            if any([T.name == self.aname for T in IQTasks]): return
            else: self.running = 0
        # Read the results of the minimization for every point that we did (stored in optinfo).
        # If not on iteration zero, also print out a job summary.
        if self.iteration == 0:
            for fy1, dnm in self.optinfo:
                self.read_minimize(fy1, dnm)
        else:
            self.lower = []
            self.new = []
            news = 0
            lowers = 0
            equals = 0
            fails = 0
            for fy1, dnm in self.optinfo:
                status = self.read_minimize(fy1, dnm)
                if status == 0:
                    news += 1
                    self.new.append(fy1)
                elif status == 1:
                    lowers += 1
                    if fy1 not in self.new:
                        self.lower.append(fy1)
                elif status == 2:
                    equals += 1
                elif status == -1:
                    fails += 1
            print "\nOptimization results for %s: %i new, %i lower, %i equal/higher, %i failed" % (self.aname, news, lowers, equals, fails)
        self.iteration += 1

    def finish(self):
        """ Print out the final energies to a file located in ALA/Scan/HF/scan.txt and also record the geometries. """
        mdnm = os.path.join(self.root, self.angles, self.method)
        xyzfin = []
        commfin = []
        o = open(os.path.join(mdnm, 'scan.txt'),'w')
        for i in self.energies.items():
            if i[1]:
                xyzfin.append(self.geoms[i[0]])
                commfin.append("Dihedral Angles = "+','.join(["%i" % j for j in i[0]]))
                print >> o, i[0][0], i[0][1], i[1]
        o.close()
        o = open(os.path.join(mdnm, 'xy.txt'),'w')
        for i in rng:
            print >> o, i
        o.close()
        o = open(os.path.join(mdnm, 'zz.txt'),'w')
        for i in self.energies.items():
            if i[1]:
                print >> o, i[1],
                if i[0][1] == 180:
                    print >> o
        o.close()
        Mfin = self.M
        Mfin.xyzs = xyzfin
        Mfin.comms = commfin
        Mfin.align(select=sorted(list(set([self.iA, self.iB, self.iC, self.iD, self.iE, self.iF, self.iG, self.iH]))))
        Mfin.write(os.path.join(mdnm, "scan.xyz"))
        # Mfin.write(os.path.join(mdnm, "scan.pdb"))

def main():
    print "Usage: %s method" % __file__
    print "where method is one of: %s" % (', '.join(MDict.keys()))
    
    method = sys.argv[1]
    if method not in MDict.keys():
        raise RuntimeError("method %s not implemented" % method)

    M = Molecule("propanol.xyz")
    d1 = [0, 1, 6, 10]
    d2 = [1, 6, 10, 11]

    DG = DihedralGrid(M, 0, 1, method, d1, d2, "propanol")

    while True:
        # Determine which new optimizations to start
        # and schedule them in the reservoir.
        if not DG.launch_iteration():
            DG.finish()
            break
        manage_wq()
        # Read the results for the sets of optimizations that are finished
        # and schedule new jobs if necessary.
        DG.read_iteration()
    print "Dihedral scan is finished!"

if __name__ == "__main__":
    main()
