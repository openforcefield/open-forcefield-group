#!/usr/bin/env python
import itertools
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import sys
import mdtraj, mdtraj.reporters
import os
from dipeptide_parameters import *

platform = mm.Platform.getPlatformByName("CUDA")
platform.setPropertyDefaultValue('CudaDeviceIndex', os.environ.get("CUDA_VISIBLE_DEVICES"))

rank = int(sys.argv[1])

for aa in amino_acids:
    
    if k != rank:
        continue

    ff = app.ForceField('%s.xml' % ff_name, '%s.xml' % water_name)
    pdb_filename = "./box/%s-%s.pdb" % (aa, capped_string)
    dcd_filename = 'dcd/%s-%s.dcd' % (a0, capped_string)
    print(k)
    print(pdb_filename)

    if os.path.exists(dcd_filename):
        continue    

    traj = mdtraj.load(pdb_filename)
    top, bonds = traj.top.to_dataframe()
    atom_indices = top.index[top.chainID == 0].values

    pdb = app.PDBFile(pdb_filename)
    
    system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(temperature, friction, timestep)
    system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))

    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)

    simulation.context.setVelocitiesToTemperature(300 * u.kelvin)
    print('Running.')
    simulation.reporters.append(mdtraj.reporters.DCDReporter(dcd_filename, output_frequency, atomSubset=atom_indices))
    simulation.step(n_steps)
