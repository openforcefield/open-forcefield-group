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

for k, (ff_name, water_name, seq) in enumerate(products):
    
    if k != rank:
        continue

    pdb_filename = "./boxes/%s_%s_%s.pdb" % (ff_name, water_name, seq)
    dcd_filename = './dcd/%s_%s_%s.dcd' % (ff_name, water_name, seq)
    print(k)
    print(pdb_filename)

    if os.path.exists(dcd_filename):
        continue    

    ff = app.ForceField('%s.xml' % ff_name, '%s.xml' % water_name)
    
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(ff, model=base_waters[water_name], padding=padding)
    topology = modeller.getTopology()
    
    
    traj = mdtraj.load(pdb_filename)
    top, bonds = traj.top.to_dataframe()
    atom_indices = top.index[top.chainID == 0].values

    pdb = app.PDBFile(pdb_filename)
    
    system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(temperature, friction, timestep)
    system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)

    simulation.context.setVelocitiesToTemperature(temperature)
    print('Running.')
    simulation.reporters.append(app.DCDReporter(dcd_filename, equilibrate_output_frequency))
    simulation.step(n_steps)



    
