#!/usr/bin/env python
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import sys
import os
from dipeptide_parameters import *

#platform = mm.Platform.getPlatformByName("CUDA")
#platform.setPropertyDefaultValue('CudaDeviceIndex', os.environ.get("CUDA_VISIBLE_DEVICES"))

rank = int(sys.argv[1])

for k, (ff_name, water_name, seq) in enumerate(products):
    
    if k != rank:
        continue

    pdb_filename = "./pdbs/%s.pdb" % (seq)
    output_pdb = './boxes/%s_%s_%s.pdb' % (ff_name, water_name, seq)
    print(k)
    print(pdb_filename)

    if os.path.exists(output_pdb):
        continue    

    ff = app.ForceField('%s.xml' % ff_name, '%s.xml' % water_name)
    
    pdb = app.PDBFile(pdb_filename)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(ff, model=base_waters[water_name], padding=padding)
    topology = modeller.getTopology()
    positions = modeller.getPositions()
    
    system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(temperature, friction, timestep)
    system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)

    simulation.context.setVelocitiesToTemperature(temperature)
    print('Running.')
    simulation.reporters.append(app.PDBReporter(output_pdb, equilibrate_output_frequency))
    simulation.step(n_steps)



    
