#!/usr/bin/env python
import itertools
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import sys
import mdtraj, mdtraj.reporters
import os

cutoff = 0.9 * u.nanometers
temperature = 298 * u.kelvin
timestep = 2.0 * u.femtoseconds
pressure = 1.0 * u.atmosphere
friction = 0.25 / u.picoseconds
n_steps = 100000000
output_frequency = 1000

ff_name = "amber99sbnmr"

ff = app.ForceField('%s.xml' % ff_name, 'tip3p.xml')

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]
amino_acids.remove("G")

sequences = ["%s-%s" % ("G", aa) for aa in amino_acids]
sequences.extend(["%s-%s" % (aa, "G") for aa in amino_acids]

capped_string = "capped"

platform = mm.Platform.getPlatformByName("CUDA")
platform.setPropertyDefaultValue('CudaDeviceIndex', os.environ.get("CUDA_VISIBLE_DEVICES"))

rank = int(sys.argv[1])

for sequence in sequences:
    
    if k != rank:
        continue
    
    pdb_filename = "./box/%s%s-%s.pdb" % (a0, a1, capped_string)
    dcd_filename = 'dcd/%s%s-%s.dcd' % (a0, a1, capped_string)
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
