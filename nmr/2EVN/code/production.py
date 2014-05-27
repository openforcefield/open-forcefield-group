import os
import time
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import mdtraj.reporters
import sys

code = "2evn"

which_forcefield = "amber99sbildn.xml"
which_water = 'tip3p-fb.xml'

platform_name = "CUDA"
timestep = 2.0 * u.femtoseconds
cutoff = 0.95 * u.nanometers
output_frequency = 25000
n_steps = 500000000
temperature = 300. 
pressure = 1.0 * u.atmospheres

rank = int(sys.argv[1])
time.sleep(rank)  # This makes sure that no two jobs run at the same time for RNG purpuses.


pdb_filename = "./%s_equil.pdb" % code
dcd_filename = "./Trajectories/%s_%d.dcd" % (code, rank)
log_filename = "./Trajectories/%s_%d.log" % (code, rank)

traj = mdtraj.load(pdb_filename)
top, bonds = traj.top.to_dataframe()
atom_indices = top.index[top.chainID == 0].values

pdb = app.PDBFile(pdb_filename)
topology = pdb.topology
positions = pdb.positions

ff = app.ForceField(which_forcefield, which_water)

platform = mm.Platform.getPlatformByName(platform_name)

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, 1.0 / u.picoseconds, timestep)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))

simulation = app.Simulation(topology, system, integrator, platform=platform)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature)


print("Using platform %s" % simulation.context.getPlatform().getName())

if os.path.exists(dcd_filename):
    sys.exit()

simulation.reporters.append(mdtraj.reporters.DCDReporter(dcd_filename, output_frequency, atomSubset=atom_indices))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True))
simulation.step(n_steps)
