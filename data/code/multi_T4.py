import time
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import mdtraj.reporters
import sys

cutoff = 0.95 * u.nanometers
output_frequency = 25000
n_steps = 500000000
temperature = 293. 
pressure = 1.0 * u.atmospheres

rank = int(sys.argv[1])
time.sleep(rank)  # This makes sure that no two jobs run at the same time for RNG purpuses.


pdb_filename = "./1am7_equil.pdb"
dcd_filename = "./Trajectories/1am7_%d.dcd" % rank
log_filename = "./Trajectories/1am7_%d.log" % rank
chk_filename = "./Trajectories/1am7_%d.chk" % rank

traj = mdtraj.load(pdb_filename)
top, bonds = traj.top.to_dataframe()
atom_indices = top.index[top.chainID == 0].values

pdb = app.PDBFile(pdb_filename)
topology = pdb.topology
positions = pdb.positions

ff = app.ForceField('amber99sbnmr.xml', 'tip3p-fb.xml')

platform = mm.Platform.getPlatformByName("CUDA")

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, 1.0 / u.picoseconds, 2.0 * u.femtoseconds)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))

simulation = app.Simulation(topology, system, integrator, platform=platform)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature)

print("Using platform %s" % simulation.context.getPlatform().getName())

dcd_file = open(dcd_filename, 'a')  # Try to append
simulation.reporters.append(mdtraj.reporters.DCDReporter(dcd_file, output_frequency, atomSubset=atom_indices))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), 5000, step=True, time=True, speed=True))
simulation.reporters.append(app.CheckPointReporter(open(chk_filename, 'w'), output_frequency, step=True, time=True, speed=True))
simulation.step(n_steps)
