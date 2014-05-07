from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import mdtraj.reporters

cutoff = 0.95 * u.nanometers

output_frequency = 25000
n_steps = 500000000
temperature = 293. 
pressure = 1.0 * u.atmospheres


pdb_filename = "./1am7_equil.pdb"
dcd_filename = "./1am7.dcd"
log_filename = "./1am7.log"

traj = mdtraj.load(pdb_filename)
top, bonds = traj.top.to_dataframe()
atom_indices = top.index[top.chainID == 0].values

pdb = app.PDBFile(pdb_filename)
topology = pdb.topology
positions = pdb.positions

ff = app.ForceField('amber99sbnmr.xml', 'tip3p-fb.xml')

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, 1.0 / u.picoseconds, 2.0 * u.femtoseconds)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))

simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature)

print("Using platform %s" % simulation.context.getPlatform().getName())

simulation.reporters.append(mdtraj.reporters.DCDReporter(dcd_filename, output_frequency, atomSubset=atom_indices))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), 1000, step=True, time=True))
simulation.step(n_steps)
