from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

code = "2evn"

padding = 1.0 * u.nanometers
cutoff = 0.95 * u.nanometers

ff = app.ForceField('amber99sbildn.xml', 'tip3p-fb.xml')

temperature = 300. 
pressure = 1.0 * u.atmospheres

pdb = app.PDBFile("./%s_fixed.pdb" % code)

modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(ff, padding=padding)

topology = modeller.topology
positions = modeller.positions

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, 1.0 / u.picoseconds, 1.0 * u.femtoseconds)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(temperature)
print('Equilibrating...')
simulation.step(50000)
state = simulation.context.getState(getPositions=True, getParameters=True)

positions = state.getPositions()
vectors = state.getPeriodicBoxVectors()
vectors = tuple([v[i] / u.nanometers for (i,v) in enumerate(vectors)])
vectors = u.Quantity(vectors, u.nanometer)
topology.setUnitCellDimensions(vectors)
app.PDBFile.writeFile(topology, positions, open("./%s__equil.pdb" % code, 'w'))
