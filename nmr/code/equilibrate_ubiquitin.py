from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

padding = 1.0 * u.nanometers
cutoff = 0.95 * u.nanometers

ff = app.ForceField('amber99sbnmr.xml', 'tip3p-fb.xml')

temperature = 298. 
pressure = 1.0 * u.atmospheres
ionicStrength = 185 * u.millimolar

pdb = app.PDBFile("./1d3z.pdb")

modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(ff, padding=padding, ionicStrength=ionicStrength)

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
simulation.step(25000)
state = simulation.context.getState(getPositions=True, getParameters=True)

positions = state.getPositions()
vectors = state.getPeriodicBoxVectors()
vectors = tuple([v[i] / u.nanometers for (i,v) in enumerate(vectors)])
vectors = u.Quantity(vectors, u.nanometer)
topology.setUnitCellDimensions(vectors)
app.PDBFile.writeFile(topology, positions, open("./1d3z_equil.pdb", 'w'))
