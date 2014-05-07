from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import pdbfixer

padding = 1.0 * u.nanometers
cutoff = 0.95 * u.nanometers

ff = app.ForceField('amber99sbnmr.xml', 'tip3p-fb.xml')

temperature = 293. 
pressure = 1.0 * u.atmospheres


fixer = pdbfixer.PDBFixer("./1am7_chain0_renamed.pdb")
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingHydrogens()
fixer.addMissingAtoms()
fixer.addSolvent(padding=padding)

topology = fixer.topology

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, 1.0 / u.picoseconds, 1.0 * u.femtoseconds)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(fixer.positions)
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
app.PDBFile.writeFile(topology, positions, open("./boxes/%s" % base_filename, 'w'))
