from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
from bpti_parameters import *

which_forcefield = "amber99sbildn.xml"
which_water = 'tip3p-fb.xml'

ff = app.ForceField(which_forcefield, which_water)
pdb = app.PDBFile("./%s_fixed.pdb" % code)
output_pdb = "./equil.pdb"

modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(ff, padding=padding, ionicStrength=ionicStrength)

topology = modeller.topology
positions = modeller.positions

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, friction, equilibration_timestep)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostat_frequency))
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(temperature)
print('Running.')
simulation.reporters.append(app.PDBReporter(output_pdb, n_equil_steps - 1))
simulation.step(n_equil_steps)
