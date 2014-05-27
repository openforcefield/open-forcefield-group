import itertools
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

padding = 1.0 * u.nanometers
cutoff = 0.9 * u.nanometers

ff = app.ForceField('amber99sbnmr.xml', 'tip3p.xml')

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]
amino_acids_nogly = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "A", "I", "L", "M", "F", "W", "Y", "V"]

sequences = ["%s-%s-%s" % ("ACE", aa, "NME") for aa in amino_acids]
sequences = []
sequences.extend(["%s-%s%s-%s" % ("ACE", "G", aa, "NH2") for aa in amino_acids_nogly])
sequences.extend(["%s-%s%s-%s" % ("ACE", aa, "G", "NH2") for aa in amino_acids_nogly])

temperatures_by_seq_length = {
1: 303. * u.kelvin,  # dipeptide X
2: 298. * u.kelvin,
}  # dipeptide X-Y
pressure = 1.0 * u.atmospheres

for sequence in sequences:
    temperature = temperatures_by_seq_length[len(sequences[25].split("-")[1])]  # Lookup the correct temperature.  
    base_filename = "%s.pdb" % sequence
    print(base_filename)
    
    pdb = app.PDBFile("./pdbs/%s" % base_filename)
    
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(ff, model="tip3p", padding=padding)
    topology = modeller.getTopology()
    system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)
    
    integrator = mm.LangevinIntegrator(temperature, 1.0/u.picoseconds, 1.0*u.femtoseconds)
    system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(modeller.getPositions())
    print('Minimizing...')
    simulation.minimizeEnergy()

    simulation.context.setVelocitiesToTemperature(temperature)
    print('Equilibrating...')
    simulation.step(15000)
    state = simulation.context.getState(getPositions=True, getParameters=True)
    
    positions = state.getPositions()
    vectors = state.getPeriodicBoxVectors()
    vectors = tuple([v[i] / u.nanometers for (i,v) in enumerate(vectors)])
    vectors = u.Quantity(vectors, u.nanometer)
    topology.setUnitCellDimensions(vectors)
    app.PDBFile.writeFile(topology, positions, open("./boxes/%s" % base_filename, 'w'))
