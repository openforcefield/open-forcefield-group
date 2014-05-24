from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import pdbfixer

padding = 1.0 * u.nanometers
cutoff = 0.95 * u.nanometers

ff = app.ForceField('amber99sbnmr.xml', 'tip3p-fb.xml')

temperature = 293. 
pressure = 1.0 * u.atmospheres


fixer = pdbfixer.PDBFixer("./1am7.pdb")

fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.removeHeterogens(True)
fixer.addMissingHydrogens()
fixer.removeChains([1, 2, 3, 4, 5])
app.PDBFile.writeFile(fixer.topology, fixer.positions, open("1am7_fixed.pdb", 'w'))
