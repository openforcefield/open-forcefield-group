from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
import pdbfixer

code = "2evn"

fixer = pdbfixer.PDBFixer("./%s.pdb" % code)
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.removeHeterogens(True)
fixer.addMissingHydrogens()
#fixer.removeChains([1, 2, 3, 4, 5])
app.PDBFile.writeFile(fixer.topology, fixer.positions, open("%s_fixed.pdb" % code, 'w'))
