from rdkit import Chem
from rdkit.Chem import AllChem
import cirpy
import pandas as pd
import bs4

data = pd.read_hdf("./data.h5", "data")

molecules = {}
has_other = {}
num_atoms = {}
max_charge = {}
for k, smiles in data.smiles.iteritems():
    try:
        m = Chem.MolFromSmiles(smiles)
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m)
        AllChem.UFFOptimizeMolecule(m)
        atoms = pd.Series([a.GetSymbol() for a in m.GetAtoms()])
        atom_counts = atoms.value_counts()
        has_other[k] = len(np.setdiff1d(atom_counts.index.values, ["C", "N", "H", "O"])) > 0
        molecules[k] = m
        num_atoms[k] = m.GetNumAtoms()
        max_charge[k] = max([abs(a.GetFormalCharge()) for a in m.GetAtoms()])
    except:
        pass

data["has_other"] = pd.Series(has_other)
data["molecules"] = pd.Series(molecules)
data["num_atoms"] = pd.Series(num_atoms)
data["max_charge"] = pd.Series(max_charge)

data = data[data.has_other == False]
data = data[data.num_atoms <= 25]
data = data[(data.Temperature >= 273) & (data.Temperature <= 320)]
data = data[data.max_charge == 0]

data.to_hdf("./subset.h5", "data")
