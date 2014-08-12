from rdkit.Chem import PandasTools
from rdkit import Chem
from rdkit.Chem import AllChem
import cirpy
import pandas as pd
import bs4

def parse_page(soup):
    for x in soup.find_all("td"):
        t = x.get("class")
        if t is not None and "term2TD" in t:
            if "CAS No." in x.next:
                a0 = x
                a1 = x.next
                a2 = x.next.next
                cas = a2.text[1:]
            if "Density" in x.next:
                a0 = x
                a1 = x.next
                a2 = x.next.next
                density = a2.text
            if "TDENL" in x.next:
                a0 = x
                a1 = x.next
                a2 = x.next.next
                temperature = a2.text
            if "Molecular Wt." in x.next:
                a0 = x
                a1 = x.next
                a2 = x.next.next
                weight = a2.text
    if density == "NA":
        density = None
    smiles = cirpy.resolve(cas, "smiles")
    return (cas, density, temperature, weight, smiles)

data = []
for i in range(1, 2000):
    f = open("./pages/page%d.html" % i).read()
    soup = bs4.BeautifulSoup(f)
    CAS, density, temperature, weight, smiles = parse_page(soup)
    data.append([i, CAS, density, temperature, weight, smiles])

data = pd.DataFrame(data, columns=["Index", "CAS", "MoleDensity", "Temperature", "weight", "smiles"])
data.set_index("Index")
data.dropna(inplace=True)
data.MoleDensity = data.MoleDensity.astype('float')
data.weight = data.weight.astype('float')
data.Temperature = data.Temperature.astype('float')
data["density"] = data.weight * data.MoleDensity
PandasTools.AddMoleculeColumnToFrame(data, smilesCol='smiles', molCol='molecule', includeFingerprints=False)
#	g-mol/cm^3 units

has_other = {}
num_atoms = {}
for k, m in data.molecule.iteritems():
    try:
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m)
        AllChem.UFFOptimizeMolecule(m)
        atoms = pd.Series([a.GetSymbol() for a in m.GetAtoms()])
        atom_counts = atoms.value_counts()
        has_other[k] = len(np.setdiff1d(atom_counts.index.values, ["C", "N", "H", "O"])) > 0
        data["molecule"][k] = m
        num_atoms[k] = m.GetNumAtoms()
    except:
        pass

data["has_other"] = pd.Series(has_other)
data["num_atoms"] = pd.Series(num_atoms)

data.to_hdf("./alldata.h5", "data")

data = data[data.has_other == False]
data = data[data.num_atoms <= 25]
data = data[(data.Temperature >= 273) & (data.Temperature <= 320)]
data.to_hdf("./data.h5", "data")
