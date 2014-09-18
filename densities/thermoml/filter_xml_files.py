import glob
from lxml import etree

def has_density(root):
    citation = root.xpath("./*[local-name()='Citation']")[0]
    keywords = citation.xpath("./*[local-name()='sKeyword']")
    if any(["Density" in keyword.text for keyword in keywords]):
        return True
    else:
        return False
        
def get_compounds(root):
    compounds = root.xpath("./*[local-name()='Compound']")
    all_data = []
    for compound in compounds:
        regnum = compound.xpath("./*[local-name()='RegNum']")[0]
        nOrgNum = int(regnum.xpath("./*[local-name()='nOrgNum']")[0].text)
        sCommonName = compound.xpath("./*[local-name()='sCommonName']")[0].text
        data = {}
        data["nOrgNum"] = nOrgNum
        data["sCommonName"] = sCommonName
        all_data.append(data)
    return all_data

def get_density_property(root):
    PureOrMixtureData = root.xpath("./*[local-name()='PureOrMixtureData']")


filenames = glob.glob("./*/*.xml")
good_filenames = []
compounds = {}
for filename in filenames:
    try:
        tree = etree.parse(filename)
    except IOError:
        continue
    root = tree.getroot()
    if not has_density(root):
        continue
    print(filename)
    good_filenames.append(filename)
    compounds[filename] = get_compounds(root)


filename = good_filenames[-1]
tree = etree.parse(filename)
root = tree.getroot()
PureOrMixtureData = root.xpath("./*[local-name()='PureOrMixtureData']")


