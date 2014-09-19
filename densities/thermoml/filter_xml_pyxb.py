import pandas as pd
import glob
import binding  # Obtained by `wget http://media.iupac.org/namespaces/ThermoML/ThermoML.xsd` and `pyxbgen ThermoML.xsd`

filename = "./10.1007/s10765-010-0742-8.xml"
root = binding.CreateFromDocument(open(filename).read())

compound_dict = {}
for Compound in root.Compound:
    nOrgNum = Compound.RegNum.nOrgNum
    sCommonName = Compound.sCommonName[0]
    sFormulaMolec = Compound.sFormulaMolec
    compound_dict[nOrgNum] = dict(sCommonName=sCommonName, sFormulaMolec=sFormulaMolec)
print(compound_dict)

for PureOrMixtureData in root.PureOrMixtureData:
    print("PureOrMixtureData")
    component_list = []
    for Component in PureOrMixtureData.Component:
        print("Component")
        nSampleNm = Component.nSampleNm
        nOrgNum = Component.RegNum.nOrgNum
        sCommonName = compound_dict[nOrgNum]
        component_list.append(dict(sCommonName=sCommonName, nOrgNum=nOrgNum))
    print(component_list)

    property_dict = {}
    for Property in PureOrMixtureData.Property:
        print("Property")
        nPropNumber = Property.nPropNumber
        try:
            ePropName = Property.Property_MethodID.PropertyGroup.VolumetricProp.ePropName
        except AttributeError:
            ePropName = None
        property_dict[nPropNumber] = dict(ePropName=ePropName)
    print(property_dict)
    
    pressure, temperature, solvent_composition = None, None, {}
    for Constraint in PureOrMixtureData.Constraint:
        print("Constraint")
        nConstraintValue = Constraint.nConstraintValue
        ConstraintType = Constraint.ConstraintID.ConstraintType
        if ConstraintType.ePressure is not None:
            pressure = nConstraintValue
        if ConstraintType.eTemperature is not None:
            temperature = nConstraintValue
        if ConstraintType.eSolventComposition is not None:
            nOrgNum = Constraint.ConstraintID.RegNum.nOrgNum
            sCommonName = compound_dict[nOrgNum]["sCommonName"]
            solvent_composition[sCommonName] = nConstraintValue
    print(pressure, temperature, solvent_composition)

    for Variable in PureOrMixtureData.Variable:
        nVarNumber = Variable.nVarNumber
        
