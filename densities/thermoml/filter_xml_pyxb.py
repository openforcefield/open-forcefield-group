import copy
import pandas as pd
import glob
import thermoml_schema  # Obtained by `wget http://media.iupac.org/namespaces/ThermoML/ThermoML.xsd` and `pyxbgen ThermoML.xsd`

def parse(filename):
    root = thermoml_schema.CreateFromDocument(open(filename).read())

    compound_dict = {}
    for Compound in root.Compound:
        nOrgNum = Compound.RegNum.nOrgNum
        sCommonName = Compound.sCommonName[0]
        sFormulaMolec = Compound.sFormulaMolec
        compound_dict[nOrgNum] = dict(sCommonName=sCommonName, sFormulaMolec=sFormulaMolec)
    print(compound_dict)

    alldata = []
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
            ePropName = Property.Property_MethodID.PropertyGroup.content()[0].ePropName  # ASSUMING LENGTH 1
            property_dict[nPropNumber] = ePropName
        print(property_dict)

        #state = dict(pressure=None, temperature=None, solvent_composition={})
        state = dict(Pressure=None, Temperature=None)
        composition = dict()
        for Constraint in PureOrMixtureData.Constraint:
            print("Constraint")
            nConstraintValue = Constraint.nConstraintValue
            ConstraintType = Constraint.ConstraintID.ConstraintType
            if ConstraintType.ePressure is not None:
                state["Pressure"] = nConstraintValue
            if ConstraintType.eTemperature is not None:
                state["Temperature"] = nConstraintValue
            if ConstraintType.eSolventComposition is not None:
                nOrgNum = Constraint.ConstraintID.RegNum.nOrgNum
                sCommonName = compound_dict[nOrgNum]["sCommonName"]
                composition[sCommonName] = nConstraintValue
                solvents = [compound_dict[x.nOrgNum]["sCommonName"] for x in Constraint.Solvent.RegNum]
                composition["other"] = solvents            

        print(state)

        variable_dict = {}
        for Variable in PureOrMixtureData.Variable:
            nVarNumber = Variable.nVarNumber
            VariableType = Variable.VariableID.VariableType
            vtype, vunit = VariableType.content()[0].split(", ")  # Assuming length one!!!
            variable_dict[nVarNumber] = vtype
        
        
        for NumValues in PureOrMixtureData.NumValues:
            current_data = copy.deepcopy(state)  # Copy in values of constraints.
            current_composition = copy.deepcopy(composition)
            for VariableValue in NumValues.VariableValue:
                nVarValue = VariableValue.nVarValue
                nVarNumber = VariableValue.nVarNumber
                vtype = variable_dict[nVarNumber]
                current_data[vtype] = nVarValue

            for PropertyValue in NumValues.PropertyValue:
                nPropNumber = PropertyValue.nPropNumber
                nPropValue = PropertyValue.nPropValue
                ptype = property_dict[nPropNumber]
                current_data[ptype] = nPropValue
            current_data["composition"] = current_composition
            alldata.append(current_data)
    return alldata, root

filename = "./10.1007/s10765-010-0742-8.xml"
alldata, root = parse(filename)

data = []
for d in alldata:
    if u'Mass density, kg/m3' in d:
        print(d)
        data.append(d)

data = pd.DataFrame(data)



filename = "./10.1007/s10765-010-0742-8.xml"
data = []
for filename in glob.glob("./*/*.xml"):
    try:
        alldata, root = parse(filename)
    except IOError:
        continue
    for d in alldata:
        if u'Mass density, kg/m3' in d:
            print(d)
            data.append(d)

data = pd.DataFrame(data)
