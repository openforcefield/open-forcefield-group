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

    alldata = []
    for PureOrMixtureData in root.PureOrMixtureData:
        component_list = []
        for Component in PureOrMixtureData.Component:
            nSampleNm = Component.nSampleNm
            nOrgNum = Component.RegNum.nOrgNum
            sCommonName = compound_dict[nOrgNum]
            component_list.append(dict(sCommonName=sCommonName, nOrgNum=nOrgNum))

        property_dict = {}
        for Property in PureOrMixtureData.Property:
            nPropNumber = Property.nPropNumber
            ePropName = Property.Property_MethodID.PropertyGroup.content()[0].ePropName  # ASSUMING LENGTH 1
            property_dict[nPropNumber] = ePropName

        state = dict(filename=filename)
        
        state["Pressure, kPa"] = None  # This is the only pressure unit used in ThermoML
        state['Temperature, K'] = None  # This is the only temperature unit used in ThermoML
        
        composition = dict()
        for Constraint in PureOrMixtureData.Constraint:
            nConstraintValue = Constraint.nConstraintValue
            ConstraintType = Constraint.ConstraintID.ConstraintType
            
            assert len(ConstraintType.content()) == 1
            constraint_type = ConstraintType.content()[0]
            state[constraint_type] = nConstraintValue
            print(filename, constraint_type)
            if ConstraintType.eSolventComposition is not None:
                nOrgNum = Constraint.ConstraintID.RegNum.nOrgNum
                sCommonName = compound_dict[nOrgNum]["sCommonName"]
                solvents = [compound_dict[x.nOrgNum]["sCommonName"] for x in Constraint.Solvent.RegNum]
            
            if constraint_type in ["Mole fraction", "Mass Fraction", "Molality, mol/kg", "Solvent: Amount concentration (molarity), mol/dm3"]:
                nOrgNum = Constraint.ConstraintID.RegNum.nOrgNum
                sCommonName = compound_dict[nOrgNum]["sCommonName"]
                solvents = [compound_dict[x.nOrgNum]["sCommonName"] for x in Constraint.Solvent.RegNum]
                solvent_string = "%s___%s" % (sCommonName, "__".join(solvents))
                state["%s metadata" % constraint_type] = solvent_string
                print(solvent_string)

        variable_dict = {}
        for Variable in PureOrMixtureData.Variable:
            nVarNumber = Variable.nVarNumber
            VariableType = Variable.VariableID.VariableType
            assert len(VariableType.content()) == 1
            vtype = VariableType.content()[0]  # Assume length 1, haven't found counterexample yet.
            variable_dict[nVarNumber] = vtype
            if vtype in ["Mole fraction", "Mass Fraction", "Molality, mol/kg", "Solvent: Amount concentration (molarity), mol/dm3"]:
                nOrgNum = Variable.VariableID.RegNum.nOrgNum
                sCommonName = compound_dict[nOrgNum]["sCommonName"]
                if Variable.Solvent is not None:
                    solvents = [compound_dict[x.nOrgNum]["sCommonName"] for x in Variable.Solvent.RegNum]
                else:
                    solvents = []
                solvent_string = "%s___%s" % (sCommonName, "__".join(solvents))
                state["%s Variable metadata" % vtype] = solvent_string
                print(solvent_string)
        
        
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

            alldata.append(current_data)
    return alldata, root

filename = "./10.1007/s10765-010-0742-8.xml"
#alldata, root = parse(filename)


data = []
for filename in glob.glob("./*/*.xml")[40:70]:
    try:
        alldata, root = parse(filename)
    except IOError:
        continue
    for d in alldata:
        if True or u'Mass density, kg/m3' in d:
            data.append(d)

data = pd.DataFrame(data)
