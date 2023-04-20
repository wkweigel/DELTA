import numpy as np 
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors

data_df = pd.read_csv(inputfile)

#Create a list of the unique members in the "MolType" column
DEL_list=data_df.MolType.unique()

def Property_Dataframe(property): 
    '''Generate a dataframe for a specfic property column in an input csv file'''
    prop_df=pd.DataFrame() #create an empty df for holding the properties
    Property=[]
    #For each type in the csv "MolType" column:
    for Del in DEL_list:
        DEL_df=pd.DataFrame()
        temp_df=data_df.loc[data_df['MolType'] == str(Del)] #create a temp dataframe for the current 'Del' iteration
        temp_properties=np.array(temp_df[str(property)].values.tolist())  #extract the property values as np array
        DEL_df[Del]=temp_properties #create the DEL_df for the current Del iteration
        prop_df=pd.concat([prop_df, DEL_df], axis=1) #Concatonate the DEL_df to the main prop_df
    return(prop_df)

def Create_Property_df():
    smiles=data_df.smiles.to_list()
    Del=data_df.DEL.to_list()
    property_dict={
        'DEL':Del,
        'smiles':smiles,
        'mol_wt':[], 
        'logp':[], 
        'h_donors':[],
        'h_acceptors':[],
        'rotatable_bonds':[],
        'atoms':[],
        'polar_surface_area':[],
        'heavy_atoms':[],
        'rings':[] }
    
    for mol in smiles:
        molecule=Chem.MolFromSmarts(mol)
        property_dict['mol_wt'].append(Descriptors.ExactMolWt(molecule))
        property_dict['logp'].append(Descriptors.MolLogP(molecule))
        property_dict['h_donors'].append(Descriptors.NumHDonors(molecule))
        property_dict['h_acceptors'].append(Descriptors.NumHAcceptors(molecule))
        property_dict['rotatable_bonds'].append(Descriptors.NumRotatableBonds(molecule))
        property_dict['atoms'].append(Chem.rdchem.Mol.GetNumAtoms(molecule))
        property_dict['polar_surface_area'].append(Chem.QED.properties(molecule).PSA)
        property_dict['heavy_atoms'].append(Chem.rdchem.Mol.GetNumHeavyAtoms(molecule))
        property_dict['rings'].append(Chem.rdMolDescriptors.CalcNumRings(molecule))

    Property_df=pd.DataFrame(property_dict)
    Property_df.to_csv('Properties.csv')
    return(Property_df)
