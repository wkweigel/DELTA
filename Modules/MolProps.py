import numpy as np 
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors


class Properties():
    '''Calculates stores molecular properties from an input csv file

    Input:
    ==========
    inputfile: Path to a csv file containing a least the following three columns: 'Smiles', 'DEL', 'ID'
        


    '''
    def __init__(self,inputfile):
        self.data_df = pd.read_csv(inputfile)
        self.DEL_list=self.data_df.DEL.unique()

        #Create the columns to be used for fingerprinting 
        self.data_df['ROMol'] = self.data_df.Smiles.apply(Chem.MolFromSmiles)
        self.data_df['Smiles'] = self.data_df.ROMol.apply(lambda x: Chem.MolToSmiles(x, kekuleSmiles=True, isomericSmiles=False))
        self.Smiles=self.data_df.Smiles.to_list()
        self.Del=self.data_df.DEL.to_list()
        self.IDs=self.data_df.ID.to_list()
        self.property_dict={
            'DEL':self.Del,
            'ID':self.IDs,
            'smiles':self.Smiles,
            'mol_wt':[], 
            'logp':[], 
            'h_donors':[],
            'h_acceptors':[],
            'rotatable_bonds':[],
            'polar_surface_area':[],
            'atoms':[],
            'heavy_atoms':[],
            'rings':[] }
        
        for mol in self.Smiles:
            #mol=rdMolStandardize.StandardizeSmiles(mol)
            molecule=Chem.MolFromSmiles(mol)
            
            self.property_dict['mol_wt'].append(Descriptors.ExactMolWt(molecule))
            self.property_dict['logp'].append(Descriptors.MolLogP(molecule))
            self.property_dict['h_donors'].append(Descriptors.NumHDonors(molecule))
            self.property_dict['h_acceptors'].append(Descriptors.NumHAcceptors(molecule))
            self.property_dict['rotatable_bonds'].append(Descriptors.NumRotatableBonds(molecule))
            self.property_dict['polar_surface_area'].append(Chem.QED.properties(molecule).PSA)
            self.property_dict['atoms'].append(Chem.rdchem.Mol.GetNumAtoms(molecule))
            self.property_dict['heavy_atoms'].append(Chem.rdchem.Mol.GetNumHeavyAtoms(molecule))
            self.property_dict['rings'].append(Chem.rdMolDescriptors.CalcNumRings(molecule))

        self.property_df=pd.DataFrame(self.property_dict)
        self.property_df.to_csv('Notebook Outputs/Properties.csv')
    
    def Check_Lipinski_Verber_Params(self):
        Lipinski_counter=0
        Lipinski_params={'mol_wt':0, 'logp':0 ,'h_donors':0,'h_acceptors':0,'rotatable_bonds':0,'polar_surface_area':0, 'All':0}
        passing_IDs=[]
        passing_smiles=[]
        passing_DEL=[]
        druglike_df=pd.DataFrame()
        for index, row in self.property_df.iterrows():
            Lipinski_counter=0
            if row['mol_wt'] <= 500:
                Lipinski_params['mol_wt']+=1
                Lipinski_counter+=1
            if row['logp']<= 5:
                Lipinski_params['logp']+=1
                Lipinski_counter+=1
            if row['h_donors'] <= 5:
                Lipinski_params['h_donors']+=1
                Lipinski_counter+=1
            if row['h_acceptors'] <= 10:
                Lipinski_params['h_acceptors']+=1
                Lipinski_counter+=1
            if row['rotatable_bonds'] <= 5:
                Lipinski_params['rotatable_bonds']+=1
                Lipinski_counter+=1
            if row['polar_surface_area'] <=140:
                Lipinski_params['polar_surface_area']+=1
                Lipinski_counter+=1
            if Lipinski_counter==6:
                Lipinski_params['All']+=1
                passing_DEL.append(row['DEL'])
                passing_IDs.append(row['ID'])
                passing_smiles.append(row['smiles'])
        druglike_df['DEL']=passing_DEL    
        druglike_df['ID']=passing_IDs
        druglike_df['Smiles']=passing_smiles
        return(Lipinski_params,druglike_df)
    
    def assemble_plot_df(self,property):
        output_df=pd.DataFrame()
        for DELtype in self.DEL_list:
            temp_df=self.property_df.loc[self.property_df['DEL'] == DELtype]
            property_vals=np.array(temp_df[property].values.tolist())
            output_df[DELtype]=property_vals
        return(output_df)
        
   





    def Property_Dataframe(self,property): 
        '''Generate a dataframe for a specfic property column in an input csv file'''
        self.prop_df=pd.DataFrame() #create an empty df for holding the properties
        Property=[]
        #For each type in the csv "DEL" column:
        for Del in self.DEL_list:
            DEL_df=pd.DataFrame()
            temp_df=self.data_df.loc[self.data_df['DEL'] == str(Del)] #create a temp dataframe for the current 'Del' iteration
            temp_properties=np.array(temp_df[str(property)].values.tolist())  #extract the property values as np array
            DEL_df[Del]=temp_properties #create the DEL_df for the current Del iteration
            self.prop_df=pd.concat([self.prop_df, DEL_df], axis=1) #Concatonate the DEL_df to the main prop_df
        return(self.prop_df)

 