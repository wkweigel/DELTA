import numpy as np
import altair as alt
from vega_datasets import data
import pandas as pd
from ugtm import pcaPreprocess
from mhfp.encoder import MHFPEncoder
from mxfp import mxfp
#from map4 import MAP4Calculator

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.AtomPairs import Pairs

from tqdm import tqdm
tqdm.pandas()

alt.data_transformers.disable_max_rows()

class PCA():
    '''Performs PCA analysis using diffrent molecular descriptors

    Input:
    ==========
    datafile: Path to a csv file containing a least the following three columns: 'Smiles', 'DEL', 'ID'
        
    descriptor: String specifying the fingerprinting method to use. Possible methods are: 'MXFP', 'MHFP6', 'ECFP6', 'MACCS', 'MQN', "AP"

    nBits: Number of bits to use. Only used for some descriptors

    '''
    def __init__(self,datafile,descriptor,nBits):
        self.nBits=nBits
        self.descriptor=descriptor
        self.Pproc_df = pd.read_csv(datafile) #Create pre-processing df from datafile
        self.Pproc_df['ROMol'] = self.Pproc_df.Smiles.apply(Chem.MolFromSmiles) #Create mols
        self.Pproc_df['Smiles'] = self.Pproc_df.ROMol.apply(lambda x: Chem.MolToSmiles(x, kekuleSmiles=True, isomericSmiles=False)) #Cleanup the smiles
        
        self.chart_df=pd.DataFrame()
        self.chart_df['Smiles'] = np.array(self.Pproc_df['Smiles'])
        self.chart_df['DEL'] = np.array(self.Pproc_df['DEL'])
        self.chart_df['ID']= np.array(self.Pproc_df['ID'])

        
        
        


        if self.descriptor=='MXFP':
            MXFP = mxfp.MXFPCalculator(dimensionality='2D')
            fp = [MXFP.mxfp_from_mol(x) for x in tqdm(self.Pproc_df['ROMol'],desc="Fingerprinting")]
            self.fp_array=np.array(fp)
            self.nBits=217

        if self.descriptor=='MHFP6':
            MHFP6 = MHFPEncoder(n_permutations=self.nBits)
            fp=[MHFP6.encode(x) for x in tqdm(self.Pproc_df['Smiles'],desc="Fingerprinting")]
            self.fp_array=np.array(fp)

        if self.descriptor=='ECFP6':
            ECFP6 = [AllChem.GetMorganFingerprintAsBitVect(x,radius=3, nBits=self.nBits) for x in tqdm(self.Pproc_df['ROMol'],desc="Fingerprinting")]
            fp_list=[list(l) for l in ECFP6]
            self.fp_array = np.array(fp_list)
        
        if self.descriptor=='ECFP4':
            ECFP6 = [AllChem.GetMorganFingerprintAsBitVect(x,radius=2, nBits=self.nBits) for x in tqdm(self.Pproc_df['ROMol'],desc="Fingerprinting")]
            fp_list=[list(l) for l in ECFP6]
            self.fp_array = np.array(fp_list)

        if self.descriptor=='MACCS':
            fp_list=[]
            MACCS=[rdMolDescriptors.GetMACCSKeysFingerprint(x) for x in tqdm(self.Pproc_df['ROMol'],desc="Fingerprinting")]
            fp_str=[fp.ToBitString() for fp in MACCS]
            for fp in fp_str:
                temp_list=[int(x) for x in fp]
                fp_list.append(temp_list)
            self.fp_array = np.array(fp_list)
            self.nBits=167

        if self.descriptor=='MQN':
            fp_list=[]
            MQN=[rdMolDescriptors.MQNs_(x)  for x in tqdm(self.Pproc_df['ROMol'],desc="Fingerprinting")]
            self.fp_array=np.array(MQN)
            self.nBits=42

        if self.descriptor=='MAP4':
            fp_list=[]
            MAP4= [MAP4Calculator(dimensions=self.nBits).calculate(x) for x in tqdm(self.Pproc_df['ROMol'],desc="Fingerprinting")]
            self.fp_array=np.array(MAP4)


        if self.descriptor=='AP':
            AP=[DataStructs.cDataStructs.FoldFingerprint(Pairs.GetAtomPairFingerprintAsBitVect(x),foldFactor=4096) for x in tqdm(self.Pproc_df['ROMol'],desc="Fingerprinting")]
            self.fp_array=np.array(AP)
            self.nBits=2048


    def make_lists_from_2D_coords(self):
        '''Takes a 2D array of coordinates and seperates it into two 1D lists of floating point values

        Parameters:
        ==========
        coord_array: An ndarray of coordinates to be split into seperate lists
        '''

        split_array=np.hsplit(self.pca_coords, self.nBits)
        PCAx1=split_array[0].tolist()
        PCAx2=split_array[1].tolist()
        PCAx1_list=[]
        PCAx2_list=[]
        for value in PCAx1:
            val_str=str(value)[1:-1]
            PCAx1_list.append(float(val_str))
        for value in PCAx2:
            val_str=str(value)[1:-1]
            PCAx2_list.append(float(val_str))
        
        return(PCAx1_list, PCAx2_list)


    def make_lists_from_3D_coords(self):
        '''Takes a 3D array of coordinates and seperates it into three 1D lists of floating point values

        Parameters:
        ==========
        coord_array: An ndarray of coordinates to be split into seperate lists
        '''
        
        split_array=np.hsplit(self.pca_coords, self.nBits)
        PCAx1=split_array[0].tolist()
        PCAx2=split_array[1].tolist()
        PCAx3=split_array[2].tolist()
        PCAx1_list=[]
        PCAx2_list=[]
        PCAx3_list=[]
        for value in PCAx1:
            val_str=str(value)[1:-1]
            PCAx1_list.append(float(val_str))
        for value in PCAx2:
            val_str=str(value)[1:-1]
            PCAx2_list.append(float(val_str))
        for value in PCAx3:
            val_str=str(value)[1:-1]
            PCAx3_list.append(float(val_str))
        return(PCAx1_list, PCAx2_list, PCAx3_list)


    def TwoDimensionalPCA(self):
        '''Takes an array for fingerprints, performs 2D pca, and adds the 2D coordinates to the specified dataframe

        Parameters:
        ==========
        fp_array: An array of fingerprints to be used for the pca
        data_df: The dataframe to append the pca coordinates to
        '''
        self.pca_coords=pcaPreprocess(self.fp_array, doPCA=True,n_components=self.nBits)
        pca_x1, pca_x2= self.make_lists_from_2D_coords()
        self.chart_df['x1']=pca_x1
        self.chart_df['x2']=pca_x2

    def ThreeDimensionalPCA(self):
        '''Takes an array for fingerprints, performs 2D pca, and adds the 2D coordinates to the specified dataframe

        Parameters:
        ==========
        fp_array: An array of fingerprints to be used for the pca
        data_df: The dataframe to append the pca coordinates to
        '''
        self.pca_coords=pcaPreprocess(self.fp_array, doPCA=True,n_components=self.nBits)
        pca_x1, pca_x2, pca_x3= self.make_lists_from_3D_coords(self.pca_coords)
        self.chart_df['x1']=pca_x1
        self.chart_df['x2']=pca_x2
        self.chart_df['x3']=pca_x3

    def PlotPCA(self):
        '''Construct plots using the PCA coords from the data_df faceted according to DEL'''
        self.PCA_charts = alt.Chart(self.chart_df).mark_circle().encode(
            x='x1', y='x2',
            color=alt.Color("DEL:N",
                            legend=alt.Legend(title="Type")),
            size=alt.value(10),
            tooltip=["x1", "x2", "DEL:N", 'ID:N']
        ).properties(title=str(self.descriptor) + " PCA Sets", width=300, height=300).facet(
            column='DEL:N'
        ).interactive()

    def PlotOverlayPCA(self):
        '''Construct plots using the PCA coords from the data_df faceted according to DEL'''
        self.PCA_Overlay_chart = alt.Chart(self.chart_df).mark_circle().encode(
            x='x1', y='x2',
            color=alt.Color("DEL:N",
                            legend=alt.Legend(title="Type")),
            size=alt.value(10),
            tooltip=["x1", "x2", "DEL:N", 'ID:N']
        ).properties(title=str(self.descriptor) + " PCA Sets", width=300, height=300).interactive()
    
    def csv_output(self):
        '''Output the PCA results to a csv file'''
        self.chart_df.to_csv('Notebook Outputs/'+str(self.descriptor)+'_PCA_Results.csv')