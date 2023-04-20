import pandas as pd
import matplotlib.pylab as plt
import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import pandas as pd
import numpy as np
import csv

#GLOBALS
global elements
global linkers
elements='ABCDEFGHIJKLMNOP'
linkers= 'abcdefgh'

# stores the vertices in the graph
vertices = []
#stores the number of vertices in the graph
vertices_no = 0
graph = []
#Stores edge list values
edge_list=[]

# Main parsing function for complex DELTA strings. 

# Imports a csv file and converts it to a dict in the form of: {id1:smiles1,...} 
def ID_Dict_from_csv(file_location):
    df=pd.read_csv(file_location) #if copying path from windows change all '\' into '/'
    frame1=df['ID']
    frame2=df['Smiles']
    global ID_Dict
    if products.prod_dict=={}:
        ID_Dict={}
        ID_Dict['Int']={} #adds an empty dictonary item to hold smiles for reaction intermediates
    for idx, i in enumerate(frame1):
        ID_Dict[frame1[idx]]=frame2[idx]


def rxn_from_BB_set(Rxn, bb1, bb2, Draw_Products=False):
    #products.MolProds=[]
    #products.SmilesProds=[]
    #products.ProdLabels=[]
    

    reaction=regular_rxn(bb1, bb2)

    if Rxn == 'R1':
        reaction.R1()
    if Rxn == 'R2':
        reaction.R2()
    if Rxn == 'R3':
        reaction.R3()
    


def find_complex_edges(X, ShowOutputs=False):  

    global E_items
    E_items=[]
    RXN_items=[]
    ALL_items=[]
    global edges
    edges={}
    
    global linkers 
    
    global ID_list
    global top_string


    add_RXN=False
    add_E=False
    temp_E=''
    temp_RXN=''
    for char in X:
        if char == '|' and add_E is False:
                add_E=True
                continue
        if char =='<':
                add_RXN=True
                continue

        # extracts a str between set of '<>' 
        if add_RXN is True:
                if char == '>':
                        add_RXN=False
                        RXN_items.append(temp_RXN)
                        ALL_items.append('<' + temp_RXN + '>')
                        temp_RXN=''
                        continue
                else:
                        temp_RXN=temp_RXN + char

        # extracts a str between set of '| |'  
        elif add_E is True:
                if char == '|':
                        add_E=False
                        E_items.append(temp_E)
                        ALL_items.append(temp_E)
                        temp_E=''
                        continue
                else:
                        temp_E=temp_E + char

    E_items.append('DNA')
    ALL_items.append('DNA')       
    c1_search=True
    for idx, item in enumerate(ALL_items):
        #Search and ID the first cycle element
        if '!' in item and c1_search is True:
            c1=item
            c1_search=False
        if '<' in item and '(' in item:
            b3=ALL_items[idx-1]
            b4=ALL_items[idx+1]
            b5=ALL_items[idx+2]
            get_rxn=False
            temp_rxn=''
            for char in item:
                if char == '(':
                    get_rxn=True
                    continue
                if char == ')':
                    r2=temp_rxn
                    get_rxn=False
                    temp_rxn=''
                    continue
                if get_rxn is True:
                    temp_rxn=temp_rxn+char
                    continue

                if char == '>':
                    r3=temp_rxn
                    break
                elif char != '<':
                        temp_rxn=temp_rxn+char
            edges[b3,b4]=r2
            edges[b3,b5]=r3
            continue

        if '<' in item and '!' in item:
            b3=ALL_items[idx-1]
            b4=ALL_items[idx+1]
            
            get_rxn=False
            cyc_rxn=False
            temp_rxn=''
            for char in item:
                #First '!'
                if char == '!' and cyc_rxn is False:
                    cyc_rxn=True
                    continue
                #Second '!'
                if char == '!' and cyc_rxn is True:
                    get_rxn=True
                    continue
                #Adds characters between the set of "!" to temp_rxn 
                if get_rxn is False and cyc_rxn is True:
                    temp_rxn=temp_rxn+char
                #Once full rxn ID is obtained, name final ID and reset 
                if get_rxn is True and cyc_rxn is True:  
                    c=temp_rxn
                    get_rxn=False
                    cyc_rxn=False
        
                    continue
                if char == '>':
                    r4=temp_rxn
                    break
                elif char != '<' and cyc_rxn is False:
                    temp_rxn=temp_rxn+char
            edges[c1,b3]=c
            edges[b3,b4]=r4
            
            continue

        if '<' in item:
            b1=ALL_items[idx-1]
            b2=ALL_items[idx+1]
            temp_rxn=''
            for char in item:
                if char == '<':
                    continue
                if char == '>':
                    r1=temp_rxn
                    break
                else:
                    temp_rxn=temp_rxn+char
            edges[b1,b2]=r1
    if ShowOutputs is True:
        print('')
        print('String Input:')
        print(X)
        print('')

        print('Elements:')
        print(E_items)
        print('')

        print('Reactions:')
        print(RXN_items)
        print('')

        print('Edge Dictionary:')
        print(edges)
        print('')

        get_bb_dict(E_items)
        complex_processing(E_items)

       # print('Final ndarray:')
        #array= np.array(graph)
        #print(array)

#A special variation of the the regular "add_vertex" function that acesses the bb_dict
def add_complex_vertex(v):
    global graph
    global vertices_no
    global vertices
    if v in vertices:
        pass
    else:
        vertices.append(v)
        vertices_no = len(vertices)
    #updates existing vertexes with appropriate number of elements (ie. expands extisting vertexes as new elements are added to list)
    if vertices_no > 1:
            for vertex in graph:
                    vertex.append(0)
    
    #creates new vertex with appropriate number of elements and adds it to ndarray
    temp = []
    for i in range(vertices_no):
        if i == vertices_no-1:
            for char in v:
                if v == 'DNA':
                    temp.append(0)
                    break
                if char in elements:
                    temp.append(bb_dict[char])
                    break       
        elif vertices_no > 1:
            temp.append(0)
    graph.append(temp)

def complex_processing(X):
         # stores the vertices in the graph matrix
        for vert in X:
                add_complex_vertex(vert)
        # Adds edge according to the 
        for (k1,k2), v1 in edges.items():
                add_edge(k1,k2,v1)

# Add an edge between vertex v1 and v2 with edge weight e
def add_edge(v1, v2, e):
    global graph
    global vertices_no
    global vertices
    if v1 not in vertices:
        pass
        #print("Vertex ", v1, " does not exist.")
    elif v2 not in vertices:
        pass
        #print("Vertex ", v2, " does not exist.")
    else:
        index1 = vertices.index(v1)
        index2 = vertices.index(v2)
        graph[index1][index2] = e

#Removes syntatical charaters from the edge list               
def get_Rxn_Dict(dict):
    global Rxn_Dict
    Rxn_Dict={}
    for (k1,k2), v in dict.items():
        temp_k1=''
        temp_k2=''
        for char in k1:
            if char == '!':
                continue
            if char == '(':
                continue
            if char == ')':
                continue
            else:
                temp_k1=temp_k1+char
        for char in k2:
            if char == '!':
                continue
            if char == '(':
                continue
            if char == ')':
                continue
            else:
                temp_k2=temp_k2+char
        Rxn_Dict[temp_k1, temp_k2]=v

def get_bb_dict(X,print=False):
    global elements
    elements= 'ABCDEFGHIJKLMNOP'
    global bb_dict
    bb_dict={}
    for element in X:
            temp_e=''
            temp_bb=''
            if element == 'DNA':
                    continue
            for char in element:
                    if char in elements:
                            temp_e=char
                    if char.isdigit():
                            temp_bb=temp_bb+char
            bb_dict[temp_e]=temp_bb
    if print is True:
        print('Building Block Dictionary:')
        print(bb_dict)
        print('')  
        

# Creates class for holding product results in Mol and Smiles formats 
class rxn_products:
    def __init__(self):
        self.MolProds=[]
        self.SmilesProds=[]
        self.ProdLabels=[]
        self.prod_dict={}
    def add_mol_prod(self, prod):
        self.MolProds.append(prod)
    def add_smiles_prod(self, prod):
        self.SmilesProds.append(prod)
    def assemble_prod_dict(self):
        self.prod_dict={}
        for i in range(len(self.SmilesProds)):
            i_smiles=self.SmilesProds[i]
            i_label=self.ProdLabels[i]
            self.prod_dict[i_label]=i_smiles
        return(self.prod_dict)

class temp_products:
    def __init__(self):
        self.SmilesProd=''
        self.ProdLabel=''


class starting_materials:
    def __init__(self):
        self.Amines={}
        self.Carb_Acids={}

    def add_amines(self, BB):
        self.Amines=ID_Dict[BB]

    def add_carb_acids(self, BB):
        self.Carb_Acids=ID_Dict[BB]

    def add_starting_materials(self, rxn_id, BB1, BB2):
        if rxn_id == "R1":
            if str(type_dict[BB1])[-2:]== 'NH':
                self.add_amines(bb1)
            if str(type_dict[BB2])[-2:]== 'NH':
                self.add_amines(bb2)
            if str(type_dict[BB1])[-2:]== 'CA':
                self.add_carb_acids(bb1)
            if str(type_dict[BB2])[-2:]== 'CA':
                self.add_carb_acids(bb2)
        if rxn_id == "R2":
            if str(type_dict[BB1])[-2:]== 'NH':
                self.add_amines(bb1)
            if str(type_dict[BB2])[-2:]== 'NH':
                self.add_amines(bb2)
            if str(type_dict[BB1])[-2:]== 'CA':
                self.add_carb_acids(bb1)
            if str(type_dict[BB2])[-2:]== 'CA':
                self.add_carb_acids(bb2)
        if rxn_id == "R3":
            if str(type_dict[BB1])[-2:]== 'NH':
                self.add_amines(bb1)
            if str(type_dict[BB2])[-2:]== 'NH':
                self.add_amines(bb2)
            if str(type_dict[BB1])[-2:]== 'CA':
                self.add_carb_acids(bb1)
            if str(type_dict[BB2])[-2:]== 'CA':
                self.add_carb_acids(bb2)



class regular_rxn:
    def __init__(self, BB1, BB2):
        self.Mol1=ID_Dict[BB1]
        self.Mol2=ID_Dict[BB2]
        self.BB1=BB1
        self.BB2=BB2
    # <R1> Amide Bond Formation (Carboxylic Acid + Amine) 
    def R1(self):
        reaction = "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
        amide_bond_formation = AllChem.ReactionFromSmarts(reaction)

        #Defines reactants using the appropriate starting_materials class attribute
        reactant1=sm.Carb_Acids 
        reactant2=sm.Amines

        reactant1_mol = Chem.MolFromSmiles(reactant1)
        reactant2_mol = Chem.MolFromSmiles(reactant2)
        amidep = Chem.MolFromSmarts('[N;$(NC=[O,S])]')
        for match1 in reactant1_mol.GetSubstructMatches(amidep):
            reactant1_mol.GetAtomWithIdx(match1[0]).SetProp('_protected','1')
        for match2 in reactant2_mol.GetSubstructMatches(amidep):
            reactant2_mol.GetAtomWithIdx(match2[0]).SetProp('_protected','1')
        product = amide_bond_formation.RunReactants ([reactant1_mol, reactant2_mol])
        smiles=Chem.MolToSmiles(product[0][0])
        label=str(self.BB1)+'_'+str(self.BB2)
        temp_prod.SmilesProd=smiles
        temp_prod.ProdLabel=label
        
    # R2 Amide Coupling > Boc Detrotection
    def R2(self):
        global temp_prod
        reaction = "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
        amide_bond_formation = AllChem.ReactionFromSmarts(reaction)
      
        #Defines reactants using the appropriate starting_materials class attribute
        reactant1=sm.Carb_Acids 
        reactant2=sm.Amines
        

        reactant1_mol = Chem.MolFromSmiles(reactant1)
        reactant2_mol = Chem.MolFromSmiles(reactant2)
        amidep = Chem.MolFromSmarts('[N;$(NC=[O,S])]')
        for match1 in reactant1_mol.GetSubstructMatches(amidep):
            reactant1_mol.GetAtomWithIdx(match1[0]).SetProp('_protected','1')
        for match2 in reactant2_mol.GetSubstructMatches(amidep):
            reactant2_mol.GetAtomWithIdx(match2[0]).SetProp('_protected','1')
        product = amide_bond_formation.RunReactants ([reactant1_mol, reactant2_mol])
        temp_smiles=Chem.MolToSmiles(product[0][0])
        label=str(self.BB1)+'_'+str(self.BB2)

        reaction2='[#7:1]-C(=O)-O-[$(C(-[CH3])(-[CH3])(-[CH3]))]>>[*:1]'
        boc_deprotection=AllChem.ReactionFromSmarts(reaction2)
        reactant1_mol = Chem.MolFromSmiles(temp_smiles)
        product = boc_deprotection.RunReactants ([reactant1_mol])
        smiles=Chem.MolToSmiles(product[0][0])
        temp_prod.SmilesProd=smiles
        temp_prod.ProdLabel=label
        



    # R3 Amide Coupling > Fmoc Detrotection
    def R3(self):
        global temp_prod
        reaction = "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]" #Defines reaction using SMARTS
        amide_bond_formation = AllChem.ReactionFromSmarts(reaction) #Creates the rdkit reaction object from the SMARTS

        reactant1=sm.Carb_Acids 
        reactant2=sm.Amines
        
        reactant1_mol = Chem.MolFromSmiles(reactant1)
        reactant2_mol = Chem.MolFromSmiles(reactant2)
        amidep = Chem.MolFromSmarts('[N;$(NC=[O,S])]')
        for match1 in reactant1_mol.GetSubstructMatches(amidep):
            reactant1_mol.GetAtomWithIdx(match1[0]).SetProp('_protected','1')
        for match2 in reactant2_mol.GetSubstructMatches(amidep):
            reactant2_mol.GetAtomWithIdx(match2[0]).SetProp('_protected','1')
        product = amide_bond_formation.RunReactants ([reactant1_mol, reactant2_mol])
        temp_smiles=Chem.MolToSmiles(product[0][0])
        label=str(self.BB1)+'_'+str(self.BB2)
    
        reaction2='[#7:1]-C(=O)-O-C-[$([C]1[cR2][cR2][cR2][cR2]1)]>>[*:1]'
        fmoc_deprotection=AllChem.ReactionFromSmarts(reaction2)
        reactant1_mol = Chem.MolFromSmiles(temp_smiles)
        product = fmoc_deprotection.RunReactants ([reactant1_mol])
        smiles=Chem.MolToSmiles(product[0][0])
        temp_prod.SmilesProd=smiles
        temp_prod.ProdLabel=label


# R4 Intramolcular Amide Bond Formation (list input)
def R4(input_list):
    reaction = "([C:1](=[O:2])-[OD1].[N!H0:3])>>[C:1](=[O:2])[N:3]"
    cyclic_amide_bond_formation = AllChem.ReactionFromSmarts(reaction)
    cyclic_products=[]

    reactant1=input_list

    for smiles in reactant1:
        reactant1_mol = Chem.MolFromSmiles(smiles)
        amidep = Chem.MolFromSmarts('[N;$(NC=[O,S])]')
        for match1 in reactant1_mol.GetSubstructMatches(amidep):
            reactant1_mol.GetAtomWithIdx(match1[0]).SetProp('_protected','1')
        product = cyclic_amide_bond_formation.RunReactants ([reactant1_mol])
        smiles=Chem.MolToSmiles(product[0][0])
        for i in range(len(product)):
            cyclic_products.append(smiles)
    return(cyclic_products)


def final_deprotection(input_df, PG):
    deprotected_smiles=[]
    if PG=='fmoc':
        reaction='[#7:1]-C(=O)-O-C-[$([C]1[cR2][cR2][cR2][cR2]1)]>>[*:1]'
        deprotection=AllChem.ReactionFromSmarts(reaction)
    if PG=='boc':
        reaction='[#7:1]-C(=O)-O-[$(C(-[CH3])(-[CH3])(-[CH3]))]>>[*:1]'
        deprotection=AllChem.ReactionFromSmarts(reaction)
    if PG=='ester':
        reaction='[CX3:1](=[OX1:2])-[OX2:3]-[$([CH3]),$([CH2]-[CH3])]>>[C:1](=[O:2])-[O:3]'
        deprotection=AllChem.ReactionFromSmarts(reaction)
    
    reactant1=input_df.values.tolist()

    for smiles in reactant1:
            reactant1_mol = Chem.MolFromSmiles(smiles)
            product = deprotection.RunReactants ([reactant1_mol])
            smiles=Chem.MolToSmiles(product[0][0])
            for i in range(len(product)):
                deprotected_smiles.append(smiles)
    return(deprotected_smiles)

def reaction_counter(string):
    rxn_counter=0
    for character in string:
        if character=='R':
            rxn_counter+=1
    return(rxn_counter)


def run_reaction():
    #Get the number of reactions used for the synthesis 
    rxn_num=complex_string.count('R')
    ID_dict={}

    #load the specified building block sets
    try:
        df_BBsets=pd.read_csv(BBset_path)
    except:
         df_BBsets=Random_df

    #parse the reaction input string
    find_complex_edges(complex_string,ShowOutputs=True)

    #create a reaction dict containing the starting materials and reactions
    get_Rxn_Dict(edges)

    #load the list of bulding block structures
    ID_Dict_from_csv(BB_IDs_path)
    global prod_dict
    prod_dict={}

    #iterate over all the building block sets and asseble them according to the rxn_dict
    for i, row in df_BBsets.iterrows():
        global bb1
        global bb2
        global sm

        r_num=0

        #iterate over all the reactions in the rxn_dict
        for [k1,k2],v in Rxn_Dict.items(): 
            r_num=r_num+1 #Reaction counter. Flow control for handling the writing of product labels and intermediates.
            Int='Z' #Specifies the placeholder BB used for the growing reaction intermediate 

            #Class object creation for starting materials. Note this class object is refreshed in each loop
            sm=starting_materials()
            
            #Assign the building blocks as the starting materials and run the reaction
            if r_num == 1: 
                Int_Dict={}
                bb1=row[k1]
                bb2=row[k2]
                sm.add_starting_materials(v,k1,k2)
                rxn_from_BB_set(v, bb1, bb2)
                ID_Dict[temp_prod.ProdLabel]=temp_prod.SmilesProd
            else:
                bb1=temp_prod.ProdLabel
                bb2=row[k2]
                sm.add_starting_materials(v,Int,k2)
                rxn_from_BB_set(v, bb1, bb2)
            if r_num==rxn_num:
                prod_dict[temp_prod.ProdLabel]=temp_prod.SmilesProd
                break
            else:
                ID_Dict[temp_prod.ProdLabel]=temp_prod.SmilesProd


   

    #create the output dataframe using the product dictionary 
    prod_df=pd.DataFrame()
    DELlist=[]
    for i in range(len(prod_dict.keys())):
        DELlist.append(DELname)
    prod_df['DEL']=DELlist
    prod_df['ID']=prod_dict.keys()
    prod_df['Smiles']=prod_dict.values()
    return(prod_df)


