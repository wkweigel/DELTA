import pandas as pd
import random

class BBset:
    r'''A class for creating and holding sets of Building Blocks used to create the individual members of a DEL

    Parameters:
    ==========
    
    DiversityElements_path: A csv containing the desired BuildingBlock IDs to use for each diversity element.
    '''
    def __init__(self, DiversityElements_path):
        self.BBgroup= pd.read_csv(DiversityElements_path)
        self.BBset=[]
        self.BBlist=[]

        self.BBtypes=len(self.BBgroup.columns)

        if self.BBtypes==1:
            self.BuildingBlock1=[bb for bb in self.BBgroup['BB1'].dropna().unique()]
        if self.BBtypes==2:
            self.BuildingBlock1=[bb for bb in self.BBgroup['BB1'].dropna().unique()]
            self.BuildingBlock2=[bb for bb in self.BBgroup['BB2'].dropna().unique()]
        if self.BBtypes==3:
            self.BuildingBlock1=[bb for bb in self.BBgroup['BB1'].dropna().unique()]
            self.BuildingBlock2=[bb for bb in self.BBgroup['BB2'].dropna().unique()]
            self.BuildingBlock3=[bb for bb in self.BBgroup['BB3'].dropna().unique()]
        if self.BBtypes==4:
            self.BuildingBlock1=[bb for bb in self.BBgroup['BB1'].dropna().unique()]
            self.BuildingBlock2=[bb for bb in self.BBgroup['BB2'].dropna().unique()]
            self.BuildingBlock3=[bb for bb in self.BBgroup['BB3'].dropna().unique()]
            self.BuildingBlock4=[bb for bb in self.BBgroup['BB4'].dropna().unique()]


    def make_combinatorial_set(self):
        r'''Assemble a full combinatorial set of compounds for the given diversity elements

        Returns:
        =======
        Pandas dataframe of all unique combinatorial BB combinations (rows) for x diversity elements (columns).
        '''
        self.BBtypes=len(self.BBgroup.columns)

        if self.BBtypes==2:
            for bb1 in self.BuildingBlock1:
                for bb2 in self.BuildingBlock2:
                        self.BBset=[bb1,bb2]
                        self.BBlist.append(self.BBset)

        if self.BBtypes==3:  
            for bb1 in self.BuildingBlock1:
                for bb2 in self.BuildingBlock2:
                    for bb3 in self.BuildingBlock3:
                        self.BBset=[bb1,bb2,bb3]
                        self.BBlist.append(self.BBset)

        if self.BBtypes==4:
            for bb1 in self.BuildingBlock1:
                for bb2 in self.BuildingBlock2:
                    for bb3 in self.BuildingBlock3:
                        for bb4 in self.BuildingBlock4:
                            self.BBset=[bb1,bb2,bb3,bb4]
                            self.BBlist.append(self.BBset)

        if self.BBtypes==2:
            df_out=pd.DataFrame(self.BBlist, columns=['A','B'])
        if self.BBtypes==3:
            df_out=pd.DataFrame(self.BBlist, columns=['A','B','C'])
        if self.BBtypes==4:
            df_out=pd.DataFrame(self.BBlist, columns=['A','B','C','D'])
        return(df_out)
    

    def get_random_DEs(self):
        '''Generates a random set of Building Blocks to use for each diversity element'''
        if self.BBtypes==1:        
            self.BBset=[(random.choice(self.BuildingBlock1))]
            return(self.BBset)
        if self.BBtypes==2:
            self.BBset=[(random.choice(self.BuildingBlock1)), (random.choice(self.BuildingBlock2))]
            return(self.BBset)
        if self.BBtypes==3:
            self.BBset=[(random.choice(self.BuildingBlock1)), (random.choice(self.BuildingBlock2)), (random.choice(self.BuildingBlock3))]
            return(self.BBset)
        if self.BBtypes==4:
            self.BBset=[(random.choice(self.BuildingBlock1)), (random.choice(self.BuildingBlock2)), (random.choice(self.BuildingBlock3)),(random.choice(self.BuildingBlock4))]
            return(self.BBset)


    def make_random_set(self, n):
        '''Assemble n number of random sets from the given Building Blocks.

        Parameters:
        ===========
        n: The desired number of random sets to generate

        Returns:
        =======
        Pandas dataframe of n random BB combinations (rows) for x diversity elements (columns).
        '''

        while len(self.BBlist) <= n-1:
            self.BBset=self.get_random_DEs()
            if self.BBset not in self.BBlist:
                self.BBlist.append(self.BBset)
        

        if self.BBtypes==1:
            df_out=pd.DataFrame(self.BBlist, columns=['A'])
        if self.BBtypes==2:
            df_out=pd.DataFrame(self.BBlist, columns=['A','B'])
        if self.BBtypes==3:
            df_out=pd.DataFrame(self.BBlist, columns=['A','B','C'])
        if self.BBtypes==4:
            df_out=pd.DataFrame(self.BBlist, columns=['A','B','C','D'])
        return(df_out)





