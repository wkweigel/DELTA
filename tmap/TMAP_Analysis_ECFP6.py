import pandas as pd
import tmap
from faerun import Faerun
from rdkit.Chem import AllChem
from timeit import default_timer as timer
from matplotlib.colors import ListedColormap
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

#Name for output html and js files
name='Combined_DELs_ECFP6'

# The number of permutations used by the MinHashing algorithm
perm = 1024

#Import csv as dataframe
df = pd.read_csv("Combined_DELs.csv") 
df.shape

#Initialize lists
fingerprint=[]
labels=[]
start = timer()

# Calculate ECFP6 decriptors for each compound and add them to the fingerprint list as a tmap VectorUnit 
for i, row in df.iterrows(): 
    ROMol=Chem.MolFromSmiles(row['Smiles'])
    ECFP6 = AllChem.GetMorganFingerprintAsBitVect(ROMol,radius=3, nBits=1024)
    fp_array = np.zeros((0, ), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(ECFP6, fp_array)
    fingerprint.append(tmap.VectorUint(fp_array))
    labels.append(row['Smiles']+'__'+row['ID'])

    if i != 0 and i % 1000 == 0:# Progress Updater. Prints the percentage of the dataset processed every 1000 rows.
            percent=round((100 * i / len(df)),2)
            print('')
            print(str(percent)+'% Complete')
            print(f"Generating the data took {round((timer() - start), 1)}s.")

#Initialize the LSH Forest
lf1 = tmap.LSHForest(perm)

# Add the Fingerprints to the LSH Forest and index
lf1.batch_add(fingerprint)
lf1.index()

#Set tmap configuration parameters with cfg class
cfg = tmap.LayoutConfiguration()
cfg.sl_scaling_type=tmap.ScalingType.RelativeToAvgLength

# Get the coordinates
x, y, s, t, _ = tmap.layout_from_lsh_forest(lf1,cfg)

# Create type catagories using the "DEL" column 
type_labels, type_data = Faerun.create_categories(df["DEL"])

# Plot the data
faerun = Faerun(view="front", clear_color='#FFFFFF', coords=True)
faerun.add_scatter(
    name,
    {   "x": x, 
        "y": y, 
        "c": [type_data], 
        "labels": labels},
    legend_labels=[type_labels],
    point_scale=10,
    max_point_size=10,
    colormap = 'tab10',
    has_legend=True,
    legend_title = ['DELs'],
    series_title = ['Type'],
    categorical= [True],
    shader = 'smoothCircle'
)

#Create an output csv of the tmap data
data_dict={"x": x, "y": y}
data_df=pd.DataFrame.from_dict(data_dict)
data_df.to_csv('tmap_data.csv')

faerun.add_tree(name+"_tree", {"from": s, "to": t}, point_helper=name)

# Choose the "smiles" template to display structure on hover
faerun.plot(name, template="smiles", notebook_height=750)
