# DELTA

DELTA (DNA Encoded Library Topological Assignment) is a system for classification of DNA encoded libraries (DELs) according to the topological arrangement of their building blocks (BBs). <br><br>
The following repository is included as supporting information for the associated publication: <br>
[Evaluation of the Topology Space of DNA-Encoded Libraries](https://doi.org/10.1021/acs.jcim.3c01008)  <br><br>
The included notebooks contain the Python scripts for the following:
  - Creation of DELs with different topologies
  - Fingerprinting and Principal Component Analysis (PCA)
  - Generative Topographic Mapping (GTM)
  - Pairwise distance calculation (Tanimoto, Euclidean, etc.)
  - Druglike property calculation
  - Principal Moment of Inertia (PMI) calculation
  - Tanimoto similarity map analysis 
  - TMAP generation 
  
The functions of each are discussed in more detail below.

### Usage and set-up
RDKit must be installed within its own virtual environment. If you don't already have RDKit installed, read the included setup instructions.
The python modules necessary to support all the notebooks can be installed within the RDKit venv with pip using the included requirements.txt file.

The contents of this this repository are intended to serve as an example workflow using only small sets of compounds. 
Everything needed for the creation of seven topologically variant libraries is included. The libraries must be created first using the GenerateLibraries notebook. The rest of the notebooks are designed to work with these library output files. 


## Library Generation
![Library generation flowchart](/assets/LibraryGeneration.png)
Running all the cells in this notebook will generate 1,000 random compounds for seven topologically variant DELs. The compounds for each DEL are output to csv files in the "Notebook Outputs" folder.

Libraries are generated from two input csv files:
  - BB file (contains smiles and IDs for all the building blocks)
  - DE file (contains the pool of BBs IDs to use for each Diversity Element)

The BB file is located in the root directory. Each DEL uses a separate DE file located in the "Library BB Groups" folder. Unique sets of BBs are first generated using the Library Set module. Each set of BBs is then assembled into the final DEL compound using the Reac-tions module. The assembly of the DEL is performed with a series of reactions using RDKit and various strings of reaction smarts. See the GenerateLibraries notebook for more information.



  

## Fingerprinting and Principal Component Analysis (PCA)
![Principle Component Analysis](/assets/PCA.png)
The PCA notebook contains scripts for calculating molecular fingerprints for the DEL compounds which are then used to perform 2D or 3D PCA. The available fingerprinting methods include ECFP4, ECFP6, AP, MACCS, MQN, MHFP6<sup>[1](https://github.com/reymond-group/mhfp)</sup>, MXFP<sup>[2](https://github.com/reymond-group/mxfp_python)</sup>, and MAP4<sup>[3](https://github.com/reymond-group/map4)</sup>. 2D plots are generated using Vega-Altair to visualize the results.

## Generative Topographic Mapping (GTM)
![Generative Topographic Mapping](/assets/GTM.png)<br>
GTM calculations are performed using the ugtm<sup>[4](https://github.com/hagax8/ugtm)</sup> library. The core GTM algorithm was modified slightly by adding a few lines of code to enable plotting using the method of "Responsibility Patterns" used by Klimenko et al.<sup>[5](https://doi.org/10.1021/acs.jcim.6b00192)</sup> These modified files are included in the "ugtm" folder and can be used to replace the main files from the ugtm package. A csv file containing a collection of 10,000 ChEMBL compounds is included in the root directory which is used as the training set for GTM calculations. With the included script, GTM is performed using a 1024-bit ECFP6 descriptor by default and calculations are typically complete within 20-30 minutes on a desktop PC. The 2D projections of the DELs are visualized using 2D plots and also may be displayed using binned heatmaps. The coordinate data from the GTM calculations can be output as a csv file in the Notebook Outputs folder for use in further analysis.

## Pairwise Distance Analysis
![Tanimoto Distributions](/assets/Tanimoto.png)<br>
Pairwise distances can be calculated from compounds or coordinate data from an input csv file. Pairwise Tanimoto similarities use smiles and can be calculated using any of the outputs from the GenerateLibraries notebook. Cartesian-based distance metrics require coordinate data that can be produced using the outputs from the PCA or GTM notebooks. Distance distributions grouped by library can be visualized using violin plots that are generated using Seaborn.

## Druglike Properties Analysis
![MW Distributions](/assets/MW.png)<br>
Properties contributing to druglikeness (Lipinski parameters) can be calculated for the smiles in an input csv file using RDKit. The note-book currenly evaluates druglikeness based on the following criteria: molecular weight < 500 g/mol, cLogP ≤ 5, H-donors ≤ 5, H-acceptors ≤ 10, rotatable bonds ≤ 10, polar surface area < 140 Å2.  An overall report for each Lipinski parameter can be generated:

```
Compounds Analyzed: 7000
mol_wt 3485 49.8%
logp 6849 97.8%
h_donors 6750 96.4%
h_acceptors 6919 98.8%
rotatable_bonds 2491 35.6%
polar_surface_area 4236 60.5%
All 874 12.5%
```

Property distributions grouped by library for each property can be shown as violin plots. The values for the calculated properties may also be output as a csv file for all the compounds or limited to only those that pass all the druglikeness checks.

## Principal Moment of Inertia Plots
![Principal Moment of Inertia](/assets/PMI.png)<br>
The PMI notebook uses the smiles in an input csv file to calculate normalized PMI ratios (NPRs) from the 3D-optimized geometry for each compound. The NPRs are then used to create ternary heatmaps<sup>[6](https://github.com/marcharper/python-ternary)</sup> based on each compound's rod-, disk-, or sphere-likeness. The NPMs and ternary coordinants can also be output to a csv file.

## TMAP
![TMAP](/assets/TMAP.png)<br>
TMAP analysis may be performed using the Python file located in the "tmap" folder. Note that the tmap library<sup>[7](https://tmap.gdb.tools/#simple-graph)</sup> must be run within linux or WSL using Python 3.9. This is best accomplished by creating a seperate RDKit venv for using this older Python version. Faerun<sup>[8](https://github.com/reymond-group/faerun-python)</sup> must also be installed for tmap visualization.
