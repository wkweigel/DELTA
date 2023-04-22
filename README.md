# DELTA

DELTA (DNA Encoded Library Topological Assignment) is a system for classificaion of DNA encoded libraries (DELs) according to the topological arrangement of thier building blocks (BBs). <br><br>
The following repository is included as supporting information for **LINK TO ARTICLE** <br><br>
The included notebooks contain the Python scripts for the following:
  - Creation of DELs with different topologies
  - Fingerprinting and Principal Component Analysis (PCA)
  - Generative Topographic Mapping (GTM)
  - Pairwise distance calculation (Tanimoto, Euclidean, etc.)
  - Druglike property calculation
  - Principal Moment of Ineria (PMI) calculation
  - Tanimoto similarity map analysis 
  
The functions of each are discussed in more detail below.

## Library Generation
![Library generation flowchart](/assets/LibraryGeneration.png)
Running all the cells in this notebook will generate 1,000 random compounds for seven topologically variant DELs. The compounds for each DEL are output to csv files in the "Notebook Outputs" folder.

Libraries are generated from two input csv files:
  - BB file (contains smiles and IDs for all the building blocks)
  - DE file (contains the pool of BBs IDs to use for each Diversity Element)

The BB file is located in the root directory. Each DEL uses a seperate DE file located in the "Library BB Groups" folder. Unique sets of BBs are first generated using the Library Set module. Each set of BBs is then assembled into the final DEL compound using the Reactions module. The assembly of the DEL is perormed with a series of reactions using RDKit and various strings of reaction smarts. See the GenerateLibraries notebook for more information.



  

## Fingerprinting and Principal Component Analysis (PCA)
![Principle Component Analysis](/assets/PCA.png)
The PCA notebook contains scripts for calculating molecular fingerprints for the DEL compounds which are then used to perform 2D or 3D PCA. The available fingerprinting methods include ECFP4, ECFP6, AP, MACCS, MQN, MHFP6<sup>[1](https://github.com/reymond-group/mhfp)</sup>, MXFP<sup>[2](https://github.com/reymond-group/mxfp_python)</sup>, and MAP4<sup>[3](https://github.com/reymond-group/map4)</sup>. 2D plots are generated using Vega-Altair to visulaize the results.

## Generative Topographic Mapping (GTM)
![Generative Topographic Mapping](/assets/GTM.png)<br>
GTM calculations are performed using the ugtm<sup>[4](https://github.com/hagax8/ugtm)</sup> library. The core GTM algoritm was modified slightly by adding a few lines of code to enable plotting using the method of "Responsibility Patterns" used by Klimenko et al.<sup>[5](https://doi.org/10.1021/acs.jcim.6b00192)</sup> These modified files are inluded in the "ugtm" folder and can be used to replace the main files from the ugtm package. A csv file containing a collection of 10,000 ChEMBL compounds is included in the root directory which is used as the training set for GTM calculations. With the included script, GTM is performed using a 1024-bit ECFP6 descriptor by defualt and calculations are typically complete within 20-30 minutes on a desktop PC. The 2D projections of the DELs are visualized using 2D plots and also may be displayed using binned heatmaps. The coordinate data from the GTM calculations can be output as a csv file in the Notebook Outputs folder for use in furthur analysis.

## Pairwise Distance Analysis
![Generative Topographic Mapping](/assets/Tanimoto.png)<br>
Pairwise distances can be calculated from compounds or coordinate data from an input csv file. Pairwise Tanimoto similarities use smiles and can be calculated using any of the outputs from the GenerateLibraries notebook. Cartesian-based distance metrics require coordinate data  
