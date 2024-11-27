This repository contains two folders:

Code Folder: Contains the code used in this project.
Data Folder: Divided into two parts:
Obtained Data
Generated Data
The code references the locations within the data folders to save and load information. However, due to GitHub's file size limitations, not all files could be uploaded. As a result, some datasets or databases were excluded from Kick Hof because they exceed this limit. These include both obtained and generated datasets:

Missing Retrieved Data (data retrievd directly from publicibly available databases):

- Biogrid's PPI:  Downloadable at BioGrids's Human Physical Interactions, "BIOGRID-MV-Physical-4.4.204.tab3.txt". The symbolical Github path for these files is "maglab/Gustavo/Data/Retrieved/Networks/BIOGRID-MV-Physical-4.4.204.tab3.txt"
  
- GeneFriends Coexpression: Downloadable at GeneFriends DropBox "GENEFRIENDS/homo_sapiens_sapiens/SRA/human_genes_correlation_matrix.tsv". The symbolical Github path for these files is "maglab/Gustavo/Data/Retrieved/Networks/human_genes_correlation_matrix.tsv"
  
- UKBB GWAS: Downloadable at "BioStudies (S- BSST407)". The symbolical Github's folder path for these files is "maglab/Gustavo/Data/Retrieved/Ukbb".


Missing Generated Data (Primarily consisting of gene distance and proximity matrices, as well as relationships between diseases and coexpression):
Files include:

- Disease self- and cross-coexpression relationships
  - maglab/Gustavo/Data/Generated/Coexpression/Coexpression_Pairwise_Frame.rds
-- maglab/Gustavo/Data/Generated/Coexpression/Coexpression_Scores_Frame.rds

Gene-Gene distanes for PPI
maglab/Gustavo/Data/Generated/PPI/Proximity_and_Distances/GenGen_DistanceFrame.rds
maglab/Gustavo/Data/Generated/PPI/Proximity_and_Distances/GenGen_ProximityFrame.rds

Gene-Gene distanes for COX90
maglab/Gustavo/Data/Generated/COX90/Proximity_and_Distances/GenGen_DistanceFrame.rds
maglab/Gustavo/Data/Generated/COX90/Proximity_and_Distances/GenGen_ProximityFrame.rds

Gene-Gene distanes for COX95
maglab/Gustavo/Data/Generated/COX95/Proximity_and_Distances/GenGen_DistanceFrame.rds
maglab/Gustavo/Data/Generated/COX95/Proximity_and_Distances/GenGen_ProximityFrame.rds

Gene-Gene distanes for KEGG
maglab/Gustavo/Data/Generated/KEGG/Proximity_and_Distances/GenGen_DistanceFrame.rds
maglab/Gustavo/Data/Generated/KEGG/Proximity_and_Distances/GenGen_ProximityFrame.rds



X1
X2
If you would like access to these files, please feel free to contact us, and we will be happy to share them.
