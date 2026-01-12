# Analysis pipeline

This document provides a detailed description of the analysis pipeline used in this study. Each section corresponds to a dedicated module in the repository and represents a logical stage of the analysis.

---

## 1. Genes and Data Retrieval

This module retrieves and preprocesses gene- and disease-related data.

### 1.1 DiseasesTree
Builds a hierarchical disease tree based on UK Biobank disease data.

### 1.2 HumanGenesFrameComplete
Uses GenEraGenes and Ensembl to generate data frames for:
- Variant-to-gene mapping
- Conversion of Ensembl identifiers to HGNC symbols

### 1.3 GenAgeGenes_Retrieval
Retrieves ageing-related genes from GenAge for both humans and model organisms.  
Orthology mapping for model organisms is performed using OMA.

### 1.4 DiseasesGenes_Retrieval
Processes GWAS summary statistics from Donertas et al., mapping variants to genes using a ±10 kb distance-based approach.  
This enables gene assignments to:
- Individual diseases (ARDs)
- Disease clusters (ARCs)

### 1.5 Coexpression_Retrieval
Retrieves and processes gene co-expression data from GeneFriends, adapting it for ageing- and disease-related gene sets.

### 1.6 KeggPathways_Retrieval
Retrieves and integrates signalling pathways from KEGG.

---

## 2. Pleiotropy and Disease Association Analysis

This module analyses gene pleiotropy and disease associations.

### 2.1 Overlap_Age_Diseases
Generates UpSet plots showing overlap between disease-associated genes and ageing-related genes.

### 2.2 Pleiotropy_Permutation
Performs permutation tests to assess pleiotropy by comparing GenAge genes against the remaining genome.

### 2.3 DiseaseAssociation_Heatmap
Displays heatmaps illustrating how genes associated with at least one disease overlap across multiple diseases.

### 2.4 Overlap_Age_Pleiotropy
Shows how ageing-related genes overlap with disease genes as a function of ARC pleiotropy, defined as the number of ARCs associated with each gene.

---

## 3. Network Analysis

This module constructs biological networks, computes network-based metrics, and generates machine-learning–ready datasets.

### 3.1 NetworksCreation
Constructs four networks:
- PPI (BioGRID)
- COX90 and COX95 (GeneFriends)
- KEGG (KEGG pathways)

### 3.2 Gene2Gene_Proximities_Closest
Computes shortest-path distance matrices between all gene pairs for each network.

### 3.3 Gene2Diseases_Proximities
Computes shortest-path distances from genes to each ARD and ARC.

### 3.4 MLdatasets_ProximityAndNeighbours
Generates datasets based on:
- Minimum distance
- Average distance
- Number of disease-associated neighbours

### 3.5 MLdatasets_RWR
Computes Random Walk with Restart (RWR) between genes and ARDs/ARCs and generates RWR-based datasets.

### 3.6 Indirect_Interactions_Frame
Computes ARC-Interactors, defined as the number of ARCs connected to a gene through at least one first-order neighbour associated with those pathologies.

### 3.7 NetworkCentralities
Computes network centrality metrics:
- Degree
- Betweenness
- Closeness
- Clustering coefficient  
across networks and gene groups (GenAge, disease genes, high-pleiotropy genes).

### 3.8 PleiotropyAndInteractors_Frame
Generates data frames relating ARC pleiotropy to ARC interactions for each gene across networks.

### 3.9 PleiotropyAndInteractors_Plot
Creates scatter plots of ARC pleiotropy vs. ARC interactions.

### 3.10 DiseaseDistances_ComputationAndHeatmap
Generates heatmaps of gene-to-ARC shortest-path distances, grouped by:
- GenAge Human
- GenAge Model
- Disease genes
- Neighbours of disease genes

### 3.11 Permutations_NetworkAssociations
Performs permutation tests for gene–disease association metrics:
- Minimum distance
- Average distance
- ARC-Interactors
- RWR scores

### 3.12 DiseaseInteractors_BoxPlot
Generates boxplots comparing ARC-Interactors across gene groups.

### 3.13 Permutations_Plots_Styling
Improves the aesthetic styling of permutation plots.

### 3.14 KEGG_Permutation
Computes shortest-path distances from genes to roots and leaves in KEGG networks and performs permutation analyses for all major gene groups.

---

## 4. Co-expression and Tissue Specificity

All analyses in this section are based on GeneFriends data.

### 4.1 Coexpression_Analysis
Analyses intra-set and inter-set co-expression for:
- ARC gene groups
- High- and low-pleiotropy genes
- GenAge Human and Model genes

### 4.2 Coexpression_Plots
Visualizes co-expression results using boxplots and heatmaps.

### 4.3 Specificity_Analysis
Computes tissue specificity (Tau values) for all gene groups.

### 4.4 Specificity_Plots
Visualizes tissue specificity results using boxplots and heatmaps.

### 4.5 CoexpressionAndSpecificity_Plot
Creates a co-expression vs. tissue specificity plot, highlighting the opposing behaviour of ageing-related genes and high-pleiotropy genes.

### 4.6 Coexpression_Permutation
Permutation tests assessing whether extreme co-expression values differ from random expectations.

### 4.7 SpecificityInNetworks
Boxplots of tissue specificity as a function of ARC-Interactions across networks.

### 4.8 Specificity_Permutation
Permutation tests assessing whether extreme tissue specificity values differ from chance.

### 4.9 Tissue_specific_gene_expression
Visualizes gene expression across GTEx tissues for all ARC and GenAge gene sets.

---

## 5. Machine Learning

This module focuses on predicting novel ageing-related genes.

### ML_Pipeline
Contains all core functions required for machine-learning execution, including:
- Data preprocessing
- Feature construction
- Nested cross-validation
- Model evaluation
- Result aggregation

### 5.1 MachineLearning_Runs
Performs predictions using Random Forest classifiers with nested cross-validation, based on features derived from:
- PPI
- KEGG
- COX90
- COX95  

and algorithms based on:
- Minimum distance
- Average distance
- Number of disease neighbours
- RWR

Predictions are performed for:
- Individual network–algorithm combinations
- All networks per algorithm
- All algorithms per network
- Multiplex network integration

### 5.2 Top_Ageing_False_Positives
Selects the best-performing method (Multiplex) and ranks all genes not originally annotated in GenAge Human by their predicted ageing association probability.

---

## 6. Enrichment Analysis

### 6.1 Enrichment
Performs GO Biological Process enrichment for:
- Human-exclusive ageing genes
- Model-organism–exclusive ageing genes
- Combined ageing genes
- High-pleiotropy genes
- Top 30 machine-learning–predicted ageing genes  

Results are compiled into structured tables for interpretation.

---

## Data sources

- UK Biobank
- GenAge
- GeneFriends
- BioGRID
- KEGG
- GTEx
- GWAS summary statistics (Donertas et al.)

---

## Reproducibility notes

- Large datasets and intermediate files are not included.
- Paths and parameters are configurable within each script.
- Random seeds are set where applicable for reproducibility.

## Data availability

All processed and generated datasets produced in this study are available through **Synapse** under project ID **synXXXXXXX**.

Due to size constraints, raw and intermediate data are not hosted directly in this repository.  
Scripts provided here allow full regeneration of all results from the retrieved data sources.
