# Data organization

This document provides a detailed description of the data organization used in this study.  
The repository is structured to clearly distinguish between externally retrieved data sources and datasets generated during the analyses.

The data directory is organized into two main sections:
- **Retrieved** – external datasets and sources used as input for all analyses
- **Generated** – data produced during the analyses presented in this work

---

## 1. Retrieved

Contains all external data sources used in this study.

### 1.1 Genes_and_diseases

Defines gene sets and gene–disease associations.

Subdirectories include:

- **Diseases**
  - Gene–ARC and gene–ARD associations
  - Represented using UK Biobank disease codes and disease names

- **Genes_List**
  - Plain-text gene lists used throughout the study

- **HAGR**
  - GenAge Human and GenAge Model Organism datasets
  - Provided as original CSV files and R objects

- **OMA**
  - Homolog mappings between human genes and model organism genes

- **Ranges**
  - `GenomicRanges` objects mapping GWAS SNPs to genes
  - Based on summary statistics from Donertas et al. (2022)

- **Showcase**
  - Public UK Biobank resources, including:
    - non-cancer self-diagnosed disease hierarchy
    - UK Biobank field definitions
    - disease age-of-onset clustering derived from Donertas et al. (2022)
  - Diseases belonging to clusters 1 (exponential increase) and 2 (linear increase) are defined as ARDs

---

### 1.2 Network_Sources

Data used to construct all networks in this study.

- **Coexpression**
  - GeneFriends coexpression matrices (original downloads)
  - Derived COX90 and COX95 matrices
  - Gene and edge objects using HGNC and Ensembl identifiers

- **KEGG**
  - Individual pathways retrieved using the `kegglincs` library
  - Integrated KEGG pathway objects stored in R format

- **PPI**
  - BIOGRID-MV-Physical-4.4.249.Tab3
  - Physical protein–protein interaction data

---

### 1.3 Specificity

- **Tau_gene_V8** (Palmer et al., 2021)
- Derived from GTEx
- Contains:
  - average gene expression values across tissues
  - tissue specificity (Tau) scores per gene

---

### 1.4 Ukbb

- UK Biobank summary statistics for self-reported diseases
- Source of GWAS SNPs used for gene mapping

---

## 2. Generated

Contains all results generated during the analyses presented in this study.

### 2.1 Genes_Enrichment

- Gene Ontology (GO) enrichment analyses (Excel format)
- Includes enrichment for:
  - human-specific ageing genes
  - model-organism-specific ageing genes
  - human–model overlap
  - high-pleiotropy genes
  - genes predicted as ageing-associated
- Each analysis uses a dedicated genetic background, as described in the main text

---

### 2.2 Networks_and_predictions

Contains all network constructions, topological analyses, and ageing gene predictions.

#### 2.2.1 Networks

Includes four core networks that constitute the pillars of this study, as well as three derived integration strategies.

##### 2.2.1.1 Core networks

- **PPI**
- **KEGG**
- **COX90**
- **COX95**

Each core network contains the following subdirectories:

- **Lists**
  - R objects (`rds`) containing:
    - node lists
    - edge lists
    - complete `igraph` network objects

- **Proximity_and_Distances**
  - Matrices of shortest-path distances and proximity scores between all gene pairs

- **Ageing_Prediction**
  - **Datasets**
    - Machine learning feature matrices
  - **Predictions**
    - Ageing gene prediction results

Machine learning features include:
- Closest Proximity (minimum distance to disease genes)
- Average Proximity (mean distance to disease genes)
- Neighbours (number of disease-associated neighbors)
- RWR Proximity (Random Walk with Restart)

All features are computed separately for:
- age-related disease clusters (ARCs)
- age-related diseases (ARDs)

Both imputed and non-imputed datasets are provided.

---

##### 2.2.1.2 Derived network integrations

Generated directly from the four core networks:

- **PER-ALGORITHM**
  - Combines all networks for a given proximity algorithm

- **PER-GRAPH**
  - Combines all proximity algorithms for a given network

- **MULTIPLEX**
  - Integrates all networks into a multiplex framework for Random Walk with Restart

These directories contain:
- **Ageing_Prediction**
  - Predictions only
  - Datasets are integrated directly within the prediction code

---

#### 2.2.2 Topological_data

- Random Walk with Restart (RWR) results
- Organized by:
  - ARCs
  - ARDs
- Includes results for:
  - PPI
  - KEGG
  - COX90
  - COX95
  - joint summaries and multiplex analyses

---

### 2.3 Permutations

Permutation-based statistical analyses.

Subdirectories include:
- **Coexpression**
- **KEGG_Hierarchy**
- **Network_Associations**
- **Specificity**

Each subdirectory contains:
- complete results stored as `rds` objects
- Excel files with summary statistics (mean, SD, p-values)

Descriptions:

- **Coexpression**
  - Tests whether ageing and high-pleiotropy genes show non-random inter-set coexpression with disease genes

- **KEGG_Hierarchy**
  - Tests whether ageing and high-pleiotropy genes occupy different hierarchical positions in KEGG signaling cascades
  - Based on distances to pathway roots and leaves

- **Network_Associations**
  - Permutation tests for:
    - minimum gene–disease distance
    - number of ARC-interactors
    - RWR scores

- **Specificity**
  - Tests whether ageing and high-pleiotropy genes exhibit distinct tissue specificity compared to random genes

---

### 2.4 Specificity_and_coexpression

Auxiliary data used for visualization and intermediate computations.

Subdirectories:
- **Coexpression**
- **Specificity**

Contains:
- intermediate R objects used for calculations
- Excel files used to generate:
  - boxplots
  - heatmaps
  - tissue-level coexpression and specificity visualizations

---

## Directory structure (overview)

```
Retrieved/
├── Genes_and_diseases
│ ├── Diseases
│ ├── Genes_List
│ ├── HAGR
│ ├── OMA
│ ├── Ranges
│ └── Showcase
├── Network_Sources
│ ├── Coexpression
│ ├── KEGG
│ └── PPI
├── Specificity
└── Ukbb

Generated/
├── Genes_Enrichment
├── Networks_and_predictions
│ ├── Networks
│ │ ├── PPI
│ │ ├── KEGG
│ │ ├── COX90
│ │ ├── COX95
│ │ ├── PER-ALGORITHM
│ │ ├── PER-GRAPH
│ │ └── MULTIPLEX
│ ├── Topological_data
│ │ └── RWR
├── Permutations
│ ├── Coexpression
│ ├── KEGG_Hierarchy
│ ├── Network_Associations
│ └── Specificity
└── Specificity_and_coexpression
├── Coexpression
└── Specificity
```
