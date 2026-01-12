# Pleiotropy and Disease Interactors

This repository contains the **analysis code** supporting the study  
**“Pleiotropy and Disease Interactors: The Dual Nature of Genes Linking Ageing and Ageing-related Diseases.”**

The project investigates how ageing-related genes and age-related disease genes are interconnected across biological networks, and how pleiotropy emerges from these interactions.

---

## Repository structure

The repository is organized into the following main components:

- **`Code/`**  
  Contains all scripts used for data processing, network construction, statistical analysis, machine learning, and figure generation.

- **`docs/`**  
  Provides detailed documentation describing the analysis pipeline, code organization, and data structure.

- **`README.md`**  
  High-level overview of the project, data access, and documentation.

---

## Documentation

Detailed technical documentation is provided in the `docs/` directory:

- **`CODE_OVERVIEW.md`**  
  High-level description of the analysis codebase, organized by major analytical modules.

- **`PIPELINE.md`**  
  Step-by-step description of the complete analysis pipeline, from data retrieval to enrichment analyses.

- **`DATA_OVERVIEW.md`**  
  Conceptual overview of retrieved and generated data types used in this study.

- **`DATA_ORGANIZATION.md`**  
  Detailed directory-level description of all retrieved and generated datasets, including their biological meaning and provenance.

Readers interested in specific analysis stages are encouraged to start with `PIPELINE.md` and `DATA_OVERVIEW.md`.

---

## Data availability

Due to GitHub file size limitations, **datasets are not hosted directly in this repository**.

All retrieved source data and generated results associated with this study are publicly available through the **Synapse data-sharing platform** under the project:

**Pleiotropy and disease interactors**  
**Synapse Project ID:** `syn72037936`  
**Direct link:** https://www.synapse.org/#!Synapse:syn72037936

Data are organized in structured subdirectories corresponding to retrieved datasets and generated analysis outputs.

Access to the data is subject to Synapse terms of use and any applicable access restrictions.

---

## Reproducibility

The scripts in the `Code/` directory are designed to work directly with the datasets hosted on Synapse.  
Paths and parameters are configurable within individual scripts.

Large intermediate files and computationally intensive results are not included in this repository but are available via Synapse to ensure full reproducibility of all reported analyses.

---

## Data sources

The analyses integrate data from multiple public resources, including:

- UK Biobank
- GenAge (Human and Model Organisms)
- GeneFriends
- BioGRID
- KEGG
- GTEx
- GWAS summary statistics (Donertas et al.)

---

## Notes

- Large datasets and intermediate files are intentionally excluded from GitHub.
- Synapse serves as the authoritative data repository for this study.
- This repository focuses on transparency, reproducibility, and clarity of the analytical workflow.
