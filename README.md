# Dynamic Structural Profiling of PINK1 Mutations (T313M and L347P)

**Computational workflow for the identification, annotation, protein–protein docking, and molecular dynamics analysis of pathogenic PINK1 variants associated with Autosomal Recessive Parkinson's Disease (ARPD).**

---

## Overview

This repository contains the Python scripts and input files used in the computational analyses accompanying the manuscript:

> **Dynamic Structural Profiling of PINK1 Mutations (T313M and L347P) Reveals a Critical Phospho-Ser65 Ubiquitin Recognition Point in Mitophagy-Mediated Autosomal Recessive Parkinson's Disease**

The study investigates how the pathogenic **T313M** and **L347P** mutations alter the interaction between **PINK1** and **ubiquitin**, affecting substrate recognition, phosphorylation-dependent signaling, and mitophagy in Autosomal Recessive Parkinson's Disease (ARPD).

The computational workflow includes variant retrieval, functional annotation, clinical interpretation, protein structure preparation, protein–protein docking, molecular dynamics simulations, and structural analyses.

---

## Repository Structure

```text
PINK1-ARPD-Computational-Analysis/
│
├── README.md
├── LICENSE
│
├── scripts/
│   ├── 01_fetch_variant_ids.py
│   ├── 02_fetch_polyphen_sift_scores.py
│   ├── 03_fetch_clinical_data.py
│   └── data_visualizer.py
│
└── input/
    ├── rs_ids.csv
    ├── variant_clinical_input.csv
    ├── variant_clinical_input_filtered.csv
    ├── WT_PINK1.pdb
    ├── T313M.pdb
    ├── L347P.pdb
    ├── Ubiquitin.pdb
    └── pUb_SEP.pdb
```

---

## Computational Workflow

```text
NCBI dbSNP
      │
      ▼
Retrieve PINK1 Missense Variants
      │
      ▼
Ensembl Variant Effect Predictor (VEP)
      │
      ▼
PolyPhen-2 & SIFT Functional Prediction
      │
      ▼
ClinVar Clinical Annotation
      │
      ▼
Selection of Pathogenic Variants
      │
      ▼
Protein Structure Preparation
      │
      ▼
Protein–Protein Docking (PINK1–Ubiquitin)
      │
      ▼
Molecular Dynamics Simulations
      │
      ▼
Trajectory Analysis
      │
      ├── RMSD
      ├── RMSF
      ├── Radius of Gyration (Rg)
      ├── SASA
      ├── Hydrogen Bonds
      ├── DSSP
      ├── PCA
      ├── DCCM
      └── MM/PBSA
      │
      ▼
Structural Interpretation
```

---

## Scripts

### 1. `01_fetch_variant_ids.py`

Retrieves human **PINK1** missense variants from the **NCBI dbSNP** database using Selenium.

**Input**

None

**Output**

```text
rs_ids.csv
```

---

### 2. `02_fetch_polyphen_sift_scores.py`

Annotates variants using the **Ensembl Variant Effect Predictor (VEP)** REST API.

Retrieved information includes:

* Amino acid substitution
* PolyPhen-2 score
* PolyPhen-2 prediction
* SIFT score
* SIFT prediction

**Input**

```text
input/rs_ids.csv
```

**Output**

```text
variant_clinical_input.csv
```

---

### 3. `03_fetch_clinical_data.py`

Retrieves clinical annotations from **ClinVar**.

Retrieved information includes:

* ClinVar Variation ID
* Clinical significance
* Review status
* Protein change
* Gene symbol

**Input**

```text
input/variant_clinical_input_filtered.csv
```

**Output**

```text
clinvar_required_fields.csv
```

---

### 4. `data_visualizer.py`

Python utility for generating publication-quality plots from molecular dynamics simulation outputs.

Supported analyses include:

* Root Mean Square Deviation (RMSD)
* Root Mean Square Fluctuation (RMSF)
* Radius of Gyration (Rg)
* Solvent Accessible Surface Area (SASA)
* Hydrogen Bond Analysis
* DSSP Secondary Structure
* Dynamic Cross-Correlation Matrix (DCCM)
* Essential RMSF (eRMSF)
* Principal Component Analysis (PCA)
* Free Energy Landscape (FEL)
* MM/PBSA Energy Analysis
* PLIP Interaction Analysis

Example:

```bash
python data_visualizer.py --folders WT T313M L347P --mode rmsd
```

---

## Input Files

| File                                  | Description                                       |
| ------------------------------------- | ------------------------------------------------- |
| `rs_ids.csv`                          | PINK1 missense variant identifiers                |
| `variant_clinical_input.csv`          | Functional annotation generated using Ensembl VEP |
| `variant_clinical_input_filtered.csv` | Filtered variants used for ClinVar annotation     |
| `WT_PINK1.pdb`                        | Wild-type PINK1 structure                         |
| `T313M.pdb`                           | T313M mutant PINK1 structure                      |
| `L347P.pdb`                           | L347P mutant PINK1 structure                      |
| `Ubiquitin.pdb`                       | Human ubiquitin structure                         |
| `pUb_SEP.pdb`                         | Ser65-phosphorylated ubiquitin structure          |

---

## Installation

Clone the repository:

```bash
git clone https://github.com/<your-github-username>/PINK1-ARPD-Computational-Analysis.git
```

Move into the repository:

```bash
cd PINK1-ARPD-Computational-Analysis
```

Install the required Python packages:

```bash
pip install pandas requests selenium matplotlib numpy
```

---

## Software Requirements

* Python 3.9 or later
* Google Chrome
* ChromeDriver

---

## Public Databases

The scripts retrieve data from publicly available databases:

* NCBI dbSNP
* Ensembl Variant Effect Predictor (VEP)
* ClinVar
* UniProt
* AlphaFold Protein Structure Database
* RCSB Protein Data Bank

---

## Expected Outputs

The workflow generates:

* Variant identifier datasets
* Functional annotation tables
* Clinical annotation tables
* Publication-quality molecular dynamics plots
* Structural analysis figures

---

## Data Availability

This repository contains the Python scripts and input files used in the computational analyses supporting the associated manuscript.

Large molecular dynamics trajectory files (`.xtc`, `.trr`, `.edr`) and intermediate simulation outputs are not included because of their size. These files are available from the corresponding author upon reasonable request.

---

## Citation

If you use the scripts, input files, or workflow provided in this repository, please cite the associated publication.

---

## Authors

**Deborah Vincent**
Department of Biotechnology
School of Bio Sciences and Technology
Vellore Institute of Technology
Vellore, Tamil Nadu, India

### Corresponding Author

**Prof. C. Sudandiradoss**
Department of Biotechnology
School of Bio Sciences and Technology
Vellore Institute of Technology

Email: [csudandiradoss@vit.ac.in](mailto:csudandiradoss@vit.ac.in)

---

## License

This project is licensed under the **MIT License**.

See the **LICENSE** file for complete details.

---

## Acknowledgements

The authors gratefully acknowledge the developers and maintainers of the public databases and open-source software used in this work, including **NCBI**, **Ensembl**, **ClinVar**, **UniProt**, **AlphaFold Protein Structure Database**, **RCSB Protein Data Bank**, **HADDOCK**, **GROMACS**, and the Python scientific computing community.
