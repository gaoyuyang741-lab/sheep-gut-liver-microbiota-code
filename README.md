[text](<Author Checklist - Full.pdf>)# Analysis code for sheep gastrointestinal microbiota and host metabolic integration study

This repository contains the R scripts used for the computational analyses and figure generation in the manuscript:

**[Rumen and colon-centered microbiota–liver interactions associated with bile acid-related metabolic regulation in sheep ]**

## Overview

This study integrated ruminal, ileal, and colonic microbiota data with hepatic gene expression profiles, blood biochemical traits, and tail fat phenotypes to investigate compartment-specific and cross-compartment microbiota-host associations in sheep.

The scripts in this repository document the computational workflow used for microbial diversity analysis, differential abundance analysis, predicted functional profiling, cross-compartment coordinated genera analysis, cross-omics correlation analysis, sparse generalized canonical correlation analysis, mediation analysis, and figure generation.

## Repository structure

```text
.
├── 01_alpha_diversity_analysis.R
├── 02_alpha_diversity_figure_generation.R
├── 03_beta_diversity_analysis.R
├── 04_beta_diversity_figure_generation.R
├── 05_ancombc2_differential_taxa.R
├── 06_ancombc2_figure_generation.R
├── 07_picrust2_functional_analysis.R
├── 08_picrust2_figure_generation.R
├── 09_cross_compartment_coordinated_genera_analysis.R
├── 10_cross_compartment_coordinated_genera_figure_generation.R
├── 11_pairwise_correlation_analysis.R
├── 12_pairwise_figure_generation.R
├── 13_sgcca_analysis.R
├── 14_sgcca_postprocessing.R
├── 15_sgcca_figure_generation.R
├── 16_mediation_analysis.R
├── 17_mediation_figure_generation.R
└── README.md
```

## Script description

| Script | Description |
|---|---|
| `01_alpha_diversity_analysis.R` | Performs alpha-diversity analysis across rumen, ileum, and colon samples. |
| `02_alpha_diversity_figure_generation.R` | Generates alpha-diversity figures. |
| `03_beta_diversity_analysis.R` | Performs beta-diversity analysis, including distance-based ordination and statistical testing. |
| `04_beta_diversity_figure_generation.R` | Generates beta-diversity PCoA figures. |
| `05_ancombc2_differential_taxa.R` | Performs differential abundance analysis using ANCOM-BC2 or related statistical procedures. |
| `06_ancombc2_figure_generation.R` | Generates figures for differential taxa and spatial pattern classification. |
| `07_picrust2_functional_analysis.R` | Performs predicted functional profiling and downstream pathway-level analysis based on PICRUSt2 outputs. |
| `08_picrust2_figure_generation.R` | Generates figures for predicted functional profiles. |
| `09_cross_compartment_coordinated_genera_analysis.R` | Identifies and classifies cross-compartment coordinated genera across rumen, ileum, and colon. |
| `10_cross_compartment_coordinated_genera_figure_generation.R` | Generates figures for cross-compartment coordinated genera. |
| `11_pairwise_correlation_analysis.R` | Performs pairwise cross-block correlation analyses among microbiota, hepatic gene expression, and host traits. |
| `12_pairwise_figure_generation.R` | Generates figures for pairwise correlation results. |
| `13_sgcca_analysis.R` | Performs sparse generalized canonical correlation analysis using multi-block datasets. |
| `14_sgcca_postprocessing.R` | Performs post-processing of sGCCA outputs and extracts representative multi-block association chains. |
| `15_sgcca_figure_generation.R` | Generates figures for sGCCA results. |
| `16_mediation_analysis.R` | Performs mediation analysis for microbiota-hepatic gene-host trait association chains. |
| `17_mediation_figure_generation.R` | Generates figures for mediation analysis results. |

## Input and output directories

The scripts use relative paths to improve portability. Input files should be placed in the corresponding `data/` subdirectories, and output files will be saved to `results/` or `figures/`.

The expected directory structure is:

```text
.
├── data/
│   ├── diversity/
│   │   ├── alpha/
│   │   │   └── input/
│   │   └── beta/
│   │       └── input/
│   ├── feature_taxa/
│   │   └── input/
│   ├── picrust2/
│   │   └── input/
│   ├── coordinated_genera/
│   │   └── input/
│   ├── pairwise/
│   │   └── input/
│   ├── sgcca/
│   │   └── input/
│   └── mediation/
│       └── input/
│
├── results/
│   ├── diversity/
│   ├── feature_taxa/
│   ├── picrust2/
│   ├── coordinated_genera/
│   ├── pairwise/
│   ├── sgcca/
│   └── mediation/
│
└── figures/
    ├── diversity/
    ├── feature_taxa/
    ├── picrust2/
    ├── coordinated_genera/
    ├── pairwise/
    ├── sgcca/
    └── mediation/
```

## Data availability

The raw sequencing data and processed datasets are described in the Data Availability Statement of the manuscript.

Input data files and large intermediate result tables are not included in this repository unless otherwise stated. The scripts use relative paths, and the required input files should be placed in the corresponding `data/` or `results/` directories before running the analyses.

## Software environment

The analyses were conducted in R.

Please install the required R packages before running the scripts. Major packages used in the analysis include:

```r
readxl
openxlsx
dplyr
tidyr
data.table
ggplot2
vegan
ANCOMBC
mixOmics
mediation
zCompositions
pheatmap
ComplexHeatmap
circlize
igraph
ggraph
ggrepel
```

The exact package list may vary among scripts. Please refer to the library-loading section of each script for details.

## Key analysis parameters

The main parameters used in the scripts include:

- Microbiome prevalence filtering threshold: `0.50` or `0.80`, depending on the analysis.
- Differential abundance analysis: ANCOM-BC2 or related non-parametric procedures.
- Multiple-testing correction: Benjamini-Hochberg or Holm correction, depending on the analysis.
- Correlation method: Spearman correlation.
- sGCCA number of components: `1`.
- sGCCA sparsity level: `0.30`.
- sGCCA network cutoff: `0.30`.
- Mediation analysis significance was evaluated using the criteria described in the manuscript.

Please refer to individual scripts for the exact parameters used in each analysis.

## Running the scripts

The scripts are numbered according to the main analytical workflow of the manuscript. In general, they should be run in the following order:

```text
01_alpha_diversity_analysis.R
02_alpha_diversity_figure_generation.R
03_beta_diversity_analysis.R
04_beta_diversity_figure_generation.R
05_ancombc2_differential_taxa.R
06_ancombc2_figure_generation.R
07_picrust2_functional_analysis.R
08_picrust2_figure_generation.R
09_cross_compartment_coordinated_genera_analysis.R
10_cross_compartment_coordinated_genera_figure_generation.R
11_pairwise_correlation_analysis.R
12_pairwise_figure_generation.R
13_sgcca_analysis.R
14_sgcca_postprocessing.R
15_sgcca_figure_generation.R
16_mediation_analysis.R
17_mediation_figure_generation.R
```

Before running the scripts, set the working directory to the root directory of this repository.

## Citation

If using this code, please cite the associated manuscript and the archived repository DOI:

**[]**

## Contact

For questions regarding this code, please contact:

**Yuyang Gao**  
Email: **gaoyuyang741@gmail.com**
