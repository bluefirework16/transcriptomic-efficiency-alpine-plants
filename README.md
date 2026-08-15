# Transcriptomic efficiency shapes growth-defense trade-offs in alpine desert plants under climate change

This repository contains the R scripts and source data used for the statistical analyses and figure generation associated with the manuscript:

**Transcriptomic efficiency shapes growth-defense trade-offs in alpine desert plants under climate change**

The study investigates physiological and transcriptomic responses of two co-dominant alpine desert shrubs, *Ajania tibetica* and *Krascheninnikovia ceratoides* (referred to as *Ceratoides compacta* in the associated manuscript), to experimental warming and watering in the alpine desert of Xizang, China.

## Repository structure

```text
R/
├── _helpers.R
├── 01_fig1_lmm_rda.R
├── 02_fig2_ptdi.R
├── 03_fig3_wgcna.R
├── 04_fig4_functional_associations.R
├── 05_fig5_hub_heatmaps.R
├── 06_fig6_climate.R
└── 10_prepare_deseq2_gene_sets.R

data/
└── Additional_file_2_Source_data.xlsx
```

The scripts are organized according to the order of the main figures in the manuscript.

## Data availability

Processed RNA-seq gene-expression matrices and functional annotations used in the transcriptomic analyses are publicly available from the OMIX repository of the National Genomics Data Center, China National Center for Bioinformation:

- *Ajania tibetica*: **OMIX010412**
- *Krascheninnikovia ceratoides*: **OMIX010413**

The OMIX datasets contain gene-level raw count matrices, FPKM-normalized expression values, and functional annotation information.

Source data underlying the main and supplementary figures are provided in:

`data/Additional_file_2_Source_data.xlsx`

The processed expression matrices are not duplicated in this repository.

## Analysis workflow

The main analysis scripts are organized as follows:

- `R/01_fig1_lmm_rda.R`  
  Linear mixed-effects models and redundancy analysis used to evaluate physiological and multivariate responses to warming and watering treatments.

- `R/02_fig2_ptdi.R`  
  Calculation and analysis of the physiological-transcriptomic decoupling index (PTDI), including transcriptomic and physiological multivariate distances.

- `R/03_fig3_wgcna.R`  
  Weighted gene co-expression network analysis (WGCNA) and module-trait association analyses.

- `R/04_fig4_functional_associations.R`  
  Functional gene-set analyses and associations between transcriptomic functional components and physiological traits.

- `R/05_fig5_hub_heatmaps.R`  
  Visualization of expression patterns for selected hub genes.

- `R/06_fig6_climate.R`  
  Analysis and visualization of long-term climatic trends at the study site.

- `R/10_prepare_deseq2_gene_sets.R`  
  Preparation of DESeq2-derived gene sets used in downstream transcriptomic analyses.

Shared functions used across scripts are stored in:

`R/_helpers.R`

## Input data

The source-data workbook included in this repository is:

`data/Additional_file_2_Source_data.xlsx`

For analyses requiring the full RNA-seq expression matrices, users should download the corresponding processed datasets from OMIX:

- OMIX010412 for *Ajania tibetica*
- OMIX010413 for *Krascheninnikovia ceratoides*

The downloaded expression files should be placed in the local directory specified in the relevant analysis scripts before running the transcriptomic analyses.

## Software

All analyses were conducted in R.

Major R packages used across the workflow include packages for:

- linear mixed-effects modelling
- multivariate analysis
- DESeq2-based differential expression analysis
- weighted gene co-expression network analysis
- data manipulation
- statistical visualization

Exact package requirements are specified within the corresponding R scripts.

## Reproducibility

The scripts are intended to reproduce the analyses and figures reported in the associated manuscript when used together with:

1. the source-data workbook provided in this repository, and
2. the processed gene-expression datasets archived in CNCB OMIX.

Users should maintain the directory structure described above when running the scripts.

## License

The analysis code in this repository is released under the MIT License.

## Citation

If you use the code or associated datasets, please cite the accompanying article and the archived release of this repository.

A permanent DOI for the archived code release will be provided through Zenodo.
