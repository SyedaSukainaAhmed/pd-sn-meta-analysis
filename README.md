# Parkinson's disease substantia nigra transcriptome meta-analysis
Archived version: https://doi.org/10.5281/zenodo.24554168 
Code, results tables and figures for the preprint "Reproducible substantia nigra changes in Parkinson's disease track dopaminergic marker expression across samples" (Ahmed SS, 2026; bioRxiv link to be added).

## Data
Public GEO cohorts: GSE7621, GSE20163, GSE20292, GSE49036, GSE8397 (Affymetrix microarray) and GSE136666 (RNA-seq). Raw files are not included; script 01 and script 06 download them from GEO.

## How to run
Open R in the repository folder and run the scripts in order:

| Script | What it does | Main output |
|---|---|---|
| scripts/01_download_preprocess.R | Downloads CEL files, RMA-normalises each cohort | data/GSE*_results.rds |
| scripts/02_main_analysis.R | limma per cohort, rank aggregation, stability | results/RRA_*.csv, results/stability_*.csv |
| scripts/03_sensitivity.R | Outlier and scan-month scenarios, 69-gene core | results/core_genes.csv |
| scripts/04_overadjustment_checks.R | Attenuation, negative controls, within-group correlation | data/overadjustment_checks.rds |
| scripts/05_pathways.R | GO enrichment and reproducibility | results/pathway_reproducibility.csv |
| scripts/06_rnaseq_check.R | RNA-seq check in GSE136666 | results/RNAseq_GSE136666_DESeq2.csv |
| scripts/07_figures.R | Figures 1-7, sessionInfo.txt | figures/*.png |

scripts/00_helpers.R holds shared functions and is sourced by the others.

## Requirements
R 4.6.2 with Bioconductor packages: GEOquery, oligo, affy, limma, hgu133a.db, hgu133plus2.db, org.Hs.eg.db, clusterProfiler, DESeq2, EnsDb.Hsapiens.v86, and CRAN packages RobustRankAggreg, ggplot2. Script 01 needs about 3 GB of disk and some time. Full package versions are in sessionInfo.txt.

## Notes
Random-control analyses use fixed seeds, but exact counts can differ slightly between R and package versions. Each script prints the values reported in the preprint as comments ("expect ...") so results can be checked.

## License
MIT
