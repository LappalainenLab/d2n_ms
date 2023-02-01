##### D2N - CIS genes FCs - sgRNA properties #####
# (1) FC differences per sgRNA type
# (2) FC vs On-target and Off-target sgRNA activity
# (3) Attenuated sgRNAs
#     - Relationship between dosage gene FC and sgRNA mismatch position?

# JDE, July 2022
# Last modified: February 2023



## Libs
library(data.table)
library(ggplot2)


## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/05-cis_fc/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


## Data
# !! temp until re-running all scripts from scratch
DE_dt <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_DGGuidePropertiesDE.RDS")
