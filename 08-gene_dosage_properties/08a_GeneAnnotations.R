##### D2N - Qualitative and quantitative gene annotations #####

# Qualitative (GA): 
# (0) Gene, transcript and biotype
# (1) GWAS blood traits associated genes
# (2) Disease genes (OMIM)
# (3) Houskeeping genes 
# (4) TFs

# Quantitative (GQ):
# (1) Constraint (gnomAD: pLI, oe_lof_upper, mis_z)
# (2) Dosage sensitivity (Collins et al. 2022)
# (3) Enhancer Domain Scores (Goldstein 2020)
# (4) TFs (Number of TFs each gene is regulated by)
# (5) PPIs (STRINGdb)



# JDE, Nov 2022
# Last modified: August 2023


## Libraries
library(data.table)
library(Seurat)
library(biomaRt)
library(foreach)
library(doParallel); registerDoParallel(detectCores()/2)
library(ggplot2); theme_set(theme_classic())
library(stringr)
library(cowplot)


## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/08-gene_dosage_properties/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


## Load data

# D2N dataset
d2n <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_SeuratObj_PostQC_dCas9_NewMD.RDS")
d2n_genes <- rownames(d2n)[!grepl("CRISPR", rownames(d2n))]
dosage_genes = c("GFI1B", "NFE2", "MYB", "TET2")


# Add ENSEMBL gene id
ensembl <- biomaRt::useEnsembl(
  biomart = "ensembl", 
  dataset = "hsapiens_gene_ensembl",
  host = "https://www.ensembl.org"
)
symbols <- as.data.table(biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "ensembl_peptide_id", "gene_biotype", "chromosome_name"),
  filters = "external_gene_name",
  values = c(d2n_genes, "CFAP210"),
  mart = ensembl
)) 
symbols[ensembl_gene_id == "ENSG00000154479", external_gene_name := "CCDC173"]
symbols <- unique(symbols)

GA.0 <- setnames(unique(symbols[, .(external_gene_name, gene_biotype)]), c("gene", "gene_biotype"))
GA.0[, value := ifelse(gene_biotype == "lncRNA", T, F)]
GA.0[, gene_set := "lncRNA"]
GA.0[, dataset := "Biotype"]

## Save to use in the future
saveRDS(file = file.path(processed_data_dir, "GA.0_IDsAnnotations.RDS"), GA.0[, .(gene,value, gene_set, dataset)])

## (1) GWAS genes - UKB study closest K562 expressed gene

d2n_lncRNAs_v <- unique(symbols[gene_biotype == "lncRNA", external_gene_name])
GWAS_dt <- rbind(fread(file.path(data_dir, "MorrisScience2023_SupplementaryTable_refGeneAnnotations.txt"))[, .(gene, Platelets, RBCs, Reticulocytes, WBCs)],
                 data.table(gene = d2n_lncRNAs_v,
                            Platelets = rep(F, length(d2n_lncRNAs_v)),
                            RBCs = rep(F, length(d2n_lncRNAs_v)), 
                            Reticulocytes = rep(F, length(d2n_lncRNAs_v)), 
                            WBCs = rep(F, length(d2n_lncRNAs_v))) )

GA.1 <- melt.data.table(GWAS_dt[gene %in% d2n_genes, ], id.vars = "gene", variable.name = "GWAS_trait", value.name = "value")
GA.1[, gene_set := paste0(GWAS_trait, "(UKB GWAS)")]
GA.1[, dataset := "UKB"]


## Save to use in the future
saveRDS(object = GA.1[, .(gene,value, gene_set, dataset)], file = file.path(processed_data_dir, "GA.1_GWASgenesUKB.RDS"))


## (2) Disease genes - OMIM database
GM2_raw <- fread(file.path(data_dir, "OMIM_genemap2.txt"))

GM2 <- GM2_raw[!(Phenotypes == "" | `Ensembl Gene ID` == ""), ]
sum(unique(GM2$`Ensembl Gene ID`) %in% symbols$ensembl_gene_id)

GM2_D2N <- GM2[`Ensembl Gene ID` %in% symbols$ensembl_gene_id, ]
omim_genes <- unique(symbols[ensembl_gene_id %in% GM2_D2N$`Ensembl Gene ID`, external_gene_name])


GA.2 <- data.table(gene = d2n_genes,
                   gene_set = "Disease_gene(OMIM)", 
                   value = d2n_genes %in% omim_genes,
                   dataset = "OMIM")


## Save to use in the future
saveRDS(file = file.path(processed_data_dir, "GA.2_DiseaseGenesOMIM.RDS"), GA.2)


## (3) Housekeeping genes

hkg_long <- fread(file.path(data_dir, "EisenbergTrendsGen2013_HK_genes.txt"), header = F,  stringsAsFactors = F, fill = T)$V1
hkg <- unique(hkg_long)[!(is.na(unique(hkg_long)))]

hkg_symbols <- as.data.table(biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = hkg,
  mart = ensembl
))
hkg_symbols[, d2n := ifelse(ensembl_gene_id %in% unique(symbols$ensembl_gene_id), T, F)]
hkg_symbols[, d2n_genename := ifelse(external_gene_name %in% d2n_genes, T, F)]

hkg_symbols[, .N, d2n]

GA.3 = data.table(gene = d2n_genes,
                   gene_set = "Houskeeping", 
                   value = d2n_genes %in% hkg_symbols[d2n == T, external_gene_name],
                   dataset = "Eisenberg_2013")

## Save to use in the future
saveRDS(file = file.path(processed_data_dir, "GA.3_HousekeepingHsiao.RDS"), GA.3)


## (4) TFs - Garcia-Alonso 2019 

# Load Garcia Alonso dataset
TF_Ga_dt <- unique(rbind(fread(file = file.path(data_dir, "GarciaAlonso_TableS3_Normal.csv"))[, .(TF, target)],
                         fread(file = file.path(data_dir, "GarciaAlonso_TableS4_Cancer.csv"))[, .(TF, target)]))
TF_Ga <- unique(TF_Ga_dt$TF)

GA.4 <- data.table(gene = d2n_genes,
               value = d2n_genes %in% TF_Ga,
               gene_set = "TF",
               dataset = "Garcia-Alonso_2019")

## Save to use in the future
saveRDS(file = file.path(processed_data_dir, "GA.4_TFs.RDS"), GA.4)



## Quantitative

# (1) Gene body constraint
gnomAD_metrics_dt <- fread(file.path(data_dir, "gnomad.v2.1.1.lof_metrics.by_gene.txt"))[gene_id %in% symbols$ensembl_gene_id,]
gnomAD_metrics_dt[, gene_gnomAD := gene]

gene_to_replace = c("MTERFD3", "FAM175A", "MESDC1", "LINC00998")
replace_names = c("MTERF2", "ABRAXAS1", "TLNRD1", "SMIM30")
gnomAD_metrics_dt[gene_gnomAD %in% gene_to_replace, gene := replace_names[match(gene_gnomAD, gene_to_replace)]]

GQ.1 <- melt.data.table(rbind(gnomAD_metrics_dt[, .(gene, pLI, oe_lof_upper, mis_z)],
                              data.table(gene = d2n_genes[!(d2n_genes %in% gnomAD_metrics_dt$gene)],
                                         pLI = NA,
                                         oe_lof_upper = NA,
                                         mis_z = NA)),
                        id.vars = "gene", variable.name = "metric")
GQ.1[, dataset := "gnomAD"]

## Save to use in the future
saveRDS(file = file.path(processed_data_dir, "GQ.1_Constraint.RDS"), GQ.1)



## (2) Dosage sensitivity
collins22 <- fread(file.path(data_dir, "CollinsCell2022_haplosufficiency.txt"))[, gene := Gene]

gene_to_replace = c("MTERFD3", "FAM175A", "MESDC1", "LINC00998", "FAM27A")
replace_names = c("MTERF2", "ABRAXAS1", "TLNRD1", "SMIM30", "FAM27C")

collins22[Gene %in% gene_to_replace, gene := replace_names[match(Gene, gene_to_replace)]]

GQ.2 <-  melt.data.table(rbind(collins22[gene %in% d2n_genes, .(gene, pHaplo, pTriplo)],
                               data.table(gene = d2n_genes[!(d2n_genes %in% collins22$gene)],
                                          pHaplo = NA,
                                          pTriplo = NA)),
                         id.vars = "gene", variable.name = "metric")
GQ.2[, dataset := "Collins et al. 2022"]


# Save
saveRDS(file = file.path(processed_data_dir, "GQ.2_DosageSensitivity.RDS"), GQ.2)





## (3) Enhancer Domain Scores
Wang20_dt <- merge.data.table(fread(file.path(data_dir, "WangAJHG2020_enhancer_scores.txt")), unique(symbols[, .(ensembl_gene_id, external_gene_name)]), by.x = "GeneSymbol", by.y = "ensembl_gene_id")
Wang20_dt[, gene := external_gene_name]

GQ.3 <- melt.data.table(rbind(Wang20_dt[gene %in% d2n_genes, .SD, .SDcols = colnames(Wang20_dt)[grepl("(Linking)|(EDS)|(^gene$)", colnames(Wang20_dt))]],
                              data.table(gene = d2n_genes[!(d2n_genes %in% Wang20_dt$gene)]), fill=TRUE),
                        id.vars = "gene", variable.name = "metric")
GQ.3[, dataset := "Wang & Goldstein 2020"]

# Save
saveRDS(file = file.path(processed_data_dir, "GQ.3_Enhancers.RDS"), GQ.3)


## (4) TFs
TF_raw <- fread(file.path(data_dir, "Minaeva2023_TF-target_PeaksMotiffs.tsv"))[ensembl_gene_id %in% unique(symbols$ensembl_gene_id),]

unique(symbols$ensembl_gene_id)[!(unique(symbols$ensembl_gene_id) %in% unique(TF_raw$ensembl_gene_id) )] #FAM27C not in there

TF_M2 <- merge.data.table(TF_raw[is_method_2 == T & tpm_total > 0.5, .(n_peaks = sum(n_peaks_method_2), n_motifs = sum(n_motifs_method2, na.rm = T)), .(ensembl_gene_id, tf)],
                          setnames(unique(symbols[, .(ensembl_gene_id, external_gene_name)]), c("ensembl_gene_id", "gene")),
                          by = "ensembl_gene_id")

TF_T_dt <- cbind(data.table(gene = d2n_genes), foreach(dg = c("GFI1B", "MYB", "NFE2"), .combine = cbind) %do% {
  new_vars = paste0(c("n_peaks", "n_motifs"), "_", dg)
  dt <- setnames(TF_M2[tf == dg, ], old = c("n_peaks", "n_motifs"), new = new_vars)
  dt2 <- merge.data.table(unique(symbols[, .(ensembl_gene_id, external_gene_name)]), dt, by = "ensembl_gene_id", all.x = T)
  dt2[is.na(gene), new_vars[1] := 0]
  dt2[is.na(gene), new_vars[2] := 0]
  dt2[match(d2n_genes, external_gene_name), .SD, .SDcols = colnames(dt2)[grepl("n_", colnames(dt2))]]
})


GQ.4 <-  melt.data.table(TF_T_dt, id.vars = "gene", variable.name = "metric")
GQ.4[, dataset := "Minaeva in prep. 2023"]


# Save
saveRDS(file = file.path(processed_data_dir, "GQ.4_TFs.RDS"), GQ.4)



## (5) PPIs
# PPIs annotations

# PPI data
STRING_info_dt <- fread(input = file.path(data_dir, "STRINGdb_protein.info.v11.5.txt.gz"), col.names = c("string_protein_id", "preferred_name", "protein_size", "annotation"))
STRING_info_dt[, ensembl_peptide_id := gsub("9606\\.", "", string_protein_id)]

# Subset d2n genes from STRING
STRING_info_d2n_dt <- STRING_info_dt[ensembl_peptide_id %in% symbols$ensembl_peptide_id | preferred_name %in% symbols$external_gene_name, ]

# Load PPI datasets
raw_string_physical <- fread(input = file.path(data_dir, "STRINGdb_protein.physical.links.detailed.v11.5.txt.gz"))
string_physical_temp <- fread(input = file.path(data_dir, "STRINGdb_protein.physical.links.detailed.v11.5.txt.gz"))[protein1 %in% STRING_info_d2n_dt$string_protein_id & protein2 %in% STRING_info_d2n_dt$string_protein_id,]

string_physical <- data.table::merge.data.table(data.table::merge.data.table(string_physical_temp, 
                                                                             data.table::setnames(STRING_info_d2n_dt[, .(string_protein_id, preferred_name)], c("protein1", "gene1")), 
                                                                             by = "protein1"), 
                                                data.table::setnames(STRING_info_d2n_dt[, .(string_protein_id, preferred_name)], c("protein2", "gene2")), 
                                                by = "protein2")

# Two genes have genes have different external gene names - change them to the d2n nomenclature 
# FAM175A (STRING) -> ABRAXAS1 (d2n)
# MESDC1  (STRING) -> TLNRD1   (d2n)
# Replace multiple strings at a time
rep_str = c('MESDC1'='TLNRD1','FAM175A'='ABRAXAS1')

# Count number of PPI per gene (within the d2n genes or within the entire proteome)
num_PPI_dt <- merge.data.table(STRING_info_d2n_dt[, .(string_protein_id, preferred_name)],
                               merge.data.table(raw_string_physical[ protein1 %in% STRING_info_d2n_dt$string_protein_id, .N, .(protein1), ], 
                                                string_physical[protein1 %in% STRING_info_d2n_dt$string_protein_id, .N, .(protein1)], 
                                                by="protein1", all.x = T, suffixes = c("_int_wp", "_int_dn")),
                               by.x = "string_protein_id", by.y = "protein1", all.x = T)
num_PPI_dt[is.na(N_int_wp), N_int_wp := 0]
num_PPI_dt[is.na(N_int_dn), N_int_dn := 0]
num_PPI_dt$gene <- str_replace_all(num_PPI_dt$preferred_name, rep_str)



GQ.5 <- melt.data.table(rbind(setnames(num_PPI_dt[, .(gene, N_int_wp, N_int_dn)], old = c("N_int_dn", "N_int_wp"), new = c("Num_PPIs_targeted_network", "Num_PPIs_whole_proteome")),
                              data.table(gene = d2n_genes[!(d2n_genes %in% num_PPI_dt$gene)]), fill=TRUE),
                        id.vars = "gene", variable.name = "metric")
GQ.5[, dataset := "String DB"]

saveRDS(file = file.path(processed_data_dir, "GQ.5_PPIs.RDS"), GQ.5)



## (6) Mean cell type expression
# Annotations to cellypes
bm_annotations <- fread(file = file.path(data_dir, "Hay2018_SupplementaryTable_CellAssociations.txt"))
bm_annotations$cell_id <- gsub("\\.BM[1-9]", "", bm_annotations$Cell)
bm_annotations[, celltype := ClusterName]
bm_annotations[ClusterName == "CD34+ Eo/B/Mast", celltype := "CD34+ Eo-B-Mast"]
bm_annotations[ClusterName == "Naive T-cell", celltype := "Naive CD8 T-cell"]

# Get processed data from the bone marrow study
bm_mtx_raw <- fread(file = file.path(data_dir, "Hay2018_GE-BM-HCA-Donor-Avg.txt"))[uid %in% d2n_genes,]
bm_mtx <- rbind(bm_mtx_raw, data.table(uid = d2n_genes[!(d2n_genes %in% bm_mtx_raw$uid)]), fill = TRUE)
BM <- melt.data.table(bm_mtx, id.vars = "uid", variable.name = "donor_celltype")
BM$celltype <- gsub("[FM]-BM[1-9]__", "", BM$donor_celltype)
BM$celltype <- factor(BM$celltype, levels = unique(BM$celltype))
BM_avgdonor <- BM[,list(mean_expr_donors = mean(value)), .(uid, celltype)]
BM_avgdonor[, metric := celltype]

GQ.6 <- setnames(BM_avgdonor[!(grepl("^CD34", celltype)), .(uid, metric, mean_expr_donors)], 
                 old = c("uid", "mean_expr_donors"), 
                 new = c("gene", "value"))
GQ.6[, dataset := "Hay et al. 2018"]

saveRDS(file = file.path(processed_data_dir, "GQ.6_BoneMarrow.RDS"), GQ.6)


