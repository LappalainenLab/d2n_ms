##### D2N - DE analyses #####
# (1) DE on NTC guides - filter out NTC with potential activity
# (2) DE on all genes
# (3) PCA on guide-to-FC or guide-to-expression matrices
# (4) Run Sceptre and calibration test

# JDE, July 2022
# last modification: October 2023


## Libraries
library(data.table)
library(Seurat)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(reshape2)
library(cowplot)
library(foreach)
library(doParallel); registerDoParallel(detectCores()/2)
library(ggfortify)



## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/04-differential_expression/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


col_crispr = list(CRISPRa = "#FF7F00", CRISPRi = "#377EB8")
dosage_genes = c("GFI1B", "NFE2", "MYB", "TET2")


### (1) DE on NTC and removal of NTCs with activity
# Seurat objects with filtered QCed cells
d2n <- readRDS("../007-QCandDemultiplex/processed_data/HTOs_GDOs_cDNA_dCas9_merged_postQC.RDS")

# Guides metadata
guides_features <- fread("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/001-library_design/processed_data/guides/D2N_ALLguides_features_sortID.csv")
guides_features[, id := gsub("_", ".", ID), ID]
guides_features[, guide_class := guide_category]
guides_features[guide_category == "weismanTSS", guide_class := "TSS"]
guides_features[guide_category == "distalCRE", guide_class := "distal CRE"]
guides_features[guide_category == "ntc", guide_class := "NTC"]

##### (1) FindMarkers approach 

# Add metadata to have extensive guide metadata
d2n$L_cell_barcode <- colnames(d2n)
new_md <- merge(d2n@meta.data, guides_features, by.x = "guide_call", by.y = "id", all.x = T, sort = F)
rownames(new_md) <- new_md$L_cell_barcode
d2n <- AddMetaData(object = d2n, metadata = new_md)
d2n$guide_crispr = paste0(d2n$short_ID, "-", d2n$cell_line)
Idents(d2n) <- "guide_crispr"


# List of all unique guide-crispr perturbations class
gene_crispr_v <- unique(d2n$guide_crispr)

# Get NT cells unique combinations
ntc_crispr_v <- gene_crispr_v[grep("NTC", gene_crispr_v)]

### (A) Run FindMarkers between all pairwise comparison of NT to find non-biased NT cells
ntc_crispr_comb_dt <- setnames(data.table(t(combn(sort(ntc_crispr_v),2))), c("ntc_1", "ntc_2"))

DE_NT_dt <- foreach(x = 1:nrow(ntc_crispr_comb_dt), .combine = rbind) %dopar% {
  id1 = ntc_crispr_comb_dt[x, ntc_1]
  id2 = ntc_crispr_comb_dt[x, ntc_2]
  FM_df <- FindMarkers(object = d2n, ident.1 = id1, ident.2 = id2,
                       min.pct = -Inf, logfc.threshold = -Inf, verbose = F)
  FM_dt <- data.table(gene = rownames(FM_df),
                      avg_log2FC = FM_df$avg_log2FC,
                      pval = FM_df$p_val,
                      pval_bfr = FM_df$p_val_adj,
                      pval_fdr = p.adjust(FM_df$p_val, method = "fdr"),
                      guide_1 = tstrsplit(id1, "-", keep = 1)[[1]],
                      guide_2 = tstrsplit(id2, "-", keep = 1)[[1]],
                      celltype_1 = tstrsplit(id1, "-", keep = 2)[[1]],
                      celltype_2 = tstrsplit(id2, "-", keep = 2)[[1]])
}

# Get summmary statistics about NTC pairwise comparisons
DE_NT_stats_dt <- DE_NT_dt[, .(ntests = .N,
                               mean_abs_log2FC = mean(abs(avg_log2FC)),
                               n_sigDE_fdr5 = sum(pval_fdr < 0.05),
                               n_sigDE_bfr5 = sum(pval_bfr < 0.05)),
                           .(guide_1, guide_2, celltype_1, celltype_2)]
DE_NT_stats_dt[, crispr_test := ifelse(celltype_1 != celltype_2, "CRISPRi \nvs.\nCRISPRa ", ifelse(celltype_1 == "CRISPRa", "CRISPRa \nvs.\nCRISPRa ", "CRISPRi \nvs.\nCRISPRi "))]
DE_NT_stats_dt[, any_DE := ifelse(n_sigDE_fdr5 > 1, "yes", "no")]
DE_NT_stats_dt[, guide_1 := factor(guide_1, levels = unique(DE_NT_stats_dt$guide_1))]
DE_NT_stats_dt[, guide_2 := factor(guide_2, levels = unique(DE_NT_stats_dt$guide_1))]


p1 <- ggplot(DE_NT_stats_dt, aes(x=crispr_test, y=n_sigDE_fdr5)) +
  geom_point(shape = 21, size=3, aes(alpha=any_DE, fill=guide_2, color=guide_1),  position = position_jitter(w = 0.3, h = 0)) +
  scale_alpha_manual(">1 DE gene?", values = c(0.2, 1)) +
  labs(x="NTC cells tested", y="number of DE genes") +
  scale_color_brewer("", palette = "Set1") +
  scale_fill_brewer("", palette = "Set1")
p1
ggsave("plots/20220721_DemuxHTOdCas9_FindMarkers_NTC_NumDEgenesCriprAvsI.pdf", p1, width = 4, height = 3.5)


p2 <- ggplot(DE_NT_stats_dt, aes(x=crispr_test, y=mean_abs_log2FC)) +
  geom_point(shape = 21, size=3, aes(alpha=any_DE, fill=guide_1, color=guide_2),  position = position_jitter(w = 0.3, h = 0)) +
  scale_alpha_manual("Any DE gene?", values = c(0.2, 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x="NTC cells tested", y="mean absolute log2FC") +
  scale_color_brewer("", palette = "Set1") +
  scale_fill_brewer("", palette = "Set1")
p2
ggsave("plots/20220721_DemuxHTOdCas9_FindMarkers_NTC_MeanAbsFCcriprAvsI.pdf", p2, width = 4, height = 3.5)

plot_dt <- merge.data.table(setnames(DE_NT_dt[pval_fdr < 0.05 & celltype_1 == "CRISPRa" & celltype_2 == "CRISPRa", .N, guide_1], c("guide", "N")), 
                            setnames(DE_NT_dt[pval_fdr < 0.05 & celltype_1 == "CRISPRa" & celltype_2 == "CRISPRa", .N, guide_2], c("guide", "N")),
                            by="guide", all.x = T, all.y = T)
plot_dt[, N := sum(N.x, N.y, na.rm = T), guide]
p1.2 <- ggplot(plot_dt, aes(x=guide, y=N)) +
  geom_bar(stat="identity", fill="grey70") +
  ylab("Num of DE genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1.2
ggsave("plots/20220805_DemuxHTOdCas9_FindMarkers_NTC_NumDE.pdf", p1.2, width = 2, height = 3.5)


# Create new identity categories in the metadata to have the NTC-core category for each cell line, excluding NTC_2-CRISPRa
d2n$ident <- d2n$guide_crispr
d2n$ident[grepl("NTC", d2n$guide_call)] <- gsub("_[1-5]", "", d2n$guide_crispr[grepl("NTC", d2n$guide_call)])
d2n$ident[d2n$guide_crispr == "NTC_2-CRISPRa"] <- NA


# Extract metadata and count how many cells we have for each guide-cell line
MD <- data.table(d2n@meta.data)
unique(d2n$ident)

plot_dt <- MD[, .N, .(ident, cell_line)]
plot_dt_median <- plot_dt[, .(N = median(N)), cell_line]
plot_dt_ntc <- plot_dt[grep("NTC", ident), N, cell_line]
p3 <- ggplot(plot_dt, aes(x=N)) +
  geom_histogram(bins = 50) +
  geom_vline(data = plot_dt_median, linetype=2, aes(xintercept=N)) +
  geom_vline(data = plot_dt_ntc, linetype=2, aes(xintercept=N), color="grey70") +
  facet_grid(cell_line ~ .) +
  geom_text(data = plot_dt_median, aes(label = N, y=13, x=N-25)) +
  geom_text(data = plot_dt_ntc, aes(label = paste0("NTC=", N), y=13, x=N-70), color="grey70") +
  labs(x="# cells per unique sgRNA")
p3
ggsave("plots/20220721_DemuxHTOdCas9_FindMarkers_DistNumCellsCond.pdf", p3, width = 4, height = 3.5)



### (2) DE on all genes
# Run FindMarkers between every guide_CRISPR condition vs well callibrated NT_CRISPR conditions
Idents(d2n) <- "ident"
guide_crispr_v <- gene_crispr_v[!(gene_crispr_v == "NTC_2-CRISPRa")]


DE_dt <- foreach(x = 1:length(guide_crispr_v), .combine = rbind) %dopar% {
  id1 = guide_crispr_v[x]
  cells1 = rownames(d2n@meta.data[d2n$guide_crispr == id1, ])
  id2 = paste0("NTC-", tstrsplit(id1, "-", keep = 2)[[1]])
  cells2 = rownames(d2n@meta.data[d2n$ident == id2 & !is.na(d2n$ident), ])
  gene = tstrsplit(id1, "_", keep = 1)[[1]]
  
  FM_df <- FindMarkers(object = d2n[["RNA"]], cells.1 = cells1, cells.2 = cells2, 
                       min.pct = -Inf, logfc.threshold = -Inf, verbose = F)
  FM_dt <- data.table(gene = rownames(FM_df), 
                      avg_log2FC = FM_df$avg_log2FC,
                      pval = FM_df$p_val,
                      pval_bfr = FM_df$p_val_adj,
                      pval_fdr = p.adjust(FM_df$p_val, method = "fdr"),
                      guide_1 = tstrsplit(id1, "-", keep = 1)[[1]],
                      guide_2 = tstrsplit(id2, "-", keep = 1)[[1]],
                      cell_line = tstrsplit(id1, "-", keep = 2)[[1]],
                      ncells_1 = length(cells1),
                      ncells_2 = length(cells2),
                      dosage_gene = gene,
                      dosage_gene_log2FC = ifelse(gene == "NTC", NA, FM_df[gene, "avg_log2FC"]),
                      dosage_gene_pval_fdr = ifelse(gene == "NTC", NA, FM_df[gene, "p_val_adj"]))
}



# Summary of DE analyses by guide
DE_stats_dt <- DE_dt[, .(ntests = .N, 
                         mean_log2FC = mean(avg_log2FC),
                         mean_abs_log2FC = mean(abs(avg_log2FC)),
                         sum_abs_log2FC = sum(abs(avg_log2FC)), 
                         n_sigDE_fdr5 = sum(pval_fdr < 0.05),
                         n_sigDE_bfr5 = sum(pval_bfr < 0.05)), 
                     .(guide_1, guide_2, cell_line, ncells_1, ncells_2, dosage_gene_log2FC, dosage_gene)]
DE_stats_dt[, dosage_gene := factor(dosage_gene, levels = c(dosage_genes, "NTC"))]

## FC depending on guide class
temp1 <- unique(DE_dt[!grepl("NTC", guide_1), .(cell_line, dosage_gene_log2FC, guide_1, dosage_gene_pval_fdr)])
temp1[, cis_gene := tstrsplit(guide_1, "_", keep = 1)]
temp2 <- setnames(DE_dt[grepl("NTC", guide_1) & gene %in% dosage_genes, .(cell_line, avg_log2FC, guide_1, pval_fdr, gene)], 
                  c("cell_line", "dosage_gene_log2FC", "guide_1", "dosage_gene_pval_fdr", "cis_gene"))

DE_guide_properties <- merge.data.table(rbind(temp1, temp2), guides_features, by.x = "guide_1", by.y = "short_ID")
DE_guide_properties[, cis_gene := factor(cis_gene, levels = c(dosage_genes, "NTC"))]
DE_guide_properties[, guide_category := factor(guide_category, levels = c("weismanTSS", "titration", "distalCRE", "attenuated", "ntc"))]
DE_guide_properties[, sig_fdr10 := ifelse(dosage_gene_pval_fdr < 0.1, T, F)]
DE_guide_properties[, guide_class := factor(guide_class, levels = c("TSS", "titration", "distal CRE", "attenuated", "NTC"))]


## Add genomic location info
TSS_gene <- fread("../001-library_design/processed_data/titration_genes_TSS_unique_site.bed", 
                  col.names = c("chr", "coord_tss", "coord_tss_1", "gene_tss", "score", "strand"))
TSS_gene[, gene := tstrsplit(gene_tss, "_", keep = 1)[[1]]]


DE_guide_properties_TSS <- merge.data.table(DE_guide_properties, 
                                            TSS_gene[, .(coord_tss, gene)], 
                                            by.x = "cis_gene", by.y = "gene", 
                                            all.x = T)
DE_guide_properties_TSS[, coord_guide := mean(c(coord_start, coord_end), na.rm = T), .(guide_1, cell_line, cis_gene)]
DE_guide_properties_TSS[, dist2tss := coord_guide - coord_tss, .(guide_1, cell_line, cis_gene)]


## SAVE DATA ##
saveRDS(DE_dt, file = "processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS")
saveRDS(DE_stats_dt, file = "processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_DGStatsDE.RDS")
saveRDS(DE_guide_properties_TSS, file = "processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_DGGuidePropertiesDE.RDS")

# Save the new Seurat object with new metadata for the future
saveRDS(d2n, file = "processed_data/after_seurat4.3/d2n_SeuratObj_PostQC_dCas9_NewMD.RDS")


### (3) PCA analyses
# On guide-by-FC matrix 

# Create a variable unique for each perturbation
DE_dt[, guide_crispr := paste0(guide_1, "-", cell_line)]

# dcast the data table with log2FC to obtain a guide-by-FC matrix
DE_dt_dcast <- dcast.data.table(DE_dt[, .(gene, avg_log2FC, guide_crispr, cell_line, dosage_gene, dosage_gene_log2FC)], 
                                formula = "guide_crispr + cell_line + dosage_gene + dosage_gene_log2FC ~ gene", value.var = "avg_log2FC")

DE_dt_dcast[,  prop_abs_dosage_gene_log2FC := abs(dosage_gene_log2FC)/max(abs(dosage_gene_log2FC), na.rm = T), dosage_gene]
DE_dt_dcast[is.na(dosage_gene_log2FC), prop_abs_dosage_gene_log2FC := 0.1]

# Run PCA
perturb_pca <- prcomp(DE_dt_dcast[, sort(unique(DE_dt$gene)),  with=FALSE], scale. = TRUE)


# Plot PCA scatterplots PC1-10
apply(matrix(1:10, ncol = 5), 2, function(R){
  p <- ggplot2::autoplot(perturb_pca, data = DE_dt_dcast[, !(colnames(DE_dt_dcast) %in% unique(DE_dt$gene)), with = FALSE], x=R[1], y=R[2], 
                         shape="cell_line", colour = "dosage_gene", alpha = "prop_abs_dosage_gene_log2FC", size = 3) +
    scale_shape_manual("cell line", values = c(16, 17)) +
    scale_color_manual("dosage gene", values = col_genes) +
    scale_alpha_continuous("Scaled cis\ngene abs(FC)")
  ggsave(paste0("plots/20220809_01_DemuxHTOdCas9_FindMarkers_PerturbatuonByGeneFCPCA_PC", R[1], "-", R[2], ".pdf"), p, width = 6, height = 4.5)
})

# Plot decay of variance explained per PC
plot_dt <- data.table(pct_var_expl = (perturb_pca$sdev^2)*100/sum(perturb_pca$sdev^2), 
                      pc = colnames(perturb_pca$x),
                      pc_num = 1:length(perturb_pca$sdev)) 

p <- ggplot(plot_dt[pc_num <= 25], aes(x = pc_num, y = pct_var_expl)) +
  geom_point() +
  geom_line() +
  labs(y = "% of variance explained", x = "PC")
p
ggsave("plots/20220809_02_DemuxHTOdCas9_FindMarkers_PerturbatuonByGeneFCPCA_PctVarExpl.pdf", p, width = 6, height = 4.5)



## PCA on guide-by-expression matrix 

# Iterate over guide_crispr and genes to get the mean expression matrix with same dimentions as the guide-by-FC matrix
genes <- colnames(DE_dt_dcast)[ colnames(DE_dt_dcast) %in% unique(DE_dt$gene) ]

Exp_dt_dcast <- cbind(
  DE_dt_dcast[, !(colnames(DE_dt_dcast) %in% unique(DE_dt$gene)), with = FALSE],
  setnames(
    data.table(
      do.call("rbind", lapply(DE_dt_dcast$guide_crispr, function(gc){
        rowMeans(d2n[["RNA"]]@data[genes, MD[guide_crispr == gc, L_cell_barcode] ])
      }))
    )
    , genes)
)

# Run PCA
expr_pca <- prcomp(Exp_dt_dcast[, sort(unique(DE_dt$gene)),  with=FALSE], scale. = TRUE)




# Plot PCA scatterplots PC1-10
apply(matrix(1:10, ncol = 5), 2, function(R){
  p <- ggplot2::autoplot(expr_pca, data = Exp_dt_dcast[, !(colnames(Exp_dt_dcast) %in% unique(DE_dt$gene)), with = FALSE], x=R[1], y=R[2], 
                         shape="cell_line", colour = "dosage_gene", alpha = "prop_abs_dosage_gene_log2FC", size = 3) +
    scale_shape_manual("cell line", values = c(16, 17)) +
    scale_color_manual("dosage gene", values = col_genes) +
    scale_alpha_continuous("Scaled cis\ngene abs(FC)")
  ggsave(paste0("plots/20220809_03_DemuxHTOdCas9_FindMarkers_PerturbatuonByGeneExprPCA_PC", R[1], "-", R[2], ".pdf"), p, width = 6, height = 4.5)
})

# Plot decay of variance explained per PC
plot_dt <- data.table(pct_var_expl = (expr_pca$sdev^2)*100/sum(expr_pca$sdev^2), 
                      pc = colnames(expr_pca$x),
                      pc_num = 1:length(expr_pca$sdev)) 

p <- ggplot(plot_dt[pc_num <= 25], aes(x = pc_num, y = pct_var_expl)) +
  geom_point() +
  geom_line() +
  labs(y = "% of variance explained", x = "PC")
p
ggsave("plots/20220809_04_DemuxHTOdCas9_FindMarkers_PerturbatuonByGeneExprPCA_PctVarExpl.pdf", p, width = 6, height = 4.5)


### (4) Run Sceptre to further validate calibration

## Data
# UMI matrices
d2n <- readRDS("processed_data/after_seurat4.3/d2n_SeuratObj_PostQC_dCas9_NewMD.RDS")
MD <- data.table(d2n@meta.data)
md <- MD[ !is.na(ident), ]

# Previous DE data
DE_dt_raw <- readRDS("processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS")
DE_dt_raw[, guide_crispr := paste0(guide_1, "-", cell_line)][]
# Remove the NFE2 outlier and the CRISPR genes
DE_dt <- merge.data.table(DE_dt_raw[!(guide_crispr == "NFE2_9-CRISPRa" | dosage_gene %in% c("NTC") | grepl("CRISPR", gene))], 
                          DE_dt_raw[gene == dosage_gene, .(avg_log2FC_dg = avg_log2FC), guide_crispr], 
                          by = "guide_crispr")
d2n_genes <- unique(DE_dt$gene)


### (1) Run Sceptre

DT <- foreach(c = c("CRISPRi", "CRISPRa"), .combine = rbind) %dopar% {
  ## Step 1: Prepare the data objects
  cells_bc <- md[ cell_line == c, L_cell_barcode]
  
  # response-by-cell matrix
  response_matrix_lowmoi <- d2n[["RNA"]]@counts[d2n_genes, cells_bc]
  
  # gRNA-by-cell matrix 
  grna_matrix_lowmoi <- d2n[["GDO"]]@counts[, cells_bc]
  
  # covariates matrix
  covariate_data_frame_lowmoi <- compute_cell_covariates(response_matrix_lowmoi)
  covariate_data_frame_lowmoi$lane_10X <- unlist(tstrsplit(colnames(response_matrix_lowmoi), "_", keep = 1)) # added 10X lane as batch
  
  # gRNA group table
  grna_group_data_frame_lowmoi <- setnames(unique(md[, .(ID, short_ID)]), c("grna_id", "grna_group"))
  grna_group_data_frame_lowmoi$grna_id <- gsub("_", "\\.", grna_group_data_frame_lowmoi$grna_id)
  grna_group_data_frame_lowmoi$grna_group[grepl("NTC", grna_group_data_frame_lowmoi$grna_group)] <- "non-targeting"
  
  # obtain the set of pairs to analyze
  response_grna_group_pairs <- generate_all_pairs(response_matrix_lowmoi,
                                                  grna_group_data_frame_lowmoi)
  
  # set the formula object
  formula_object <- formula(~log(n_umis) + 
                              log(n_nonzero) +
                              lane_10X)
  
  # run the calibration check analysis (NOTE: `calibration_check` set to `TRUE`)
  calibration_result <- run_sceptre_lowmoi(
    response_matrix = response_matrix_lowmoi,
    grna_matrix = grna_matrix_lowmoi,
    covariate_data_frame = covariate_data_frame_lowmoi,
    grna_group_data_frame = grna_group_data_frame_lowmoi,
    formula_object = formula_object,
    response_grna_group_pairs = response_grna_group_pairs,
    calibration_check = TRUE
  )
  
  # plot the calibration result to ensure adequate calibration of the p-values
  p <- plot_calibration_result(calibration_result)
  ggsave(paste0("plots/20230517_01_Sceptre_Calibration_", c, ".pdf"), p, width = 7, height = 6)
  
  # run the discovery analysis (NOTE: `calibration_check` set to `FALSE`)
  discovery_result <- run_sceptre_lowmoi(
    response_matrix = response_matrix_lowmoi,
    grna_matrix = grna_matrix_lowmoi,
    covariate_data_frame = covariate_data_frame_lowmoi,
    grna_group_data_frame = grna_group_data_frame_lowmoi,
    formula_object = formula_object,
    response_grna_group_pairs = response_grna_group_pairs,
    calibration_check = FALSE,
    n_nonzero_trt_thresh = 0
  )
  
  # compare discovery p-values to the negative control p-values; make a volcano plot
  compare_calibration_and_discovery_results(calibration_result, discovery_result)
  p <- make_volcano_plot(discovery_result)
  ggsave(paste0("plots/20230517_01_Sceptre_Volcano_", c, ".pdf"), p, width = 4, height = 4)
  
  dt <- data.table(discovery_result)
  dt[, cell_line := c]
  dt
}

DT[, guide_crispr := paste0(grna_group, "-", cell_line)]
DT[!(is.na(p_value)), p_value_fdr := p.adjust(p_value, method = "fdr")][]


### (2) Compare Sceptre FC and p-values to previous DE results
DT_merge <- merge.data.table(DE_dt, DT[, -c("cell_line"), with = FALSE], by.x = c("gene", "guide_crispr"), by.y = c("response_id", "guide_crispr"))

# Plot all FC together for all genes
cor_dt <- DT_merge[!(is.na(log_2_fold_change) | is.infinite(log_2_fold_change)), .(r = cor(log_2_fold_change, avg_log2FC)), dosage_gene]
p <- ggplot(DT_merge, aes(x = avg_log2FC, y = log_2_fold_change)) + 
  geom_point(alpha=0.2) +
  facet_grid(. ~  dosage_gene) +
  geom_abline(linetype = 2, color ="grey50") +
  geom_text(data=cor_dt, aes(x=-1.2, y=2, label=paste0("r = ", round(r,2)))) +
  labs(x = "Current log2(FC)", y= "Sceptre log2(FC)")
ggsave("plots/20230517_02_Sceptre_CorFC_ScatterAllGenes.pdf", p, width = 9, height = 3)

# Dist of correlations per trans gene
FC_cor_dt <- DT_merge[!(is.na(log_2_fold_change) | is.infinite(log_2_fold_change)), .(r = cor(log_2_fold_change, avg_log2FC)), .(dosage_gene, gene)]
p <- ggplot(FC_cor_dt, aes(x=r)) +
  geom_histogram() +
  facet_grid(. ~  dosage_gene) +
  labs(x="Cor FCs between methods", y = "Num genes")
ggsave("plots/20230517_02_Sceptre_CorFC_DistPerGene.pdf", p, width = 9, height = 2.5)

# Correlation of pvalues
cor_dt <- DT_merge[!(is.na(log_2_fold_change) | is.infinite(log_2_fold_change)), .(r = cor(p_value,pval)), dosage_gene]
p <- ggplot(DT_merge, aes(x = -log10(pval), y = -log10(p_value))) + 
  geom_point(alpha=0.2) +
  facet_grid(. ~  dosage_gene) +
  geom_abline(linetype = 2, color ="grey50") +
  geom_text(data=cor_dt, aes(x=15, y=60, label=paste0("r = ", round(r,2)))) +
  labs(x = "Current -log10(pval)", y= "Sceptre -log10(pval)")
ggsave("plots/20230517_02_Sceptre_CorPval_ScatterAllGenes.pdf", p, width = 9, height = 3)


### (3) Check how FC or num. significant DE genes are biases towards number of cells per perturbation

# Load guide property table
DE_pt <- readRDS("processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_DGGuidePropertiesDE.RDS")
DE_pt[, guide_crispr := paste0(guide_1, "-", cell_line)]

# Get cis gene FC and number of trans DE
DT_stats <- merge.data.table(DT_merge[gene == dosage_gene, .(guide_crispr, log_2_fold_change, p_value_fdr, ncells_1, dosage_gene, cell_line)],
                             DT_merge[, .(n_sigDE_fdr5 = sum(p_value_fdr <= 0.05, na.rm = T)), .(guide_crispr)],
                             by = "guide_crispr")


# Plot FC and num of DE trans genes vs number of cells
cor_dt <- DT_stats[dosage_gene != "NTC", .(r = round(cor(ncells_1, abs(log_2_fold_change)), 2), 
                                           pval = round(cor.test(ncells_1, abs(log_2_fold_change))$p.value, 3),
                                           ncells_1 = quantile(ncells_1, 0.9)), dosage_gene]
pA <- ggplot(DT_stats[dosage_gene != "NTC"], aes(y = abs(log_2_fold_change), x = ncells_1)) +
  geom_point(aes(fill=cell_line),color="white", size=4, shape=21, alpha=0.8) +
  scale_fill_manual("dosage gene", values = col_crispr, guide = "none") +
  labs(x="Number of cells with unique CRISPR perturbation", y = "Absolute cis gene log2(FC)") +
  facet_grid(. ~ dosage_gene, scales = "free_x") +
  geom_smooth(method = "lm", formula = 'y ~ x', linetype=1, color = "grey40", show.legend = F, alpha=0.2, linewidth=0.5) +
  geom_text(data=cor_dt, aes(label=paste0("r=", r, "\npval=", pval), y = 1.6)) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
pA

cor_dt <- DT_stats[dosage_gene != "NTC", .(r = round(cor(ncells_1, n_sigDE_fdr5), 2), 
                                           pval = round(cor.test(ncells_1, n_sigDE_fdr5)$p.value, 3),
                                           ncells_1 = quantile(ncells_1, 0.5)), dosage_gene]
pB <- ggplot(DT_stats[dosage_gene != "NTC"], aes(y = n_sigDE_fdr5, x = ncells_1)) +
  geom_point(aes(fill=cell_line),color="white", size=4, shape=21, alpha=0.8) +
  scale_fill_manual("dosage gene", values = col_crispr, guide = "none") +
  labs(x="Number of cells with unique CRISPR perturbation", y = "Number of DE trans genes\n(FDR < 0.05)") +
  facet_grid(. ~ dosage_gene, scales = "free_x") +
  geom_smooth(method = "lm", formula = 'y ~ x', linetype=1, color = "grey40", show.legend = F, alpha=0.2, linewidth=0.5) +
  geom_text(data=cor_dt, aes(label=paste0("r=", r, "\npval=", pval), y = -12)) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
pB
p <- plot_grid(pA, pB, align = "v", axis = "rl", nrow = 2, rel_heights = c(5/10, 5/10))

ggsave("plots/20230517_03_Sceptre_NcellsVsCisTransFC.pdf", p, width = 10.5, height = 6.5)


