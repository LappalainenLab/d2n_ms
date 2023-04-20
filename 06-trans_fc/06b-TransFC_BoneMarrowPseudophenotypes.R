##### D2N - Relationship trans genes FC to transcriptional states of bone marrow cell types #####
# (1) Z-score mean expression per donor per cell type of all trans genes (Supplementary plot)
# (2) Correlation of trans effects of dosage perturbation to the transcriptional state of each cell type
#     - Example correlation between dosage trans effects vs. mean expression in specific cell type
#     - Heat map of correlation 
#     - Example of dosage-to-phenotype curves


# JDE, July 2022
# Last modified: April 2023


## Libs
library(data.table)
library(ggplot2); theme_set(theme_bw());
library(foreach)
library(doParallel); registerDoParallel(detectCores()/2)
library(cowplot)
library(gridExtra)


## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/06-trans_fc/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


## Color palette
dosage_genes = c("GFI1B", "MYB", "NFE2", "TET2")
col_genes = RColorBrewer::brewer.pal(5, "Set1")
names(col_genes) <- c("MYB", "TET2", "NFE2", "GFI1B", "NTC")
col_crispr = list(CRISPRa = "#FF7F00", CRISPRi = "#377EB8")
control_genes = c("GAPDH", "LHX3")
crispr_genes = c("CRISPRi", "CRISPRa")



## Data

# FC data
DE_dt_raw <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS") 
DE_dt_raw[, guide_crispr := paste0(guide_1, "-", cell_line)][]
# Remove the NFE2 outlier and the CRISPR genes
DE_dt <- merge.data.table(DE_dt_raw[!(guide_crispr == "NFE2_9-CRISPRa" | dosage_gene %in% c("NTC") | grepl("CRISPR", gene))], 
                          DE_dt_raw[gene == dosage_gene, .(avg_log2FC_dg = avg_log2FC), guide_crispr], 
                          by = "guide_crispr")
d2n_genes <- unique(DE_dt$gene)

# Guides properties
DG_dt <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_DGGuidePropertiesDE.RDS")
DG_dt[, final_guide_class := guide_class]
DG_dt[guide_class == "titration", final_guide_class := "Tiling"]
DG_dt[guide_class == "distal CRE", final_guide_class := "Enhancer"]
DG_dt[guide_class == "attenuated", final_guide_class := "Attenuated"]
DG_dt[, final_guide_class := factor(final_guide_class, levels = c("TSS", "Tiling", "Enhancer", "Attenuated", "NTC"))]
DG_dt[, guide_crispr := paste0(guide_1, "-", cell_line)]


# Bone Marrow data
# Annotations to cell types
bm_annotations <- fread(file = file.path(data_dir, "Hay2018_SupplementaryTable_CellAssociations.txt"))
bm_annotations$cell_id <- gsub("\\.BM[1-9]", "", bm_annotations$Cell)
bm_annotations[, celltype := ClusterName]
bm_annotations[ClusterName == "CD34+ Eo/B/Mast", celltype := "CD34+ Eo-B-Mast"]
bm_annotations[ClusterName == "Naive T-cell", celltype := "Naive CD8 T-cell"]

# Get processed data from the bone marrow study
bm_mtx <- fread(file = file.path(data_dir, "GE-BM-HCA-Donor-Avg.txt"))[uid %in% d2n_genes,]
BM <- melt.data.table(bm_mtx, id.vars = "uid", variable.name = "donor_celltype")
BM$celltype <- gsub("[FM]-BM[1-9]__", "", BM$donor_celltype)
BM$celltype <- factor(BM$celltype, levels = unique(BM$celltype))
BM_avgdonor <- BM[,list(mean_expr_donors = mean(value)), .(uid, celltype)]

# All genes in the bone marrow dataset?
d2n_genes[!(d2n_genes %in% unique(BM$uid))] # Two genes not found in the BM data
d2n_bm_genes <- unique(BM$uid)



### (1) Heat map d2n genes mean expression in each cell type

# Plot Heat map of d2n genes mean expression in each cell type (using supplementary processed data)
# Scale expression values across cell types
BM_avgdonor[, expr_zscore := scale(mean_expr_donors), by="uid"]

# cluster genes for similar expression patterns across cell types
# Use euclidean distance as metric of similarity for hierarchical clustering
zscor_dt <- dcast(BM_avgdonor[, .(uid, celltype, expr_zscore)], value.var = "expr_zscore", formula = "uid ~ celltype")
dist.mat <- dist(data.matrix(zscor_dt[, 2:ncol(zscor_dt)]), method = "euclidean")
clust_obj <- hclust(dist.mat)
gene_ord <- zscor_dt$uid[clust_obj$order]
BM_avgdonor[, uid := factor(uid, levels = gene_ord)]

plot_dt <- copy(BM_avgdonor)
plot_dt[, uid := factor(uid, levels = gene_ord)]
p <- ggplot(plot_dt, aes(y=celltype, x=uid)) +
  geom_raster(aes(fill=expr_zscore)) +
  scale_fill_gradient2("Expression z-score", low = "#984EA3", high = "#4DAF4A") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", axis.title = element_blank()) +
  guides(fill = guide_colourbar(barheight = 0.5, barwidth = 10))
p
ggsave(file.path(plots_dir, "06b_01_BoneMarrowNormExpCellTypeMatrix.pdf"), p, width = 12, height = 5)




#### (2) Correlation of trans effects of dosage perturbation to the transcriptional state of each cell type

# Heat map of correlation between trans FCs and the z-score expression of each cell type

# Cis dosage effects 
Cis_dt <- DG_dt[ cis_gene %in% dosage_genes[dosage_genes != "TET2"] & !grepl("NTC", guide_1), ]

# Get ordering of dosage perturbations
order_dosage <- Cis_dt[order(dosage_gene_log2FC), guide_crispr]
Cis_dt[, guide_crispr := factor(guide_crispr, levels = order_dosage)]

# Plot dosage effect
pA <- ggplot(Cis_dt, aes(x=guide_crispr, y=dosage_gene_log2FC)) +
  facet_grid(. ~ cis_gene, scales = "free_x") +
  geom_bar(stat = "identity", aes(alpha=sig_fdr10)) +
  theme_classic() +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  scale_alpha_manual("FDR < 0.1?", values = c(0.3, 0.8)) +
  theme(legend.position = "") +
  ylab("Cis gene log2(FC)") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
pA

Cor_dt <- foreach (dg = dosage_genes[dosage_genes != "TET2"], .combine = "rbind" ) %dopar% {
  
  # Trans genes effects
  sel_trans_genes <- d2n_bm_genes
  DG_DE <- DE_dt[dosage_gene == dg & gene %in% sel_trans_genes, ]
  DG_DE[, guide_crispr := paste0(guide_1, "-", cell_line)]
  DG_DE[, guide_crispr := factor(guide_crispr, levels = order_dosage)]
  
  # Correlation between cell type mean expression and the trans gene's FC due to dosage perturbations
  BM_DE <- merge.data.table(BM_avgdonor, DG_DE[, .(gene, avg_log2FC, guide_crispr, dosage_gene_log2FC)], by.x="uid", by.y="gene", allow.cartesian=TRUE)
  BM_DE_cor <- BM_DE[, .(r = cor(mean_expr_donors, avg_log2FC), 
                         r_pval = cor.test(mean_expr_donors, avg_log2FC)$p.value,
                         r_sp = cor(mean_expr_donors, avg_log2FC, method = "spearman"), 
                         r_sp_pval = cor.test(mean_expr_donors, avg_log2FC, method = "spearman")$p.value,
                         zscore_r = cor(expr_zscore, avg_log2FC),
                         zscore_r_pval = cor.test(expr_zscore, avg_log2FC)$p.value),
                     .(guide_crispr, celltype, dosage_gene_log2FC)]

  BM_DE_cor[, r_fdr := p.adjust(r_pval, method = "fdr")]
  BM_DE_cor[, r_sp_fdr := p.adjust(r_sp_pval, method = "fdr")]
  BM_DE_cor[, zscore_r_fdr := p.adjust(zscore_r_pval, method = "fdr")]
  BM_DE_cor[, cis_gene := dg]
  BM_DE_cor
}

Cor_dt[, guide_crispr := factor(guide_crispr, levels = order_dosage)]

# Plot heat map of correlations
pB <- ggplot(Cor_dt, aes(x=guide_crispr, y = celltype, fill=zscore_r)) +
  facet_grid(. ~ cis_gene, scales = "free_x") +
  geom_tile() +
  scale_fill_gradient2("Cor fold changes with mean gene expression", low = "#377EB8", high = "#FF7F00", midpoint = 0) +
  theme_classic() +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank(), legend.position = "bottom") +
  geom_point(data = Cor_dt[zscore_r_fdr < 0.1,], shape=8, color="grey50", size=0.5) +
  guides(fill = guide_colourbar(barheight = 0.5, barwidth = 10)) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white"))

p <- cowplot::plot_grid(pA, pB, align = "v", nrow = 2, rel_heights = c(1/5, 3/5))
p
ggsave(file.path(plots_dir, "06b_02_CorHeatmapFCvsNormExprCelltypes.pdf"), p, width = 12, height = 6)


## Examples of dosage-to-phenotype curves

# Identify the cell types with significant correlations to trans genes FCs
cis_celltypes_v <- Cor_dt[, sum(zscore_r_fdr < 0.1) >= 5, .(celltype, cis_gene)][V1 == T, paste0(cis_gene, " - ", celltype)]


# Dose response curves vs z-score expression acros donors per cell type
min_r = min(Cor_dt$r)-0.1
max_r = max(Cor_dt$r)+0.1

plot_cis_celltypes_v <- cis_celltypes_v[cis_celltypes_v !=  "NFE2 - Erythroblast" ]

Cor_dt[, cis_celltype := paste0(cis_gene, " - ", celltype)]

plot_dt <- Cor_dt[ cis_celltype %in% plot_cis_celltypes_v, ]
plot_dt[, cis_celltype := factor(cis_celltype, levels = c("GFI1B - Erythroblast","MYB - Erythroblast","GFI1B - Pre-Dendritic","NFE2 - Platelet" ))]
p <- ggplot(plot_dt, aes(x= dosage_gene_log2FC, y=zscore_r)) +
  geom_vline(xintercept = 0, linetype=2, color="grey50") +
  geom_hline(yintercept = 0, linetype=2, color="grey50") +
  geom_point() +
  geom_smooth(method = "loess", formula = "y ~ x", color="#FF7F00", alpha=0.2) +
  facet_wrap( . ~ cis_celltype, ncol = 2, scales = "free_x") +
  labs(x="Cis gene log2(FC)", y="Pearson corr of trans genes FC to z-scored celltype mean expression") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white"))
ggsave(file.path(plots_dir, "06b_03_ExampleDosage2PhenCurves.pdf"), p, width = 5, height = 5)





