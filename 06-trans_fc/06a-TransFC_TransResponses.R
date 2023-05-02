##### D2N - TRANS genes FCs #####
# (1) Number of differential expressed trans gens with cis gene dosage changes
# (2) Trans FCs heat maps (and GFI1B one with example dosage responses highlighted)


# JDE, July 2022
# Last modified: April 2023


## Libs
library(data.table)
library(ggplot2); theme_set(theme_bw());
library(foreach)
library(doParallel); registerDoParallel(detectCores()/2)
library(cowplot)
library(vegan)


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
rect_trans_GFI1B <- c("KLK1", "SLC22A4", "RHD", "GATA2", "PITX1", "LHX3", "TUBB1", "FOXP1", "MYB")


# FC data
DE_dt_raw <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS") 
DE_dt_raw[, guide_crispr := paste0(guide_1, "-", cell_line)][]
# Remove the NFE2 outlier and the CRISPR genes
DE_dt <- merge.data.table(DE_dt_raw[!(guide_crispr == "NFE2_9-CRISPRa" | dosage_gene %in% c("NTC") | grepl("CRISPR", gene))], 
                          DE_dt_raw[gene == dosage_gene, .(avg_log2FC_dg = avg_log2FC), guide_crispr], 
                          by = "guide_crispr")
DE_dt[, sig_fdr10 := pval_fdr <= 0.1]

DG_stats <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_DGStatsDE.RDS")
DG_stats[, guide_crispr := paste0(guide_1, "-", cell_line)]


DG_dt <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_DGGuidePropertiesDE.RDS")
DG_dt[, final_guide_class := guide_class]
DG_dt[guide_class == "titration", final_guide_class := "Tiling"]
DG_dt[guide_class == "distal CRE", final_guide_class := "Enhancer"]
DG_dt[guide_class == "attenuated", final_guide_class := "Attenuated"]
DG_dt[, final_guide_class := factor(final_guide_class, levels = c("TSS", "Tiling", "Enhancer", "Attenuated", "NTC"))]
DG_dt[, guide_crispr := paste0(guide_1, "-", cell_line)]



### (1) Number of differential expressed trans gens with cis gene dosage changes
plot_dt <- merge(DG_stats, unique(DG_dt[, .(guide_crispr, final_guide_class)]), by = "guide_crispr")

p <- ggplot(plot_dt[dosage_gene != "NTC",], aes(x = dosage_gene_log2FC, y = n_sigDE_fdr5)) +
  geom_smooth(aes(color=dosage_gene), formula =  'y ~ x', method = "loess", alpha=0.2) +
  geom_point(aes(fill=dosage_gene, shape=final_guide_class),color="white", size=2, alpha=0.8) +
  geom_vline(xintercept = 0, linetype=2) +
  facet_grid(dosage_gene ~ ., scales = "free_y") +
  scale_color_manual("dosage gene", values = col_genes, guide = "none") +
  scale_fill_manual("dosage gene", values = col_genes, guide = "none") +
  scale_shape_manual("sgRNA class", values = c(22,21,23,24)) +
  guides(shape = guide_legend(override.aes = list(size=3, color="black"))) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x="Cis gene log2(FC)", y = "Num. of DE trans genes (FDR < 0.05)") +
  theme(legend.direction = "vertical") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
p
ggsave(file.path(plots_dir, "06a_01_CisFCvsTransNumDE.pdf"), p, width = 2.25, height = 7)


p <- ggplot(plot_dt[dosage_gene != "NTC",], aes(x = dosage_gene_log2FC, y = mean_abs_log2FC)) +
  geom_smooth(aes(color=dosage_gene), formula =  'y ~ x', method = "loess", alpha=0.2) +
  geom_point(aes(fill=dosage_gene, shape=final_guide_class),color="white", size=2, alpha=0.8) +
  geom_vline(xintercept = 0, linetype=2) +
  facet_grid(dosage_gene ~ .) +
  scale_color_manual("dosage gene", values = col_genes, guide = "none") +
  scale_fill_manual("dosage gene", values = col_genes, guide = "none") +
  scale_shape_manual("sgRNA class", values = c(22,21,23,24)) +
  guides(shape = guide_legend(override.aes = list(size=3, color="black"))) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x="Cis gene log2(FC)", y = "Mean abs(log2(FC)) of trans genes") +
  theme(legend.direction = "vertical") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
p
ggsave(file.path(plots_dir, "06a_02_CisFCvsTransMeanAbsFC.pdf"), p, width = 2.25, height = 7)



### (2) Trans FCs heat maps


foreach(dg = dosage_genes) %do% {
  
  G_DG <- DE_dt[ dosage_gene == dg & gene == dosage_gene, ]
  order_dosage <- G_DG[order(dosage_gene_log2FC), guide_crispr]
  G_DG[, guide_crispr := factor(guide_crispr, levels = order_dosage)]
  
  pA <- ggplot(G_DG, aes(x=guide_crispr, y=dosage_gene_log2FC)) +
    geom_bar(stat = "identity", aes(fill = cell_line, alpha=sig_fdr10)) +
    theme_classic() +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual("cell line", values = col_crispr) +
    scale_alpha_manual("FDR < 0.1?", values = c(0.3, 0.8)) +
    theme(legend.position = "") +
    ylab(paste0(dg, " log2(FC)")) +
    scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)))
  
  DG_DE <- DE_dt[dosage_gene == dg, ]
  
  # Cluster trans genes based on similar response to DG effect
  # generate a matrix of genes by dosage effects and cluster genes to get the ordered list of genes
  gxd_df <- dcast(DG_DE[, .(gene, guide_crispr, avg_log2FC)], gene ~ guide_crispr, value.var = "avg_log2FC")
  mat_g <- as.matrix(gxd_df[, 2:ncol(gxd_df)], labels = TRUE)
  rownames(mat_g) <- gxd_df$gene
  
  # Use eucledian distance as metric of similarity for hyerarchical clustering
  dist.mat <- dist(mat_g, method = "euclidean")
  clust_obj <- hclust(dist.mat)
  
  # Use log2(FC) as the weights to order the clustering of the dosage perturbations
  # Using the spearman correlation of FC agains dosage gene FC
  TRANS_dt <- setnames(DG_DE[, cor(avg_log2FC, dosage_gene_log2FC, method = "spearman"), gene], old = "V1", new ="r_dg_FC")[]
  wgts_g <- TRANS_dt[match(clust_obj$labels, TRANS_dt$gene), r_dg_FC]
  
  # Reorder the dosage perturbation cluster object based on FC as longer as possible
  clust_obj_reordered <- reorder(clust_obj, wts = wgts_g)
  
  gene_ord <- gxd_df$gene[clust_obj_reordered$order]
  DG_DE[, gene := factor(gene, levels = gene_ord)]
  
  # Set max min limits in heatmap to see smaller variation
  plot_dt <- copy(DG_DE)
  transFC_quant <- quantile(DG_DE$avg_log2FC, c(0.005, 0.995))
  plot_dt[avg_log2FC <= transFC_quant[1], avg_log2FC :=  transFC_quant[1]]
  plot_dt[avg_log2FC >= transFC_quant[2], avg_log2FC :=  transFC_quant[2]]
  plot_dt[, guide_crispr := factor(guide_crispr, levels = order_dosage)]
  
  pB <- ggplot(plot_dt, aes(x = guide_crispr, y = gene)) +
    geom_raster(aes(fill = avg_log2FC)) +
    scale_fill_gradient2("Log2(FC)", low = "#377EB8", high = "#FF7F00", midpoint = 0, na.value = "grey80") +
    theme_classic() +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
    theme(axis.line.y = element_blank(), axis.title.y = element_blank()) +
    scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
    theme(legend.position = "bottom") +
    guides(fill = guide_colourbar(barheight = 0.5,  barwidth = 6)) 
  
  if (dg == "GFI1B") {
    
    p <- plot_grid(pA, pB, align = "v", nrow = 2, rel_heights = c(1/7, 6/7))
    ggsave(file.path(plots_dir, paste0("06a_03_TransFC_HeatMap_", dg, "_NoRect.pdf")), p, width = 6, height = 12)
    
    rec_dt <- plot_dt[gene %in% rect_trans_GFI1B, .(xmin = 1, 
                                                     xmax = length(unique(plot_dt$guide_crispr)),
                                                     ymin = which(gene == gene_ord),
                                                     ymax = which(gene == gene_ord)), 
                      gene]
    
    pB2 <- pB +
      geom_rect(inherit.aes = F, data=rec_dt,
                aes(xmin = xmin-0.5,xmax=xmax+0.5,ymin=ymin-0.5,ymax = ymax+0.5),
                fill="white", alpha=0, show.legend = F, color="grey10") 
    
    p <- plot_grid(pA, pB2, align = "v", nrow = 2, rel_heights = c(1/7, 6/7))
    ggsave(file.path(plots_dir, paste0("06a_03_TransFC_HeatMap_", dg, "_Rect.pdf")), p, width = 6, height = 12)
  }
  
  p <- plot_grid(pA, pB, align = "v", nrow = 2, rel_heights = c(1/7, 6/7))
  ggsave(file.path(plots_dir, paste0("06a_03_TransFC_HeatMap_", dg, ".pdf")), p, width = 5, height = 12)
}
