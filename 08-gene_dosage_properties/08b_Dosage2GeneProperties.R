#### D2N - Dosage to gene properties - Qualitative and quantitative gene annotations ####
# (1) Get qualitative gene annotations
# (2) Get quantitative gene annotations
#     - Scale 0-1 quantitative features
# (3) Plot mega-heatmap with annotations
# (4) Systematic correlations between quantitative measures and sigmoid parameters
# (5) Systematic differences between categorical features on sigmoid parameters
# (6) Simplified and non-redundant gene annotations with dosage heatmap


# JDE, April 2023
# Last modified: August 2023



## Libs
library(ggplot2); theme_set(theme_bw())
library(data.table)
library(foreach)
library(doParallel); registerDoParallel(detectCores()/2)
library(vegan)
library(scales)



## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/08-gene_dosage_properties/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


## Load data
concat_raw_dt <- fread("../../processed_data/07-nonlinear_fit/D2N_S4LoessModelFits.txt")
range_fc = concat_raw_dt[gene == "ABRAXAS1" & dosage_gene == "GFI1B", dosage_gene_log2FC]
concat_dt <- rbind(concat_raw_dt[dosage_gene != "TET2", ], 
                   data.table(gene = c(rep("GFI1B", length(range_fc)), rep("NFE2", length(range_fc)), rep("MYB", length(range_fc))),
                              dosage_gene = c(rep("GFI1B", length(range_fc)), rep("NFE2", length(range_fc)), rep("MYB", length(range_fc))), 
                              dosage_gene_log2FC = rep(range_fc, 3),
                              avg_log2FC = rep(range_fc, 3)))


# Fold changes of all genes
DE_dt <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS")


# Sigmoid params after CV
PT_singles <- fread("../../processed_data/07-nonlinear_fit/D2N_S4Param_10fCV.txt")
PT_singles[, abs_slope_IF := abs(slope_IF)]
PT_singles[, responsive := ifelse(unresponsive == T, F, T)]
PT_singles[, cis_trans := paste0(dosage_gene, "_", gene)]




## Concatenated Log2FC heatmap across dosage genes

# Generate matrices with trans genes 
concat_mat_dt <- dcast.data.table(concat_dt[gene != dosage_gene, .(gene, dosage_gene_log2FC, avg_log2FC)], formula = gene ~ dosage_gene_log2FC, value.var = "avg_log2FC")

mat_dt <- as.matrix(concat_mat_dt[, 2:ncol(concat_mat_dt)])
rownames(mat_dt) <- concat_mat_dt$gene

# Calculate distances between rows
dist.mat_tg <- dist(mat_dt, method = "euclidean")
clust_obj_tg <- hclust(dist.mat_tg)

# Reorder gene clustering using the spearman correlation of trans genes FC agains dosage GFI1B FC
GFI1B_rdg_dt <- setnames(DE_dt[dosage_gene == "GFI1B", cor(avg_log2FC, dosage_gene_log2FC, method = "spearman"), gene], old = "V1", new ="r_dg_FC")[]
wgts_g <- GFI1B_rdg_dt[match(clust_obj_tg$labels, GFI1B_rdg_dt$gene), r_dg_FC]

# Reorder the genes based on the GFI1B trans genes FC as long as possible
clust_obj_tg_reordered <- reorder(clust_obj_tg, wts = wgts_g)
ord_tg <- clust_obj_tg_reordered$labels[clust_obj_tg_reordered$order]
concat_dt[, gene := factor(gene, levels = ord_tg)]

# Limit the color range to a threshold so that the figure is more visual
plot_dt <- copy(concat_dt)
setnames(plot_dt, old = "avg_log2FC", new = "transFC")
transFC_quant <- quantile(concat_dt$transFC, c(0.005, 0.995))
plot_dt[transFC <= transFC_quant[1], transFC :=  transFC_quant[1]]
plot_dt[transFC >= transFC_quant[2], transFC :=  transFC_quant[2]]


pA <- ggplot(plot_dt[gene != dosage_gene,], aes(x = dosage_gene_log2FC, y = gene)) +
  geom_raster(aes(fill = transFC)) +
  scale_fill_gradient2("Pred log2(FC)", low = "#377EB8", high = "#FF7F00", midpoint = 0, na.value = "grey80") +
  theme_classic() +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.line.y = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  facet_grid(. ~ dosage_gene, scales = "free_x") +
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barheight = 0.5,  barwidth = 6)) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white"))
pA

# Add dendogram to the plot
dd_d <- as.dendrogram(clust_obj_tg_reordered, type = "rectangle")
dd_d_data <- ggdendro::dendro_data(dd_d)

pB <- ggplot(ggdendro::segment(dd_d_data)) + 
  geom_segment(aes(x = y, y = x, xend = yend, yend = xend)) +
  scale_x_reverse(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  theme_classic() +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title= element_blank()) +
  ylab("") 
p <- cowplot::plot_grid(pB, pA, align = "h", axis = "tb", nrow = 1, rel_widths = c(2/10, 8/10))
# ggsave(file.path(plots_dir, "08a_HeatmapFC_GeneProperties.pdf"), p, width = 5,  height = 12)

pA2 <- ggplot(plot_dt[gene == dosage_gene,], aes(x = dosage_gene_log2FC, y = gene)) +
  geom_raster(aes(fill = transFC)) +
  scale_fill_gradient2("Pred log2(FC)", low = "#377EB8", high = "#FF7F00", midpoint = 0, na.value = "grey80") +
  theme_classic() +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.line.y = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  facet_grid(. ~ dosage_gene, scales = "free_x") +
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barheight = 0.5,  barwidth = 6)) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white"))
pA2
ggsave(file.path(plots_dir, "08a_DosageGeneFCrange.pdf"), pA2, width = 4,  height = 1.2)

## Qualitative values

qualitative_v <- list.files(processed_data_dir)[grepl("^GA.", list.files(processed_data_dir))]
GA <- do.call("rbind",
              lapply(qualitative_v, function(quantfile){
                DT <- readRDS(file.path(processed_data_dir, quantfile))
              })
)



# Order categorical annotations based on similarity 
# Generate matrix with 1/0 for each
GA[, num_value := ifelse(value == T, 1, 0)]

plot_dt <- GA
ga_mat_dt <- dcast.data.table(GA[, .(gene, gene_set, num_value)], formula = gene_set ~ gene, value.var = "num_value")

ga_mat <- as.matrix(ga_mat_dt[, 2: ncol(ga_mat_dt)])
rownames(ga_mat) <- ga_mat_dt$gene_set

ga_dist <- dist(ga_mat, method = "euclidean")
ga_clust_obj <- hclust(ga_dist)
ord_ga <- ga_clust_obj$labels[ga_clust_obj$order]

plot_dt[, gene_set := factor(gene_set, levels = ord_ga)]
plot_dt[, gene := factor(gene, levels = ord_tg)]

pC <- ggplot(plot_dt, aes(x = gene_set, y = gene)) +
  geom_tile(aes(fill=value), colour="white") +
  scale_fill_manual("", values = list(`TRUE` = "black", `FALSE` = "grey95")) +
  facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
  theme_classic() +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  theme(legend.position = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank()) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(strip.text.x.top = element_text(angle = 90))
pC



## Quantitative annotation values
quantfile_v <- list.files(processed_data_dir)[grepl("^GQ.", list.files(processed_data_dir))]
GQ <- do.call("rbind",
              lapply(quantfile_v, function(quantfile){
                DT <- readRDS(file.path(processed_data_dir, quantfile))
              })
)

GQ[, scaled_value := rescale(value, to = c(0, 1)), metric]




gq_mat_dt <- dcast(GQ[, .(gene, metric, scaled_value)], metric ~ gene, value.var = "scaled_value")
gq_mat <- as.matrix(gq_mat_dt[, 2: ncol(gq_mat_dt)])
rownames(gq_mat) <- gq_mat_dt$metric

gq_dist <- dist(gq_mat, method = "euclidean")
gq_clust_obj <- hclust(gq_dist)
ord_gq <- gq_clust_obj$labels[gq_clust_obj$order]

GQ[, metric := factor(metric, levels = ord_gq)]
GQ[, gene := factor(gene, levels = ord_tg)]
GQ[, dataset := factor(dataset, levels = c("gnomAD", "Collins et al. 2022", "Wang & Goldstein 2020", "Minaeva in prep. 2023", "String DB","Hay et al. 2018"))]

pD <- ggplot(GQ, aes(x = metric, y = gene)) +
  geom_tile(aes(fill=scaled_value), colour="white") +
  scale_fill_gradient("scaled\nvalue", low = "white", high = "#377EB8", na.value = "grey80") +
  facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
  theme_classic() +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank()) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(strip.text.x.top = element_text(angle = 90))
pD

p <- cowplot::plot_grid(pB, pA, pC, pD, align = "h", axis = "tb", nrow = 1, rel_widths = c(1/20, 7.75/20, 1.75/15, 10/20))
p
ggsave(file.path(plots_dir, "08b_PredDosage_GeneProperties_Heatmap.pdf"), p, width = 16, height = 15)

sel_gq = c("pLI", "mis_z", "pHaplo", "EDS", "n_peaks_GFI1B", "n_peaks_MYB", "n_peaks_NFE2", "Num_PPIs_whole_proteome", "Platelet", "Erythroblast", "Dendritic Cell", "Monocyte")

pD <- ggplot(GQ[metric %in% sel_gq, ], aes(x = metric, y = gene)) +
  geom_tile(aes(fill=scaled_value), colour="white") +
  scale_fill_gradient("scaled\nvalue", low = "white", high = "#377EB8", na.value = "grey80") +
  facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
  theme_classic() +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank()) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(strip.text.x.top = element_text(angle = 90))
pD

p <- cowplot::plot_grid(pB, pA, pC, pD, align = "h", axis = "tb", nrow = 1, rel_widths = c(1/15, 7/15, 2.25/15, 4.75/15))
p
ggsave(file.path(plots_dir, "08b_PredDosage_GeneProperties_Heatmap.pdf"), p, width = 16, height = 15)
ggsave("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/temp/08b_PredDosage_GeneProperties_Heatmap.pdf", p, width = 9, height = 6.5)




## (4) Correlation between quantitative gene annotations and all sigmoid features

sig_feat <-  c("slope_IF", "abs_slope_IF","min_assmp","max_assmp", "min_max_range")

Cor_dt <- foreach(dg = c("GFI1B", "MYB", "NFE2"), .combine = rbind) %do% {
  dt1 <- merge.data.table(GQ, PT_singles[dosage_gene == dg & non_monotonic == F & responsive == T, ], by = "gene")
  dt2 <- dt1[!is.na(value), {
    lapply(.SD[, sig_feat, with = FALSE], function(x) {
      cor.test(value, x)$p.value
    })
  }, by = .(metric, dataset)]
  dt3 <- dt1[!is.na(value), {
    lapply(.SD[, sig_feat, with = FALSE], function(x) {
      cor.test(value, x)$estimate
    })
  }, by = .(metric, dataset)]
  
  dt4 <- merge.data.table(melt(dt2, id.vars = c("metric", "dataset"), variable.name = "sigmoid_param", value.name = "pval"),
                          melt(dt3, id.vars = c("metric", "dataset"), variable.name = "sigmoid_param", value.name = "r"),
                          by = c("metric", "dataset", "sigmoid_param"))
  dt4[, dosage_gene := dg]
}

Cor_dt[, log10_pval := -log10(pval)]

# Order metrics same as heatmap
Cor_dt[, metric := factor(metric, levels = ord_gq)]

p <- ggplot(Cor_dt, aes(y = sigmoid_param, x = metric)) +
  geom_point(aes(color=r, size=log10_pval)) +
  facet_grid(dosage_gene ~ dataset, scales = "free_x", space = "free_x") +
  scale_color_gradient2(low = "#377EB8", high = "#FF7F00", midpoint = 0) +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.title.x = element_blank())  +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(strip.text.y.right = element_text(angle = 0)) +
  labs(y = "Sigmoid parameter") +
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barheight = 0.5,  barwidth = 6))
p
ggsave(file.path(plots_dir, "08c_SigmoidParams_GeneProperties_Cor.pdf"), p, width = 9, height = 6.5)


p <- ggplot(Cor_dt[metric %in% sel_gq, ], aes(y = sigmoid_param, x = metric)) +
  geom_point(aes(color=r, size=log10_pval)) +
  facet_grid(dosage_gene ~ dataset, scales = "free_x", space = "free_x") +
  scale_color_gradient2(low = "#377EB8", high = "#FF7F00", midpoint = 0) +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.title.x = element_blank())  +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(strip.text.y.right = element_text(angle = 0)) +
  labs(y = "Sigmoid parameter") +
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barheight = 0.5,  barwidth = 6))
p
ggsave(file.path(plots_dir, "08c_SigmoidParams_GeneProperties_Cor.pdf"), p, width = 5, height = 6.5)
ggsave("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/temp/08c_SigmoidParams_GeneProperties_Cor.pdf", p, width = 4, height = 6.5)



### (5) Systematic differences between categorical features on sigmoid parameters

Diff_dt <- foreach(dg = c("GFI1B", "MYB", "NFE2"), .combine = rbind) %do% {
  dt1 <- merge.data.table(GA, PT_singles[dosage_gene == dg & non_monotonic == F, ], by = "gene")
  dt1_melt <- melt.data.table(dt1, id.vars = colnames(dt1)[!(colnames(dt1) %in% sig_feat)], value.name = "param_value", variable.name = "param")
}
p <- ggplot(Diff_dt, aes(x=value, y=param_value)) + 
  facet_grid(param ~ gene_set, scales = "free_y") +
  geom_violin(scale = "width", fill="grey90") + 
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  geom_boxplot(width = 0.2, outlier.colour = NA) +
  labs(y = "Sigmoid parameter value", x = "Is title category?")
p

qual_feat = unique(Diff_dt$gene_set)
Diff_test_dt <- foreach(dg = c("GFI1B", "MYB", "NFE2"), .combine = rbind) %do% {
  dt1 <- Diff_dt[!is.na(value) & dosage_gene == dg,]
  dt2 <- foreach(gs = qual_feat, .combine = rbind) %do% {
    foreach(sigp = sig_feat, .combine = rbind) %do% {
      test <- wilcox.test(dt1[gene_set == gs & param == sigp & value == F, param_value], dt1[gene_set == gs & param == sigp & value == T, param_value])
      data.table(dosage_gene = dg,
                 gene_set = gs,
                 sig_param = sigp,
                 pval = test$p.value,
                 mean_diff = mean(dt1[gene_set == gs & param == sigp & value == F, param_value])- mean(dt1[gene_set == gs & param == sigp & value == T, param_value]),
                 dosage_gene = dg,
                 dataset = unique(dt1[gene_set == gs, dataset])
      )
    }
  }
}  

Diff_test_dt[, sign := ifelse(mean_diff >= 0, 1, -1)]
Diff_test_dt[, log10_pval := -log10(pval)]
Diff_test_dt[, gene_set := factor(gene_set, levels = ord_ga)]

plot_dt <- unique(Diff_test_dt[, .(gene_set, sig_param, dosage_gene, mean_diff, log10_pval,dataset)])
p <- ggplot(plot_dt, aes(x= gene_set, y = sig_param)) +
  geom_point(aes(color=mean_diff, size=log10_pval)) +
  facet_grid(dosage_gene ~ dataset , scales = "free_x", space = "free_x") +
  scale_color_gradient2(low = "#377EB8", high = "#FF7F00", midpoint = 0, mid = "grey90") +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.title.x = element_blank())  +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(strip.text.y.right = element_text(angle = 0)) +
  labs(y = "Sigmoid parameter") +
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barheight = 0.5,  barwidth = 6))
p
ggsave(file.path(plots_dir, "08d_SigmoidParams_GeneProperties_QualTest.pdf"), p, width = 4, height = 6.5)
ggsave("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/temp/08d_SigmoidParams_GeneProperties_QualTest.pdf", p, width = 3.4, height = 6.5)

# Test if un-responsive genes are enriched for houskeeping or TFs
Fisher_dt <- unique(Diff_dt[, .(gene, gene_set, dataset, dosage_gene, value, unresponsive)])
Fisher_test_dt <- Fisher_dt[, list(odds_ratio = fisher.test(table(value, unresponsive))$estimate,
                                p_value = fisher.test(table(value, unresponsive))$p.value),
                         by = .(gene_set, dataset)]
Fisher_test_dt[, fdr := p.adjust(p_value, "fdr")]
