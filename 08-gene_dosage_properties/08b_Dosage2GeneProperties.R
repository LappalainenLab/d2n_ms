#### D2N - Dosage to gene properties - Qualitative and quantitative gene annotations ####
# (1) Get qualitative gene annotations
# (2) Get quantitative gene annotations
#     - Scale 0-1 quantitative features
# (3) Plot mega-heatmap with annotations
# (4) Systematic correlations between quantitative measures and sigmoid parameters


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
# concat_dt <- concat_raw_dt[dosage_gene != "TET2",]

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
ggsave(file.path(plots_dir, "08a_HeatmapFC_GeneProperties.pdf"), p, width = 5,  height = 12)




## Qualitative values

# Biotype
symbols <- readRDS("../016-gene_annotations/processed_data/GA.0_IDsAnnotations.RDS")
GA.0 <- setnames(unique(symbols[, .(external_gene_name, gene_biotype)]), c("gene", "gene_biotype"))
GA.0[, value := ifelse(gene_biotype == "lncRNA", T, F)]
GA.0[, gene_set := "lncRNA"]
GA.0[, dataset := "Biotype"]

# # Mohammadi
# GA.2 <- readRDS("../016-gene_annotations/processed_data/GA.2_Mohammadi2019.RDS")

# TFs
GA.7 <- readRDS("../016-gene_annotations/processed_data/GA.7_TFs.RDS")[data_set == "Garcia-Alonso_2019",]
GA.7[, dataset := "TFs"]

# # TF targets of DGs
# GA.8 <- readRDS("../016-gene_annotations/processed_data/GA.8_TFtargets.RDS")

# OMIM
GA.9 <- readRDS("../016-gene_annotations/processed_data/GA.9_DiseaseGenesOMIM.RDS")

# GWAS
GA.10 <- readRDS("../016-gene_annotations/processed_data/GA.10_GWASgenesUKB.RDS")

# Housekeeping
GA.11 <- readRDS("../016-gene_annotations/processed_data/GA.11_HousekeepingHsiao.RDS")


cols = c("gene", "gene_set", "value", "dataset")
ga_l <- list(GA.0, GA.7, GA.9, GA.10, GA.11)

GA <- foreach(G = ga_l, .combine = rbind) %do% {
  G[, .SD, .SDcols = cols]
}



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
quantfile_v <- list.files("../016-gene_annotations/processed_data/")[grepl("^GQ.[1-6]", list.files("../016-gene_annotations/processed_data/"))]
GQ <- do.call("rbind",
              lapply(quantfile_v, function(quantfile){
                DT <- readRDS(file.path("../016-gene_annotations/processed_data/", quantfile))
              })
)

GQ[, scaled_value := rescale(value, to = c(0, 1)), metric]

gq_mat_dt <- dcast.data.table(GQ[, .(gene, metric, scaled_value)], formula = metric ~ gene, value.var = "scaled_value")
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

p <- cowplot::plot_grid(pB, pA, pC, pD, align = "h", axis = "tb", nrow = 1, rel_widths = c(1/20, 8.5/20, 1.5/15, 9/20))
p
ggsave("plots/20230404_01_Megaheatmap_QualQuantAnnot.pdf", p, width = 15, height = 12)
ggsave("plots/20230404_01_Megaheatmap_QualQuantAnnot.png", p, width = 15, height = 12)


## (4) Correlation between quantitative gene annotations and all sigmoid features

sig_feat <-  c("slope_IF", "abs_slope_IF","min_assmp","max_assmp", "min_max_range", "x_IF")

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

Cor_dt[, fdr := p.adjust(pval, method = "fdr")]
Cor_dt[, log10_pval := -log10(pval)]
Cor_dt[, log10_fdr := -log10(fdr)]


p <- ggplot(Cor_dt, aes(x= sigmoid_param, y = metric)) +
  geom_point(aes(color=r, size=log10_pval)) +
  facet_grid(dataset ~ dosage_gene, scales = "free_y", space = "free_y") +
  scale_color_gradient2(low = "#377EB8", high = "#FF7F00", midpoint = 0) +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(strip.text.y.right = element_text(angle = 0))
p
ggsave("plots/20230405_01_CorSigmoidParamsVsQuantGeneAnnot.pdf", p, width = 9, height = 9)

## Aggregated values
PT_singles[non_monotonic == F & responsive == T, .N, dosage_gene]
PT_singles_agg <- PT_singles[non_monotonic == F & responsive == T & dosage_gene != "TET2", .(slope_IF = mean(slope_IF), 
                                                                                             abs_slope_IF = mean(abs_slope_IF),
                                                                                             min_assmp = mean(min_assmp),
                                                                                             max_assmp = mean(max_assmp),
                                                                                             min_max_range = mean(min_max_range),
                                                                                             x_IF = mean(x_IF)),
                             .(gene)]

dt1 <- merge.data.table(GQ, PT_singles_agg, by = "gene")
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

Cor_dt_agg <- merge.data.table(melt(dt2, id.vars = c("metric", "dataset"), variable.name = "sigmoid_param", value.name = "pval"),
                               melt(dt3, id.vars = c("metric", "dataset"), variable.name = "sigmoid_param", value.name = "r"),
                               by = c("metric", "dataset", "sigmoid_param"))
Cor_dt_agg[, fdr := p.adjust(pval, method = "fdr")]
Cor_dt_agg[, log10_pval := -log10(pval)]

p <- ggplot(Cor_dt_agg, aes(x= sigmoid_param, y = metric)) +
  geom_point(aes(color=r, size=log10_pval)) +
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradient2(low = "#377EB8", high = "#FF7F00", midpoint = 0) +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(strip.text.y.right = element_text(angle = 0))
p
ggsave("plots/20230405s_01_CorSigmoidParamsVsQuantGeneAnnot_AggregCisGene.pdf", p, width = 7, height = 9)




