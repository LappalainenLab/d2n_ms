##### D2N - CIS genes FCs - sgRNA properties #####
# (1) FC differences per sgRNA type
# (2) FC vs On-target and Off-target sgRNA activity
# (3) Attenuated sgRNAs
#     - Relationship between dosage gene FC and sgRNA mismatch position?


# JDE, July 2022
# Last modified: April 2023



## Libs
library(data.table)
library(ggplot2); theme_set(theme_bw())
library(cowplot)


## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/05-cis_fc/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


## Color palette
dosage_genes = c("GFI1B", "MYB", "NFE2", "TET2")
col_genes = RColorBrewer::brewer.pal(5, "Set1")
names(col_genes) <- c("MYB", "TET2", "NFE2", "GFI1B", "NTC")

col_crispr = list(CRISPRa = "#FF7F00", CRISPRi = "#377EB8")



## Data
# !! temp until re-running all scripts from scratch
DE_dt <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS")
DG_dt <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_DGGuidePropertiesDE.RDS")
DG_stats <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_DGStatsDE.RDS")

DG_dt[, final_guide_class := guide_class]
DG_dt[guide_class == "titration", final_guide_class := "Tiling"]
DG_dt[guide_class == "distal CRE", final_guide_class := "Enhancer"]
DG_dt[guide_class == "attenuated", final_guide_class := "Attenuated"]
DG_dt[, final_guide_class := factor(final_guide_class, levels = c("TSS", "Tiling", "Enhancer", "Attenuated", "NTC"))]

## (1) FC distributions and differences per sgRNA type

plot_dt <- copy(DG_dt)
plot_dt[guide_class == "NTC", cis_gene := "NTC"]
plot_dt[, cis_gene := factor(cis_gene, levels = c(dosage_genes, "NTC"))]
pA <- ggplot(plot_dt, aes(y = dosage_gene_log2FC)) +
  geom_hline(yintercept = 0, linetype=2) +
  geom_histogram(bins = 50, fill= "grey95", color="grey25", linewidth = 0.25) +
  scale_x_reverse() +
  scale_y_continuous(position = "right") +
  facet_grid(. ~ cis_gene, scales = "free_x") +
  theme(axis.title.y.right = element_blank()) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white"))+
  theme(panel.grid.minor = element_blank())


pB <- ggplot(DG_dt, aes(x=cis_gene, dosage_gene_log2FC)) +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(aes(alpha=sig_fdr10, fill=cell_line), size=2.5, color="white", shape=21, position = position_jitter(w = 0.2, h = 0)) +
  facet_grid(. ~ final_guide_class) +
  scale_fill_manual("cell line", values = col_crispr) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_alpha_manual("FDR < 0.1?", values = c(0.3, 0.8)) + 
  guides(alpha = guide_legend(override.aes = list(size=3, color="black"))) +
  theme(axis.title.x = element_blank()) + ylab("dosage gene log2FC") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 

p <- plot_grid(pA, pB, align = "h", axis = "tb", nrow = 1, rel_widths = c(3/10, 7/10))

ggsave(file.path(plots_dir, "05a_01_CisGenes_DistFC_FCbyGuideType.pdf"), p, width = 8, height = 3)


# (3) Attenuated sgRNAs

DG_att <- merge.data.table(DG_stats, DG_dt, 
                           by.x = c("guide_1", "cell_line", "dosage_gene_log2FC", "dosage_gene"), 
                           by.y = c("guide_1", "cell_line", "dosage_gene_log2FC", "cis_gene"))[guide_class == "attenuated", ]

cor_dt <- DG_att[, .(r = cor(MM_POS, dosage_gene_log2FC), 
                     pval = format(cor.test(MM_POS, dosage_gene_log2FC)$p.value, scientific = T, digits = 2), 
                     dosage_gene_log2FC = ifelse(cell_line == "CRISPRi", quantile(dosage_gene_log2FC, 0.11), quantile(dosage_gene_log2FC, 0.98)),
                     MM_POS = 15), 
                 .(cell_line)]
p <- ggplot(DG_att, aes(x = MM_POS, y = dosage_gene_log2FC)) +
  geom_hline(yintercept = 0, linetype=2) +
  facet_grid(cell_line ~ ., scales = "free_y") +
  geom_smooth(formula = 'y ~ x', method = "lm", color = "grey40", linewidth = 0.75) +
  geom_point(aes(color = dosage_gene), size = 2) +
  scale_color_manual("Cis gene", values = col_genes) +
  labs(x = "sgRNA missmatch position", y = "Cis gene log2(FC)") +
  geom_text(data = cor_dt, aes(label = paste0("r = ", round(r, 2), "\np = ", pval))) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
p

