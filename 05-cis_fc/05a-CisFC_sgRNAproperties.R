##### D2N - CIS genes FCs - sgRNA properties #####
# (1) FC differences per sgRNA type
# (2) FC vs On-target and Off-target sgRNA activity
# (3) Attenuated sgRNAs
#     - Relationship between dosage gene FC and sgRNA mismatch position?
#     - Number of trans genes and sgRNA mismatch position?

# JDE, July 2022
# Last modified: February 2023



## Libs
library(data.table)
library(ggplot2); theme_set(theme_bw())


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



## Data
# !! temp until re-running all scripts from scratch
DE_dt <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS")
DG_dt <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_DGGuidePropertiesDE.RDS")
DG_stats <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_DGStatsDE.RDS")

DG_att <- merge.data.table(DG_stats, DG_dt, 
                           by.x = c("guide_1", "cell_line", "dosage_gene_log2FC", "dosage_gene"), 
                           by.y = c("guide_1", "cell_line", "dosage_gene_log2FC", "cis_gene"))[guide_class == "attenuated", ]

cor_dt <- DG_att[, .(r = cor(MM_POS, dosage_gene_log2FC), 
                     pval = format(cor.test(MM_POS, dosage_gene_log2FC)$p.value, scientific = T, digits = 2), 
                     dosage_gene_log2FC = ifelse(cell_line == "CRISPRi", quantile(dosage_gene_log2FC, 0.11), quantile(dosage_gene_log2FC, 0.98)),
                     MM_POS = 15), 
                 .(cell_line)]
ggplot(DG_att, aes(x = MM_POS, y = dosage_gene_log2FC)) +
  geom_hline(yintercept = 0, linetype=2) +
  facet_grid(cell_line ~ ., scales = "free_y") +
  geom_smooth(formula = 'y ~ x', method = "lm", color = "grey40", linewidth = 0.75) +
  geom_point(aes(color = dosage_gene), size = 2) +
  scale_color_manual("Cis gene", values = col_genes) +
  labs(x = "sgRNA missmatch position", y = "Cis gene log2(FC)") +
  geom_text(data = cor_dt, aes(label = paste0("r = ", round(r, 2), "\np = ", pval)))


ggplot(DG_att[dosage_gene != "TET2", ], aes(x = MM_POS, y = n_sigDE_fdr5)) +
  geom_hline(yintercept = 0, linetype=2) +
  facet_grid(cell_line ~ ., scales = "free_y") +
  geom_smooth(formula = 'y ~ x', method = "lm", color = "grey40", linewidth = 0.75) +
  geom_point(aes(color = dosage_gene), size = 2) +
  scale_color_manual("Cis gene", values = col_genes) +
  labs(x = "sgRNA missmatch position", y = "Num. trans effects (FDR < 0.05)") 
