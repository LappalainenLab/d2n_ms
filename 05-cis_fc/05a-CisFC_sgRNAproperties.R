##### D2N - CIS genes FCs - sgRNA properties #####
# (1) FC differences per sgRNA type
# (2) FC vs. On-target and Off-target sgRNA activity
# (3) FC vs. number of cells per perturbation
# (4) Attenuated sgRNAs - Dosage gene FC and sgRNA mismatch position from PAM


# JDE, July 2022
# Last modified: April 2023



## Libs
library(data.table)
library(ggplot2); theme_set(theme_bw());
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




### (1) FC distributions and differences per sgRNA type

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



### (2) ON-/OFF-target activity of gRNAs
plot_dt <- melt.data.table(DG_dt[!grepl("NTC", guide_1), .(dosage_gene_log2FC, Doench2014OnTarget, otCount, guide_class, guide_1, cell_line)],
                           id.vars = c("dosage_gene_log2FC", "guide_class", "guide_1", "cell_line"), 
                           variable.name = "on_off_param", 
                           value.name = "value")
plot_dt[, parameter := on_off_param]
plot_dt[, parameter := ifelse(on_off_param == "Doench2014OnTarget", "ON-target activity", "OFF-target counts")]
plot_dt_cor = plot_dt[!is.na(value), .(r = round(cor(abs(dosage_gene_log2FC), value), 2),
                                       pval = round(cor.test(abs(dosage_gene_log2FC), value)$p.value, 3)),
                      .(cell_line, parameter)]
plot_dt_cor[, ycrisp := c(0.78, 0.92, 350, 420)]
p <- ggplot(plot_dt[guide_class != "NTC", ], aes(x=abs(dosage_gene_log2FC), y=value)) +
  geom_point(aes(fill=cell_line), size=2.5, color="white", shape=21, alpha=0.75, show.legend = F) +
  facet_grid(parameter ~ ., scales="free_y") +
  scale_fill_manual("cell line", values = col_crispr) +
  labs(x= "Absolute cis gene log2(FC)", y="Property value") +
  geom_smooth(method = "lm", formula = 'y ~ x', linetype=1, aes(color=cell_line), show.legend = F, alpha=0.1, linewidth=0.5) +
  scale_color_manual("cell line", values = col_crispr) +
  geom_text(data = plot_dt_cor, aes(x=0.75, y=ycrisp, color=cell_line, label = paste0(cell_line, " r = ", r, "\npval = ", pval)), show.legend = F, size =3) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
p
ggsave(file.path(plots_dir, "05a_02_OnOffActivityVsFC.pdf"), p, width = 3.3, height = 5)



### (3) FC vs. number of cells per perturbation

cor_dt <- DG_stats[dosage_gene != "NTC", .(r = round(cor(ncells_1, abs(dosage_gene_log2FC)), 2), 
                                 pval = round(cor.test(ncells_1, abs(dosage_gene_log2FC))$p.value, 3),
                                 ncells_1 = quantile(ncells_1, 0.9)), dosage_gene]
pA <- ggplot(DG_stats[dosage_gene != "NTC"], aes(y = abs(dosage_gene_log2FC), x = ncells_1)) +
  geom_point(aes(fill=cell_line),color="white", size=4, shape=21, alpha=0.8) +
  scale_fill_manual("dosage gene", values = col_crispr, guide = "none") +
  labs(x="Number of cells with unique CRISPR perturbation", y = "Absolute cis gene log2(FC)") +
  facet_grid(. ~ dosage_gene, scales = "free_x") +
  geom_smooth(method = "lm", formula = 'y ~ x', linetype=1, color = "grey40", show.legend = F, alpha=0.2, linewidth=0.5) +
  geom_text(data=cor_dt, aes(label=paste0("r=", r, "\npval=", pval), y = 1.6)) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 


cor_dt <- DG_stats[dosage_gene != "NTC", .(r = round(cor(ncells_1, n_sigDE_fdr5), 2), 
                                           pval = round(cor.test(ncells_1, abs(dosage_gene_log2FC))$p.value, 3),
                                           ncells_1 = quantile(ncells_1, 0.5)), dosage_gene]
pB <- ggplot(DG_stats[dosage_gene != "NTC"], aes(y = n_sigDE_fdr5, x = ncells_1)) +
  geom_point(aes(fill=cell_line),color="white", size=4, shape=21, alpha=0.8) +
  scale_fill_manual("dosage gene", values = col_crispr, guide = "none") +
  labs(x="Number of cells with unique CRISPR perturbation", y = "Number of DE trans genes\n(FDR < 0.05)") +
  facet_grid(. ~ dosage_gene, scales = "free_x") +
  geom_smooth(method = "lm", formula = 'y ~ x', linetype=1, color = "grey40", show.legend = F, alpha=0.2, linewidth=0.5) +
  geom_text(data=cor_dt, aes(label=paste0("r=", r, "\npval=", pval), y = -12)) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 


p <- plot_grid(pA, pB, align = "v", axis = "rl", nrow = 2, rel_widths = c(5/10, 5/10))

ggsave(file.path(plots_dir, "05a_03_NcellsVsCisTransFC.pdf"), p, width = 10.5, height = 6.5)


### (4) Attenuated sgRNAs

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
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(legend.position = "bottom", legend.direction = "vertical")
p
ggsave(file.path(plots_dir, "05a_04_AttenuatedDistFromPamVsFC.pdf"), p, width = 3.3, height = 6.5)

