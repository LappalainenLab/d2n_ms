##### D2N - CIS genes FCs - sgRNA properties and possible biases #####
# (1) FC differences per sgRNA type
# (2) FC vs. On-target and Off-target sgRNA activity
# (3) FC vs. number of cells per perturbation
# (4) Attenuated sgRNAs - Dosage gene FC and sgRNA mismatch position from PAM
# (5) Distribution of cis genes UMIs per sgRNA and correlation of trans effects when binarizing cells into strong/mild cis gene effect

# JDE, July 2022
# Last modified: April 2023



## Libs
library(data.table)
library(ggplot2); theme_set(theme_bw());
library(cowplot)
library(foreach)
library(doParallel); registerDoParallel(detectCores()/2)
library(ggridges)
library(SeuratObject)


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
DG_stats[, guide_crispr := paste0(guide_1, "-", cell_line)]

DG_dt[, final_guide_class := guide_class]
DG_dt[guide_class == "titration", final_guide_class := "Tiling"]
DG_dt[guide_class == "distal CRE", final_guide_class := "Enhancer"]
DG_dt[guide_class == "attenuated", final_guide_class := "Attenuated"]
DG_dt[, final_guide_class := factor(final_guide_class, levels = c("TSS", "Tiling", "Enhancer", "Attenuated", "NTC"))]

#UMI matrix
d2n <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_SeuratObj_PostQC_dCas9_NewMD.RDS")
MD <- as.data.table(d2n@meta.data)


### (0) Stats and numbers for manuscript 

## Min and Max FC across all target genes
DG_dt[, .(min_log2fc = min(dosage_gene_log2FC), max_log2fc = max(dosage_gene_log2FC))]
# min_log2fc max_log2fc
# 1:   -1.82795  0.7986941

## Min and Max FC for each target gene
DG_dt[, .(min_log2fc = min(dosage_gene_log2FC), max_log2fc = max(dosage_gene_log2FC)), .(cis_gene)]
# cis_gene min_log2fc max_log2fc
# 1:    GFI1B -1.8279497  0.5068685
# 2:      MYB -0.9441177  0.1774883
# 3:     NFE2 -0.8146770  0.2810703
# 4:     TET2 -0.4309872  0.7986941

## Accuracy of target gene FC direction and cell line of origin
acc_dt <- DG_dt[sig_fdr10 == T, .(cell_line, dosage_gene_log2FC)]
acc_dt[, fc_dir_cell_line := ifelse(dosage_gene_log2FC < 0, "CRISPRi", "CRISPRa")]
acc_dt[, sum(cell_line == fc_dir_cell_line) / .N * 100]
# [1] 98.8764



### (1) FC distributions and differences per sgRNA type

plot_dt <- copy(DG_dt)
plot_dt[guide_class == "NTC", cis_gene := "NTC"]
plot_dt[, cis_gene := factor(cis_gene, levels = c(dosage_genes, "NTC"))]

pA <- ggplot(DG_dt, aes(x=cis_gene, dosage_gene_log2FC)) +
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(aes(alpha=sig_fdr10, fill=cell_line), size=2.5, color="white", shape=21, position = position_jitter(w = 0.2, h = 0)) +
  facet_grid(. ~ final_guide_class) +
  scale_fill_manual("cell line", values = col_crispr) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_alpha_manual("FDR < 0.1?", values = c(0.3, 0.8)) + 
  guides(alpha = guide_legend(override.aes = list(size=3, color="black"))) +
  theme(axis.title.x = element_blank()) + ylab("dosage gene log2FC") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white"), legend.position = "left") +
  coord_cartesian(ylim = c(min(plot_dt$dosage_gene_log2FC)- 0.1, max(plot_dt$dosage_gene_log2FC)+ 0.1))

pB <- ggplot(plot_dt, aes(y = dosage_gene_log2FC)) +
  geom_hline(yintercept = 0, linetype=2) +
  geom_histogram(bins = 50, fill= "grey95", color="grey25", linewidth = 0.25) +
  facet_grid(. ~ cis_gene, scales = "free_x") +
  theme(axis.title.y.left = element_blank()) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white"))+
  theme(panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = c(min(plot_dt$dosage_gene_log2FC)- 0.1, max(plot_dt$dosage_gene_log2FC)+ 0.1))

p <- plot_grid(pA, pB, align = "h", axis = "tb", nrow = 1, rel_widths = c(7/10, 3/10))

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

# For trans effects remove guide 9 for NFE2 since trans response is complete outlier
cor_dt <- DG_stats[!(dosage_gene == "NTC" | guide_crispr == "NFE2_9-CRISPRa"), .(r = round(cor(ncells_1, n_sigDE_fdr5), 2), 
                                           pval = round(cor.test(ncells_1, n_sigDE_fdr5)$p.value, 3),
                                           ncells_1 = quantile(ncells_1, 0.5)), dosage_gene]
pB <- ggplot(DG_stats[dosage_gene != "NTC"], aes(y = n_sigDE_fdr5, x = ncells_1)) +
  geom_point(aes(fill=cell_line),color="white", size=4, shape=21, alpha=0.8) +
  scale_fill_manual("dosage gene", values = col_crispr, guide = "none") +
  labs(x="Number of cells with unique CRISPR perturbation", y = "Number of DE trans genes\n(FDR < 0.05)") +
  facet_grid(. ~ dosage_gene, scales = "free_x") +
  geom_smooth(method = "lm", formula = 'y ~ x', linetype=1, color = "grey40", show.legend = F, alpha=0.2, linewidth=0.5) +
  geom_text(data=cor_dt, aes(label=paste0("r=", r, "\npval=", pval), y = -12)) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 

p <- plot_grid(pA, pB, align = "v", axis = "rl", nrow = 2, rel_heights = c(5/10, 5/10))

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


# (5) Distribution of cis genes UMIs per sgRNA and correlation of trans effects when binarizing cells into strong/mild cis gene effect
col_cell_line = list(CRISPRa = "#FF7F00", CRISPRi = "#377EB8", CRISPRi_NTC = "#afcbe2", CRISPRa_NTC = "#ffcb99")

# Plot the UMI distribution of each dosage gene in each unique perturbation
DE_cis <- rbind(DE_dt[gene == dosage_gene, ], DE_dt[dosage_gene == "NTC" & gene %in% dosage_genes,])
DE_cis[, guide_crispr := paste0(guide_1, "-", cell_line)]

# Iterate over dosage genes and get UMIs of the different cell categories
pL <- foreach(dg = dosage_genes) %do% {
  
  DG_dt <- DE_cis[gene == dg,]
  gRNA_ord_vec <- DG_dt[order(avg_log2FC), guide_crispr]
  
  dt <- foreach(gRNA = gRNA_ord_vec, .combine = rbind) %dopar% {
    cells_indx <- MD[guide_crispr == gRNA, L_cell_barcode]
    data.table(dosage_gene_expr = d2n[["RNA"]]@data[dg, cells_indx],
               guide_crispr = gRNA)
  }
  dt[, cell_line := tstrsplit(guide_crispr, "-", keep = 2),]
  dt[, guide_class := ifelse(grepl("NTC", guide_crispr), "NT", "T")]
  dt[, guide_crispr := factor(guide_crispr, levels = gRNA_ord_vec)]
  dt[, line_ntc := cell_line]
  dt[guide_class == "NT", line_ntc := paste0(cell_line, "_", "NTC")]
  
  ggplot(dt, aes(x=dosage_gene_expr, y=guide_crispr)) +
    geom_density_ridges(scale = 2, size = 0.5, rel_min_height = 0.001, aes(fill = line_ntc), quantile_lines = TRUE, quantiles = 2, show.legend = F) +
    labs(x=paste0(dg, " norm. expression"), y="Unique guides") +
    scale_fill_manual("", values = col_cell_line) +
    facet_wrap(. ~ cell_line, ncol = 1, scales = "free_y") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) + 
    theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white"))
}

p <- do.call("plot_grid",  c(pL, nrow=1, align = "v"))
p
ggsave(file.path(plots_dir, "05a_05_NormUMIDistCisGene.pdf"), p, width = 8, height = 6.5)



# Run differential expression analyses when splitting guide_crispr from cell that have 0 or non-zero norm expr for that dosage gene

# Create Metadata table with cells that have only dosage genes, including the UMI value of the dosage gene
MD_tg <- MD[gene != "ntc", ]
MD_tg$dosage_gene_umi <- foreach (x = 1:nrow(MD_tg), .combine = c) %dopar% {
  d2n[["RNA"]]@data[MD_tg[x, gene], MD_tg[x, L_cell_barcode]]
}

# Classify low vs high dosage gene UMI 
MD_tg[, UMI_class := ifelse(dosage_gene_umi > 0, "high", "low")]

# New MD with umi classification
MD_new <- merge.data.table(MD, MD_tg[, .(L_cell_barcode, UMI_class, dosage_gene_umi)], all.x = T, by = "L_cell_barcode")
MD_new[, guide_crispr_umi := ifelse(gene == "ntc", ident, paste0(ident, "__", UMI_class)) ]

# Add new metadata into Seurat object
d2n@meta.data$guide_crispr_umi <- MD_new[match(rownames(d2n@meta.data), L_cell_barcode), guide_crispr_umi]


# Run FindMarkers between every guide_CRISPR condition vs well callibrated NT_CRISPR conditions
Idents(d2n) <- "guide_crispr_umi"
guide_crispr_v <- unique(d2n$guide_crispr_umi)[!(grepl("NTC", unique(d2n$guide_crispr_umi)) | is.na(unique(d2n$guide_crispr_umi)))]


low_cell_ids <- unique(MD_new[, .N, guide_crispr_umi][N <= 10, tstrsplit(guide_crispr_umi, "__", keep = 1)]$V1)

guide_crispr_test_v <- guide_crispr_v[!(guide_crispr_v %in% c(paste0(low_cell_ids, "__high"), paste0(low_cell_ids, "__low")))]


DE_lh_dt <- foreach(x = 1:length(guide_crispr_test_v), .combine = rbind) %dopar% {
  
  md <- as.data.table(d2n@meta.data)
  
  id1 = guide_crispr_test_v[x]
  cells1 = md[guide_crispr_umi == id1, L_cell_barcode]
  cell_line = tstrsplit(tstrsplit(id1, "-", keep = 2)[[1]], "__", keep = 1)[[1]]
  id2 = paste0("NTC-", cell_line)
  cells2 = md[guide_crispr_umi == id2, L_cell_barcode]
  gene = tstrsplit(id1, "_", keep = 1)[[1]]
  
  FM_df <- FindMarkers(object = d2n[["RNA"]], cells.1 = cells1, cells.2 = cells2, min.pct = -Inf, logfc.threshold = -Inf, verbose = F)
  FM_dt <- data.table(gene = rownames(FM_df), 
                      avg_log2FC = FM_df$avg_log2FC,
                      pval = FM_df$p_val,
                      pval_fdr = p.adjust(FM_df$p_val, method = "fdr"),
                      idnt = id1, 
                      guide = tstrsplit(id1, "-", keep = 1)[[1]],
                      cell_line = cell_line,
                      UMI_class = tstrsplit(id1, "__", keep = 2)[[1]],
                      ncells = length(cells1),
                      dosage_gene = gene,
                      dosage_gene_log2FC = ifelse(gene == "NTC", NA, FM_df[gene, "avg_log2FC"]),
                      dosage_gene_pval_fdr = ifelse(gene == "NTC", NA, FM_df[gene, "p_val_adj"]))
}

# Data table to compare High vs Low
DE_LvH_dt <- merge.data.table(DE_lh_dt[UMI_class == "high" & !grepl("CRISPR", gene), ], DE_lh_dt[UMI_class == "low" & !grepl("CRISPR", gene), .(gene, avg_log2FC, guide, cell_line, ncells)], 
                              by = c("gene", "guide", "cell_line"), suffixes = c("_highUMI", "_lowUMI"))

# Calculate correlations between trans effects
DE_LvH_cor <- DE_LvH_dt[, .(r = cor(avg_log2FC_lowUMI, avg_log2FC_highUMI)), .(guide, cell_line, dosage_gene, ncells_highUMI, ncells_lowUMI)][]

# Include the difference in number of cells
DE_LvH_cor[, diff_ncells := abs(ncells_highUMI - ncells_lowUMI)]

# Include the cis gene fold change estimation when using all cells
DE_LvH_cor_fc <- merge.data.table(DE_LvH_cor, DE_cis[gene == dosage_gene, .(guide_1, cell_line, dosage_gene, avg_log2FC)],
                                  by.x = c("dosage_gene", "guide", "cell_line"), by.y = c("dosage_gene", "guide_1", "cell_line"))

pA <- ggplot(DE_LvH_cor, aes(x = r)) +
  geom_histogram(aes(fill=cell_line), bins = 30, show.legend = F) +
  facet_grid(. ~ dosage_gene, scales = "free_y") +
  scale_fill_manual("cell line", values = col_crispr) +
  geom_vline(xintercept = 0, linetype=2, color="grey50") +
  coord_cartesian(xlim=c(-1, 1)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(panel.grid.minor.y = element_blank())

pB <- ggplot(DE_LvH_cor_fc, aes(y = avg_log2FC, x = r)) +
  geom_vline(xintercept = 0, linetype=2, color="grey50") +
  geom_hline(yintercept = 0, linetype=2, color="grey50") +
  geom_point(aes(size=diff_ncells, fill=cell_line), color="white", shape=21) +
  facet_grid(. ~ dosage_gene) +
  theme_bw() +
  scale_fill_manual("cell line", values = col_crispr) +
  xlab("Correlation log2(FC) trans effects between cells with High (Norm. UMI > 0) vs. Low cis gene expresion (Norm. UMI = 0)") +
  ylab("Cis gene log2(FC) (All cells)") +
  guides(size = guide_legend(override.aes = list(color="black"))) +
  theme(legend.position = "bottom") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white"))+
  coord_cartesian(xlim=c(-1, 1))

p <- plot_grid(pA, pB, align = "v", axis = "rl", nrow = 2, rel_heights = c(2.5/10, 7.5/10))

ggsave(file.path(plots_dir, "05a_05_CorLowsHighNormUMIvsCisGeneFC.pdf"), p, width = 9, height = 4.5)




