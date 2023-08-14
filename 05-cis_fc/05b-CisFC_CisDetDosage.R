##### D2N - CIS genes FCs - Cis determinants of dosage #####
# (1) CRISPRi vs. CRISPRa
# (2) FC vs. distance from TSS
# (3) FC vs. presence/absence of epigenetic marks/open peaks


# JDE, July 2022
# Last modified: April 2023



## Libs
library(data.table)
library(ggplot2); theme_set(theme_bw());
library(foreach)
library(doParallel); registerDoParallel(detectCores()/2)
library(ggridges)
library(cowplot)
library(RColorBrewer)
library(rtracklayer)
library(ggupset)
library(GenomicRanges)


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

DG_dt[, final_guide_class := guide_class]
DG_dt[guide_class == "titration", final_guide_class := "Tiling"]
DG_dt[guide_class == "distal CRE", final_guide_class := "Enhancer"]
DG_dt[guide_class == "attenuated", final_guide_class := "Attenuated"]
DG_dt[, final_guide_class := factor(final_guide_class, levels = c("TSS", "Tiling", "Enhancer", "Attenuated", "NTC"))]

DG_dt[, guide_crispr := paste0(guide_1, "-", cell_line)]


### (0) Stats and numbers for manuscript

# max effect of CRISPRi and CRISPRa relative to TSS distance
x_values = seq(from = -500, to = 500, by = 1)

# Ca
l_a <- loess(DG_dt[guide_category %in% c("titration", "weismanTSS") & cell_line == "CRISPRa", .(dosage_gene_log2FC, dist2tss)], formula = "dosage_gene_log2FC ~ dist2tss" )
pred_a <- predict(l_a, newdata = data.frame(dist2tss = x_values))
x_values[which.max(pred_a)]
# [1] -99

# Ci
l_i <- loess(DG_dt[guide_category %in% c("titration", "weismanTSS") & cell_line == "CRISPRi", .(dosage_gene_log2FC, dist2tss)], formula = "dosage_gene_log2FC ~ dist2tss" )
pred_i <- predict(l_i, newdata = data.frame(dist2tss = x_values))
x_values[which.min(pred_i)]
# [1] 238

### (1) CRISPRi vs. CRISPRa cis gene FCs

DE_merged <- merge.data.table(DG_dt[cell_line == "CRISPRa" & final_guide_class != "NTC", .(guide_1, dosage_gene_log2FC, dosage_gene_pval_fdr)],
                              DG_dt[cell_line == "CRISPRi" & final_guide_class != "NTC", .(guide_1, dosage_gene_log2FC, dosage_gene_pval_fdr, cis_gene, final_guide_class)],
                              suffixes = c("_Ca", "_Ci"), by = "guide_1", all.x = T)

p <- ggplot(DE_merged, aes(x=dosage_gene_log2FC_Ci, y=dosage_gene_log2FC_Ca)) +
  geom_point(aes(fill=cis_gene, shape=final_guide_class), color="white", size=3, alpha=0.85) +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  labs(x="cis gene log2(FC) (CRISPRi)", y="cis gene log2(FC) (CRISPRa)") +
  scale_fill_manual("Cis gene", values = col_genes[dosage_genes]) +
  scale_shape_manual("Guide class", values = c(22,21,23,24)) +
  guides(shape = guide_legend(override.aes = list(size=3, color="black")),
         fill = guide_legend(override.aes = list(size=3, color=col_genes[c("GFI1B", "MYB", "NFE2", "TET2")])))
p
ggsave(file.path(plots_dir, "05b_01_CisGenes_CiVsCa.pdf"), p, width = 4.5, height = 3)


## (2) Distance from TSS
max_fc = max(DG_dt[!(guide_category %in% c("attenuated", "NTC")), dosage_gene_log2FC])+0.1
min_fc = min(DG_dt[!(guide_category %in% c("attenuated", "NTC")), dosage_gene_log2FC])-0.1



plot_dt <- DG_dt[!(guide_category %in% c("attenuated", "distalCRE", "NTC")), ]
pB <- ggplot(plot_dt, aes(x=dist2tss, y=dosage_gene_log2FC)) +
  geom_smooth(aes(color=cell_line), show.legend = F, alpha=0.1, method = "loess" ,  formula = "y ~ x") +
  scale_color_manual("cell line", values = col_crispr) +
  geom_hline(yintercept = 0, linetype=2) + geom_vline(xintercept = 0, linetype=1) +
  geom_point(aes(alpha=sig_fdr10, fill=cell_line), size=2.5, color="white", shape=21) +
  scale_fill_manual("cell line", values = col_crispr) +
  scale_alpha_manual("FDR < 0.1?", values = c(0.4, 0.8)) + 
  guides(alpha = guide_legend(override.aes = list(size=3, color="black"))) +
  labs(y = "dosage gene log2FC", x = "distance to TSS (bp)") +
  coord_cartesian(ylim = c(min_fc, max_fc))
pB

yrange = ggplot_build(pB)$layout$panel_params[[1]]$y.range
xrange = ggplot_build(pB)$layout$panel_params[[1]]$x.range



plot_dt <- DG_dt[!(guide_category %in% c("attenuated", "NTC")), ]
pA <- ggplot(plot_dt, aes(x=dist2tss, y=dosage_gene_log2FC)) +
  geom_rect(mapping=aes(ymin = yrange[1], ymax = yrange[2], xmin = xrange[1], xmax = xrange[2]), fill = NA, color="grey20", linewidth = 0.1, linetype=3) +
  geom_hline(yintercept = 0, linetype=2) + geom_vline(xintercept = 0, linetype=1) +
  geom_point(aes(alpha=sig_fdr10, fill=cell_line), size=2.5, color="white", shape=21) +
  scale_fill_manual("cell line", values = col_crispr) +
  scale_alpha_manual("FDR < 0.1?", values = c(0.4, 0.8)) + 
  guides(alpha = guide_legend(override.aes = list(size=3, color="black"))) +
  labs(y = "dosage gene log2FC", x = "") +
  coord_cartesian(ylim = c(min_fc, max_fc))
pA

p <- plot_grid(pA, pB, align = "v", axis = "rl", nrow = 2, rel_heights = c(4/10, 6/10))
ggsave(file.path(plots_dir, "05b_02_CisGenes_FCvsTSSdist.pdf"), p, width = 6, height = 5)



## (3) Epigenetic/open chromatin marks determinants of dosage

# BED narrowPeak files from ENCODE
# To import narrowPeak files
extraCols_narrowPeak <- c(signalValue = "numeric", 
                          pValue = "numeric", 
                          qValue = "numeric", 
                          peak = "integer")
H3K4me3_gr <- import("/gpfs/commons/groups/lappalainen_lab/jdomingo/data/ENCODE/K562/ChIPseq-H3K4me3_GRCh38_narrowPeak.bed", format = "BED",extraCols = extraCols_narrowPeak)
H3K4me1_gr <- import("/gpfs/commons/groups/lappalainen_lab/jdomingo/data/ENCODE/K562/ChIPseq-H3K4me1_GRCh38_narrowPeak.bed", format = "BED",extraCols = extraCols_narrowPeak)
H3K27ac_gr <- import("/gpfs/commons/groups/lappalainen_lab/jdomingo/data/ENCODE/K562/ChIPseq-H3K27ac_GRCh38_narrowPeak.bed", format = "BED",extraCols = extraCols_narrowPeak)
H3K9ac_gr <- import("/gpfs/commons/groups/lappalainen_lab/jdomingo/data/ENCODE/K562/ChIPseq-H3K9ac_GRCh38_narrowPeak.bed", format = "BED",extraCols = extraCols_narrowPeak)
ATAC_gr <- import("/gpfs/commons/groups/lappalainen_lab/jdomingo/data/ENCODE/K562/ATACseq_GRCh38_narrowPeak.bed", format = "BED",extraCols = extraCols_narrowPeak)
DNase_gr <- import("/gpfs/commons/groups/lappalainen_lab/jdomingo/data/ENCODE/K562/DNaseseq_GRCh38_narrowPeak.bed", format = "BED",extraCols = extraCols_narrowPeak)

EPI_l <- list("H3K4me3" = H3K4me3_gr, "H3K4me1" = H3K4me1_gr, "H3K27ac" =  H3K27ac_gr, "H3K9ac" = H3K9ac_gr, "ATAC" = ATAC_gr, "DNase" = DNase_gr)


# Overlap gRNA target site with epigenetic features

# Remove NTC and attenuated guides 
GDO_df <- rbind.data.frame(
  setnames(DG_dt[guide_category %in% c("titration", "weismanTSS", "distalCRE") & GUIDE_STRAND == "+", .(chr, coord_start, coord_end, GUIDE_STRAND, guide_1, guide_class)],
           old = c("GUIDE_STRAND"), new=c("strand")),
  setnames(DG_dt[guide_category %in% c("titration", "weismanTSS", "distalCRE") & GUIDE_STRAND == "-", .(chr, coord_end, coord_start, GUIDE_STRAND, guide_1, guide_class)], 
           old = c("coord_start", "coord_end", "GUIDE_STRAND"), new=c("coord_end", "coord_start", "strand")))
GDO_dt <- data.table(GDO_df)

# Create GRanges for guides coordinates
GDO_gr <- makeGRangesFromDataFrame(GDO_df, keep.extra.columns = T)


# In loop find overlaps between epi gr with guides gr
for (epi in names(EPI_l)) {
  
  # Find indexes of overlaping guides to peaks and 
  idx_dt <- as.data.table(GenomicRanges::findOverlaps(GDO_gr, EPI_l[[epi]], ignore.strand=TRUE, type = "any"))
  
  # Extract guide IDs and peak's signal value
  guide_1s_overlap <- unique(GDO_gr[idx_dt$queryHits]$guide_1)
  
  # If multiple peaks per genomic range (e.g. ATAC bed has mulitiple peaks per same coord range) in the epi GR, find the signal of the most significant peak
  if (length(guide_1s_overlap) < length(idx_dt$queryHits)) {
    signal_overlap <- foreach(i=unique(idx_dt$queryHits), .combine = c) %dopar% {
      peak_ids <- idx_dt[ queryHits == i, subjectHits]
      if (length(peak_ids) > 1) {
        dt <- as.data.table((EPI_l[[epi]][peak_ids, c("signalValue", "qValue")]))
        dt[qValue == dt[, max(qValue)], signalValue]
      } else{
        EPI_l[[epi]][peak_ids]$signalValue
      }
    } 
  } else {
    signal_overlap <- EPI_l[[epi]][idx_dt$subjectHits]$signalValue
  }
  
  # Incorporate as metadata into GDO dt
  GDO_dt[, paste0("inpeak_", epi) := ifelse(guide_1 %in% guide_1s_overlap, T, F), guide_1]
  GDO_dt[guide_1 %in% guide_1s_overlap,  paste0("peaksignal_", epi) := signal_overlap]
}

# Melt GDO_dt to have the 'inpeak' and 'peaksignal' variables on a single column
constcol <- c("chr", "coord_start", "coord_end", "strand", "guide_1", "guide_class")

tmp1 <- melt.data.table(GDO_dt[ , .SD, .SDcols = c(constcol, colnames(GDO_dt)[grepl("inpeak", colnames(GDO_dt))])], id.vars = constcol, value.name = "inpeak", variable.name = "epi")
tmp1[, peak_type := gsub("inpeak_", "", epi)]
tmp2 <- melt.data.table(GDO_dt[ , .SD, .SDcols = c("guide_1", colnames(GDO_dt)[grepl("peaksignal", colnames(GDO_dt))])], id.vars = "guide_1", value.name = "peaksignal", variable.name = "epi")
tmp2[, peak_type := gsub("peaksignal_", "", epi)]

GDO_melt_dt <- merge.data.table(tmp1[, !"epi", with=F], tmp2[, !"epi", with=F], by = c("guide_1", "peak_type"))
GDO_melt_dt[, dosage_gene := tstrsplit(guide_1, "_", keep = 1)][]

## Plot Frequency of overlap with any combination of those (ggupset plot)
plot_dt <- GDO_melt_dt[, .N, .(inpeak, peak_type)]
pA <- ggplot(plot_dt, aes(x = peak_type, y = N)) +
  geom_bar(stat = "identity", aes(fill=inpeak)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  ylab("number of sgRNAs") +
  scale_fill_manual("Overlap in peak?", values = list(`TRUE` = "grey30", `FALSE`= "grey80" )) +
  theme(legend.position = "bottom", legend.direction = "vertical")
pA


# Relationship btw the cis gene FC with the overlap between categories

# Add dosage gene's fold change info
GDO_melt_a_dt <- copy(GDO_melt_dt)
GDO_melt_i_dt <- copy(GDO_melt_dt)
GDO_melt_a_dt[, guide_crispr := paste0(guide_1, "-CRISPRa")]
GDO_melt_i_dt[, guide_crispr := paste0(guide_1, "-CRISPRi")]

GDO_FC_melt_dt <- unique(merge.data.table(x = rbind(GDO_melt_a_dt, GDO_melt_i_dt),
                                   y = DG_dt[, .(guide_crispr, cis_gene, cell_line, dosage_gene_log2FC)],
                                   by = "guide_crispr"))
GDO_FC_melt_dt[, inpeak_yn := ifelse(inpeak == T, "Yes", "No")]


wcxtest_dt <- foreach(cr = c("CRISPRa", "CRISPRi"), .combine=rbind) %do% {
  foreach(epi = names(EPI_l), .combine=rbind) %do% {
    pval_long <- wilcox.test(GDO_FC_melt_dt[cell_line == cr & peak_type == epi & inpeak == T, dosage_gene_log2FC],
                             GDO_FC_melt_dt[cell_line == cr & peak_type == epi & inpeak == F, dosage_gene_log2FC])$p.value
    data.table(cell_line = cr, 
               peak_type = epi,
               pval = format(pval_long, scientific = TRUE, digits = 3),
               col_test = ifelse(pval_long <= 0.05, T, F),
               dosage_gene_log2FC = ifelse(cr == "CRISPRa", GDO_FC_melt_dt[cell_line == "CRISPRa", max(dosage_gene_log2FC)], GDO_FC_melt_dt[cell_line == "CRISPRi", min(dosage_gene_log2FC)]))
    
    
  }
}

# Plot difference in FC when guide overlaps or not with epigenetic peak
pB <- ggplot(GDO_FC_melt_dt, aes(x = inpeak_yn, y=dosage_gene_log2FC)) +
  geom_violin(aes(fill=cell_line), alpha=0.5, show.legend = F) +
  geom_jitter(width = 0.1, alpha=0.2, color="grey20", shape = 20, size = 1) +
  geom_boxplot(width = 0.2, outlier.shape = NA, notch = T, aes(fill=cell_line), show.legend = F) +
  scale_fill_manual("cell line", values = col_crispr) +
  labs(x=paste0("Overlap with peak?"), y = "Cis gene log2(FC)") +
  facet_grid(cell_line ~ peak_type, scales = "free_y") +
  geom_text(data=wcxtest_dt, aes(x = 1.5, color = col_test, label=paste0("p = ", pval)), size=3, show.legend = F) +
  scale_color_manual("", values = list(`TRUE` = "grey30", `FALSE`= "grey70" )) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
pB


p <- plot_grid(pA, pB, align = "h", axis = "tb", nrow = 1, rel_widths = c(1/5, 4/5))
ggsave(file.path(plots_dir, "05b_03_FCvsEpiPeaks.pdf"), p, width = 8.75, height = 5)




