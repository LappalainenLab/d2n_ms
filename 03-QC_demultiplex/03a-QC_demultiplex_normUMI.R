##### D2N - QC, demultiplexing and UMI normalization #####
# (1) QC cDNA
# (2) QC GDOs
# (3) QC HTOs
# (4) Demultiplex cell lines
# (5) Merge lanes and normalize UMIs
# (6) Plot distribution of UMIs per TSS sgRNAs


# JDE, May 2022
# Last modified: May 2023


library(Seurat)
library(data.table)
library(ggridges)
library(foreach)
library(doParallel); registerDoParallel(detectCores()/2)
library(ggplot2); theme_set(theme_bw())


## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/03-QC_demultiplex/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


col_crispr = list(CRISPRa = "#FF7F00", CRISPRi = "#377EB8")




### (6) Plot distribution of UMIs per TSS sgRNAs

# !!! TEMP until rest of the code in place #
d2n <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_SeuratObj_PostQC_dCas9_NewMD.RDS") 
MD <- data.table(d2n@meta.data)
MD_clean <- MD[!is.na(ident), ]

ntc_guides = MD_clean[grepl("NTC", short_ID), unique(short_ID)]
control_guides <- MD_clean[guide_class == "TSS", unique(short_ID)]

plot_dt <- foreach(cnt_guide = control_guides, .combine = rbind) %do% {
  gRNA_vec <- c(cnt_guide, ntc_guides)
  dosage_gene = unlist(strsplit(cnt_guide, "_"))[[1]]
  dt <- do.call("rbind", lapply(gRNA_vec, function(gRNA){
    cells_idx <- MD_clean[short_ID == gRNA, L_cell_barcode]
    data.table(dosage_gene_expr = d2n[["RNA"]]@data[dosage_gene, cells_idx],
               cell_barcodes = cells_idx,
               cell_line = MD[L_cell_barcode %in% cells_idx, cell_line],
               guide = gRNA,
               perturbation = cnt_guide)
  }))
  dt[, is_cnt_guide := ifelse(guide == cnt_guide, T, F)]
  dt[guide == cnt_guide, guide := "Cis gene"]
}
ncells_dt <- plot_dt[, .N, .(guide, cell_line, perturbation)]

p <- ggplot(plot_dt, aes(x = dosage_gene_expr, y = guide)) +
  geom_density_ridges(scale = 2, size = 0.25, rel_min_height = 0.001, aes(fill = cell_line, alpha = is_cnt_guide), quantile_lines = TRUE, quantiles = 2, show.legend = F) +
  facet_grid(cell_line ~ perturbation, scales = "free_y") +
  geom_text(data = ncells_dt, aes(label = paste0("n = ", N), x=-2), size=2.5) +
  labs(x="Cis gene normalized expression", y="Cells with guides:") +
  scale_fill_manual("", values = col_crispr) +
  scale_alpha_manual("", values = c(0.5, 1)) +
  scale_x_continuous(limits = c(-3, ceiling(max(dt$dosage_gene_expr)))) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
p
ggsave(file.path(plots_dir, "03a_XX_DistUMItss.pdf"), p, width = 15, height = 4)

