##### D2N - QC, demultiplexing and UMI normalization #####
# (1) QC cDNA & targeted enrichment
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
dosage_genes = c("GFI1B", "NFE2", "MYB", "TET2")

### (1) QC cDNA


## Basic cDNA QC
# 1. Set lower UMI count and feature count thresholds (diff number of cells depending on 10x lane)
# 2. Remove cells with excess mito/genomic DNA ratio (remove dead cells)
nCount_low_thrs = 500
nFeature_low_thr = 50

# gene classes
crispr_genes <- c("CRISPRi", "CRISPRa")
dosage_genes <- c("GFI1B", "NFE2", "MYB", "TET2")


## Filtering cells by nCounts and nFeatures (use only filtered matrices)

cDNA_tg_L <- list() # keep full cDNA datasets for later
# Iterate over 10X lanes
for (L in 1:2) {
  # Load data and create Seurat object
  data_dir = paste0("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/006-fastq2features/cellranger/d2n_cDNA_CRISPR_L", 
                    L, "/outs/filtered_feature_bc_matrix")
  data <- Read10X(data.dir = data_dir)
  cDNA <- CreateSeuratObject(data[!(rownames(data) %in% crispr_genes),])
  
  # Add aditional metadata
  # complexity of cells
  cDNA$GenesPerUMI <- cDNA$nFeature_RNA/cDNA$nCount_RNA
  # cells passing filters based on minimum total TG UMI counts and feature count
  cDNA$passing_filter_1 <- cDNA$nCount_RNA >= nCount_low_thrs & cDNA$nFeature_RNA >= nFeature_low_thr 
  # remove cells that have to high total content of RNA
  nCount_high_thrs = quantile(cDNA$nCount_RNA[cDNA$passing_filter_1], 0.99)
  cDNA$passing_filter_2 <- cDNA$nCount_RNA <= nCount_high_thrs
  # cells passing all filters
  cDNA$passing_filter <- cDNA$passing_filter_1 & cDNA$passing_filter_2
  
  # Save object
  cDNA_tg_L[[L]] <- cDNA
  n_cells = sum(cDNA$passing_filter)
  pct_cells <- n_cells*100/nrow(cDNA@meta.data)
  
  p1 <- ggplot(cDNA@meta.data, aes(x=nCount_RNA)) + geom_histogram(bins = 100) +
    geom_vline(xintercept = nCount_low_thrs, color= "red") +
    geom_vline(xintercept = nCount_high_thrs, color= "gold") +
    ggtitle(paste0("L", L))
  ggsave(paste0("plots/20220411_cDNA_L", L, "_nCountDist.pdf"), p1, width = 6, height = 4)
  
  p2 <- ggplot(cDNA@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) +
    stat_bin_hex(bins = 50) + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    geom_vline(xintercept = nCount_low_thrs, color= "red") +
    geom_vline(xintercept = nCount_high_thrs, color= "gold") +
    geom_hline(yintercept = nFeature_low_thr, color = "red") +
    ggtitle(paste0("L", L, ", cells passing filters = ", n_cells, " (", round(pct_cells, 1), "%)"))
  ggsave(paste0("plots/20220411_cDNA_L", L, "_nCountVsnFeature.pdf"), p2, width = 6, height = 6)
  
  p2b <- ggplot(cDNA@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) +
    stat_bin_hex(bins = 50) + 
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "grey20", size=0.5) +
    geom_vline(xintercept = nCount_low_thrs, color= "grey20", linetype=2) +
    geom_vline(xintercept = nCount_high_thrs, color= "grey20", linetype=3) +
    geom_hline(yintercept = nFeature_low_thr, color = "grey20", linetype=2) +
    scale_fill_gradient("# cells", low = "grey90", high = "grey20") +
    theme(legend.position = c(0.8, 0.38)) +
    labs(x="Total RNA UMI counts", y="Number genes")
  p2b
  ggsave(paste0("plots/20220411_cDNA_L", L, "_nCountVsnFeature_Poster.pdf"), p2b, width = 3.2, height = 3)
  
  p3 <- ggplot(cDNA@meta.data, aes(x=GenesPerUMI)) + 
    geom_histogram(bins = 50, aes(y = ..density..*2), color="grey", fill="grey") +
    geom_density(aes(color=passing_filter)) +
    ggtitle(paste0("L", L)) +
    scale_x_log10() + ylab("density")
  ggsave(paste0("plots/20220411_cDNA_L", L, "_GenesPerUMI.pdf"), p3, width = 6, height = 4)
  
}

# Genes classes
target_genes <- rownames(cDNA_tg_L[[1]])[!grepl("CRISPR", rownames(cDNA_tg_L[[1]]))]


## Load raw data to filter based on mito DNA and analyse enrichment
# Iterate over lanes
cDNA_raw_L <- list()
for (L in 1:2) {
  # Load data
  data_dir = paste0("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/006-fastq2features/cellranger/d2n_cDNA_CRISPR_L", 
                    L, "/outs/raw_feature_bc_matrix")
  data <- Read10X(data.dir = data_dir)
  cDNA <- CreateSeuratObject(data)
  cDNA_raw_L[[L]] <- cDNA
  
  p1 <- ggplot(cDNA@meta.data, aes(x=nCount_RNA)) + geom_histogram(bins = 100) +
    geom_vline(xintercept = nCount_low_thrs, color= "red") +
    ggtitle(paste0("L", L)) +
    scale_x_log10()
  p1
  ggsave(paste0("plots/20220411_cDNAraw_L", L, "_nCountDist.pdf"), p1, width = 6, height = 4)
  
}


# Us non-targeted genes for calculating the
remaining_genes <- rownames(cDNA_raw_L[[1]])[!(rownames(cDNA_raw_L[[1]]) %in% target_genes)]




mitoRatio_high_thr = 0.2
cDNA_ntg_L <- list()
for (L in 1:2) {
  counts <- GetAssayData(cDNA_raw_L[[L]], slot="counts", assay="RNA")[remaining_genes, ]
  cDNA <- CreateSeuratObject(counts = counts)
  cDNA$GenesPerUMI <- cDNA$nFeature_RNA/cDNA$nCount_RNA
  cDNA$mitoRatio <- PercentageFeatureSet(object = cDNA, pattern = "^MT-", assay = "RNA")/100
  cDNA$passing_filter <- cDNA$mitoRatio <= mitoRatio_high_thr | is.na(cDNA$mitoRatio)
  
  cDNA_ntg_L[[L]] <- cDNA
  
  p1 <- ggplot(cDNA_ntg_L[[L]]@meta.data[!is.na(cDNA$mitoRatio),], aes(x=mitoRatio)) +
    geom_density() +
    geom_vline(xintercept = mitoRatio_high_thr, color="red") +
    scale_x_log10() +
    ggtitle(paste0("L", L, ", filtered cells = ", sum(!cDNA$passing_filter, na.rm = T) ))
  suppressWarnings(ggsave(paste0("plots/20220411_cDNAntg_L", L, "_mitoRatio.pdf"), p1, width = 6, height = 4))
}



## Analyse enrichment due to targeted sequencing

# Total UMIs captured per gene
DT_L <- list()
for (L in 1:2) {
  # get cells that pass all filters (TG RNA and NTG mitoRatio)
  pass_count_QC <- colnames(cDNA_tg_L[[L]][, cDNA_tg_L[[L]]$passing_filter])
  pass_count_mito_QC <- pass_count_QC[!(pass_count_QC %in% colnames(cDNA_ntg_L[[L]])[!(cDNA_ntg_L[[L]]$passing_filter)])]
  
  counts_tg <- GetAssayData(cDNA_tg_L[[L]], slot = "counts", assay = "RNA")[, pass_count_mito_QC]
  counts_ntg <- GetAssayData(cDNA_ntg_L[[L]], slot = "counts", assay = "RNA")[, pass_count_mito_QC]
  print(paste0("cells (L", L, ") = ", dim(counts_tg)[2]))
  
  dt <- rbind(data.table(gene = rownames(counts_tg), total_UMI = rowSums(counts_tg), gene_class = "targeted"),
              data.table(gene = rownames(counts_ntg), total_UMI = rowSums(counts_ntg), gene_class = "non-targeted"))
  DT_L[[L]] <- dt
}

DT <- merge.data.table(DT_L[[1]], DT_L[[2]][, .(gene, total_UMI)], suffixes = c("_L1", "_L2"), by="gene")

DT[, gene_class := factor(gene_class), levels(c("non-targeted", "targeted"))]


p <- ggplot(DT, aes(x=total_UMI_L1, y=total_UMI_L2)) +
  geom_abline(linetype=2) +
  geom_point(aes(color=gene_class, alpha=gene_class)) +
  scale_x_log10() + scale_y_log10() +
  theme_bw() + theme(legend.position = c(0.25, 0.8)) +
  scale_alpha_manual(values = c(0.25, 0.75), guide = "none") +
  ggrepel::geom_label_repel(data = DT[gene %in% c(dosage_genes, crispr_genes, "LHX3", "GAPDH"),], 
                            aes(label=gene, color=gene_class), 
                            size = 3,
                            min.segment.length = 0, 
                            seed = 2, 
                            box.padding = 1.5,
                            nudge_x = 0.5,
                            nudge_y = 0.1,
                            arrow = arrow(length = unit(0.010, "npc")),
                            segment.color = "grey30",
                            show.legend = F) +
  ggtitle(paste0("R = ", round(cor(DT$total_UMI_L1, DT$total_UMI_L2), 3)))
suppressWarnings(ggsave("plots/20220412_GeneTotalUMIscatter.pdf", p, width = 5, height = 4.8))

p_poster <- ggplot(DT, aes(x=total_UMI_L1, y=total_UMI_L2)) +
  geom_abline(linetype=2) +
  geom_point(aes(color=gene_class, alpha=gene_class)) +
  scale_x_log10() + scale_y_log10() +
  theme_bw() + theme(legend.position = c(0.25, 0.65), legend.background = element_blank()) +
  scale_alpha_manual(values = c(0.25, 0.75), guide = "none") +
  ggrepel::geom_label_repel(data = DT[gene %in% c(dosage_genes, "LHX3", "GAPDH"),],
                            aes(label=gene, color=gene_class),
                            size = 3,
                            min.segment.length = 0,
                            seed = 2,
                            box.padding = 1.5,
                            nudge_x = 0.5,
                            nudge_y = 0.1,
                            arrow = arrow(length = unit(0.010, "npc")),
                            segment.color = "grey30",
                            show.legend = F) +
  scale_color_manual("gene type", values = c("grey20", "#FF7F00")) +
  annotate("text", label = paste0("R = ", round(cor(DT$total_UMI_L1, DT$total_UMI_L2), 3)), x = 10, y = 10e6) +
  labs(x = "Total RNA UMI (Lane 1)", y = "Total RNA UMI (Lane 2)")
p_poster
suppressWarnings(ggsave("plots/20220412_GeneTotalUMIscatter_Poster.pdf", p_poster, width = 3.2, height = 3))

## Get final cDNA Seurat objects, run dimensionality reduction and a quick visualizationand before saving them
n_var_feat = 40
n_pcs = 20
n_dim_umap = 10
n_dim_tsne = 5

cDNA_final_L <- list()
for (L in 1:2) {
  pass_count_QC <- colnames(cDNA_tg_L[[L]][, cDNA_tg_L[[L]]$passing_filter])
  pass_count_mito_QC <- pass_count_QC[!(pass_count_QC %in% colnames(cDNA_ntg_L[[L]])[!(cDNA_ntg_L[[L]]$passing_filter)])]
  cDNA_tg_L[[L]]$passing_filter_final <- rownames(cDNA_tg_L[[L]]@meta.data) %in% pass_count_mito_QC
  
  d2n <- subset(cDNA_tg_L[[L]], subset = passing_filter_final == T)
  
  d2n <- NormalizeData(d2n)
  d2n <- FindVariableFeatures(d2n, selection.method = "vst", nfeatures = n_var_feat)
  # VariableFeaturePlot(d2n)
  d2n <- ScaleData(d2n, features = rownames(d2n))
  d2n <- RunPCA(d2n, features = VariableFeatures(object = d2n), npcs = n_pcs, approx = FALSE, verbose = F)
  # DimPlot(d2n, reduction = "pca")
  # DimHeatmap(d2n, dims = 1:20, cells = 500, balanced = TRUE)
  # d2n <- JackStraw(d2n, num.replicate = 100)
  # d2n <- ScoreJackStraw(d2n, dims = 1:20)
  # JackStrawPlot(d2n, dims = 1:20)
  # ElbowPlot(d2n, ndims = 20)
  d2n <- RunUMAP(d2n, dims = 1:n_dim_umap)
  p1 <- FeaturePlot(d2n, reduction = "umap", features = "nCount_RNA")
  ggsave(paste0("plots/20220413_cDNA_L", L, "_UMAP_nCount.pdf"), p1, width = 6, height = 4.8)
  p2 <- FeaturePlot(d2n,  reduction = "umap",features = dosage_genes)
  ggsave(paste0("plots/20220413_cDNA_L", L, "_UMAP_DosageGenes.pdf"), p2, width = 12, height = 9.6)
  d2n <- RunTSNE(d2n, dims = 1:n_dim_tsne)
  p3 <- FeaturePlot(d2n, reduction = "tsne", features = "nCount_RNA")
  ggsave(paste0("plots/20220413_cDNA_L", L, "_TSNE_nCount.pdf"), p3, width = 6, height = 4.8)
  p4 <- FeaturePlot(d2n, reduction = "tsne", features = dosage_genes)
  ggsave(paste0("plots/20220413_cDNA_L", L, "_TSNE_DosageGenes.pdf"), p4, width = 12, height = 9.6)
  
  cDNA_final_L[[L]] <- d2n
}

for (L in 1:2) {
  saveRDS(cDNA_final_L[[L]], file = paste0("processed_data/cDNA_L", L, "_postQC.RDS"))
}


### (2) QC GDOs

## Load data
cDNA_L1 <- readRDS("processed_data/cDNA_dCas9_L1_postQC.RDS")
cDNA_L2 <- readRDS("processed_data/cDNA_dCas9_L2_postQC.RDS")
cDNA_L <- list(cDNA_L1, cDNA_L2)
guides_features <- fread("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/001-library_design/processed_data/guides/D2N_ALLguides_features_sortID.csv")
guides_features[, id := gsub("_", ".", ID), ID]


## Load GDO data, get guide calls, update metadata and merge with cDNA assay
d2n_L <- list()
for (L in 1:2) {
  # load only the GDO data
  GDO_raw <- Read10X(data.dir = paste0("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/006-fastq2features/cellranger/d2n_GDOs_L", L, "/outs/filtered_feature_bc_matrix"),
                     gene.column = 1, strip.suffix = T)[["CRISPR Guide Capture"]]
  # change guide id so that they don't have underscores "_" (not compatible with Seurat)
  old_guide_id <- rownames(GDO_raw)
  rownames(GDO_raw) <- gsub("_", ".", old_guide_id)
  
  # subset those cells that passed cDNA QC
  cells_cDNA <- colnames(cDNA_L[[L]])
  GDO <- CreateSeuratObject(GDO_raw[, cells_cDNA], assay = "GDO") 
  GDO <- NormalizeData(GDO, normalization.method = "CLR")
  
  # Create new metadata with the guide calling information
  guide_per_cell <- as.data.frame(fread(paste0("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/006-fastq2features/cellranger/d2n_GDOs_L", L, "/outs/crispr_analysis/protospacer_calls_per_cell.csv")))
  rownames(guide_per_cell) <- gsub("-1", "", guide_per_cell$cell_barcode)
  guide_per_cell$guide_call <- gsub("_", ".", guide_per_cell$feature_call)
  GDO <- AddMetaData(GDO, guide_per_cell[, c("num_features", "guide_call", "num_umis")])
  GDO$num_features[is.na(GDO$num_features)] = 0
  GDO$num_guides <- GDO$num_features
  GDO$num_guides[GDO$num_features >= 3] = ">2" 
  GDO$num_guides <- factor(GDO$num_guides, levels = c("0", "1", "2", ">2"))
  GDO$cell_barcode = colnames(GDO)
  
  d2n <- GDO
  d2n[["RNA"]] <- CreateAssayObject(cDNA_L[[L]][["RNA"]]@counts)
  d2n <- NormalizeData(d2n, assay = "RNA")
  d2n_L[[L]] <- d2n
}

## GDOs stats

# Cells with guides and number of guides per cell
for (L in 1:2) {
  DT <- data.table(d2n_L[[L]]@meta.data)
  
  # (1) % of cells w/o, w/ 1 or more guides
  plot_dt <- data.table(table(DT$num_features))
  plot_dt[, num_guides_cell := as.numeric(V1)]
  
  p1 <- ggplot(plot_dt, aes(x=num_guides_cell, y=N)) +
    geom_bar(stat = "identity", fill = "grey60") +
    geom_text(aes(label = N), nudge_y = 750, angle=60, nudge_x = 0.2) +
    ylab("num. of cells") + xlab("num. guides per cell") +
    scale_y_continuous(limits = c(0,(max(plot_dt$N)+1000)))
  if (max(plot_dt$num_guides_cell) > 11){
    p1 <- p1 + scale_x_continuous(breaks = c(0:10, seq(from = 12, to = max(plot_dt$num_guides_cell), by = 2)))
  } else {
    p1 <- p1 + scale_x_continuous(breaks = c(0:max(plot_dt$num_guides_cell)))
  }
  p1
  ggsave(paste0("plots/20220708_GDO_dCas9_L", L, "_GuidesPerCell.pdf"), p1, width = 3.5, height = 4)
  
  # (2) Number of cells per guide
  # Only singlets
  plot_dt_all <- setnames(data.table(table(DT[num_features > 0, guide_call])), c("guide_call", "N"))
  plot_dt_all[, num_guides_cell := str_count(pattern = "\\|", guide_call)+1, guide_call]
  plot_dt_singlets <- plot_dt_all[num_guides_cell == 1, ]
  plot_dt_singlets[, dosage_gene := tstrsplit(guide_call, "\\.", keep = 1)]
  plot_dt_singlets[grepl("^NTC-", dosage_gene), dosage_gene := "NTC"]
  plot_dt_singlets[, dosage_gene := factor(dosage_gene, levels = c(dosage_genes, "NTC"))]
  
  p2 <- ggplot(plot_dt_singlets, aes(x=N)) +
    geom_histogram(aes(color = dosage_gene, fill=dosage_gene), show.legend = F, bins = 50) +
    facet_wrap(dosage_gene ~ ., scales = "fixed", ncol = 1) +
    labs(x = "# cells per guide", y = "# guides") +
    scale_color_brewer("gRNA gene", palette = "Set1") +
    scale_fill_brewer("gRNA gene", palette = "Set1") +
    ggtitle(paste0("L", L, ", cells with 1 guide only"))
  p2
  ggsave(paste0("plots/20220708_GDO_dCas9_L", L, "_CellsPerGuide.pdf"), p2, width = 6, height = 12)
  
  # (3) GDO counts and normalized expression given the number of guides
  # Categorize those >2 guides into single category
  
  p3 <- ggplot(DT, aes(x=factor(num_features), y=nCount_GDO)) +
    geom_violin(scale = "width", fill="grey90") +
    geom_boxplot(outlier.colour = NA, width = 0.2) +
    scale_y_log10() +
    labs(y="Total GDO UMI counts", x="Number guides per cell")
  p3
  suppressWarnings(ggsave(paste0("plots/20220708_GDO_dCas9_L", L, "_nCountsGDOperNumGuides.pdf"), p3, width = 6, height = 4))
  
}


## Visually check guide calling

col_genes = RColorBrewer::brewer.pal(5, "Set1")
names(col_genes) <- c("MYB", "TET2", "NFE2", "GFI1B", "NTC")
# (3) Generate RidgePlots for all guides expression where one lead guide has been selected

# Iterate over lane and single guide RNA called cells
for (L in 1:2) {
  
  d2n_singlets <- subset(d2n_L[[L]], subset = num_features == 1)
  all_guides <- rownames(d2n_singlets) 
  
  for (lead_guide in all_guides) {
    
    # select all cells called as having that single lead guide
    MD <- as.data.table(d2n_singlets@meta.data)
    idx_cells <- MD[guide_call == lead_guide, cell_barcode]
    
    # Get normalized expression for all guides on those cells subset
    plot_dt <- merge.data.table(do.call("rbind", lapply(all_guides, function(x){
      data.table(guide_ident = x,
                 norm_gdo_expression = d2n_singlets[["GDO"]]@data[x, idx_cells])
    })), guides_features[, .(id, short_ID, guide_category)], by.x = "guide_ident", by.y = "id")
    plot_dt[,  is_called := ifelse(guide_ident == lead_guide, T, F) ]
    plot_dt[, gene := tstrsplit(short_ID, "_", keep = 1)]
    plot_dt[, gene := factor(gene, levels = c(dosage_genes, "NTC"))]
    
    
    lead_guide_short_id = unique(plot_dt[is_called == T, short_ID])
    p <- ggplot(plot_dt, aes(x = norm_gdo_expression, y = short_ID)) +
      geom_density_ridges(aes(fill=is_called, color=gene), show.legend = T,
                          scale = 6, size = 0.5, rel_min_height = 0.01) +
      scale_color_manual("sgRNA gene", values = col_genes) +
      scale_fill_manual("called sgRNA", values = c("grey80", "grey20")) +
      labs(x="sgRNA norm. expression", y="sgRNAs") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      ggtitle(paste0(lead_guide_short_id, " (", length(idx_cells), " cells), L", L, " - ", plot_dt[is_called == T, unique(guide_category)]))
    suppressMessages(ggsave(paste0("plots/ridge_plots_per_guide_cas9/", lead_guide_short_id, "_L", L, ".pdf"), p, width = 3.5, height = 5))
  }
}

### Total UMI across cells vs initial guide library reads distribution

# load read distribution data of the sgRNA library for d2n 
lib_dt <- fread("../005-sgRNA_library_distribution/processed_data/d2n-J.counts_h.txt", col.names = c("ID", "read_counts"))
lib_dt[, guide := gsub("_", ".", ID), ID]

dt_L <- list()
for (L in 1:2) {
  # Usin all cells, suming all the UMI counts per guide
  dt <- merge.data.table(lib_dt, data.table(guide = rownames(d2n_L[[L]]),
                                            sc_counts = rowSums(d2n_L[[L]]),
                                            lane = paste0("L", L),
                                            cells = "All"),
                         by = "guide")
  dt_L[[L]] <- dt
}
plot_dt <- rbindlist(dt_L)
plot_dt[, gene := tstrsplit(guide, "\\.", keep = 1), guide]
plot_dt[grepl("NTC", gene), gene := "NTC", guide]
plot_dt[, gene := factor(gene, levels = c(dosage_genes, "NTC"))]


dt_cor <- setnames(plot_dt[, cor(read_counts, sc_counts), by="lane"], c("lane", "r"))
p6 <- ggplot(plot_dt, aes(x=read_counts, y=sc_counts)) +
  geom_abline(linetype=2, size=0.5, color="grey75") +
  scale_x_log10() + scale_y_log10() +
  geom_smooth(formula = 'y ~ x', method = "lm", size=0.75, aes(group=gene, color=gene), show.legend = F, se = T, alpha = 0.1) +
  geom_point(aes(color=gene)) +
  labs(x = "Read counts per guide (gRNA plasmid library)", y = "Total UMI counts per guide (scRNA-seq GDOs)") +
  facet_wrap(lane ~ ., ncol = 1) +
  geom_text(data = dt_cor, aes(label = paste0("r = ", round(r, 2)), y=6000, x=750)) +
  scale_color_brewer("gRNA gene", palette = "Set1")
p6
ggsave("plots/20220708_GDO_dCas9_LibDistVsTotalUMI.pdf", p6, height = 8, width = 6)


### Expression dCas9 across cells with different number of guides


plot_dt <- do.call("rbind", lapply(1:2, function(L){
  dt <- rbind(data.table(norm_RNA = d2n_L[[L]][["RNA"]]@data["CRISPRi",],
                         cell_barcode = colnames(d2n_L[[L]]),
                         class = "CRISPRi"),
              data.table(norm_RNA = d2n_L[[L]][["RNA"]]@data["CRISPRa",],
                         cell_barcode = colnames(d2n_L[[L]]),
                         class = "CRISPRa"),
              data.table(norm_RNA = colMeans(d2n_L[[L]][["RNA"]]@data),
                         cell_barcode = colnames(d2n_L[[L]]),
                         class = "Mean all genes")
  )
  mdt <- data.table(d2n_L[[L]]@meta.data) 
  DT <- merge.data.table(mdt[, .(cell_barcode, num_features)], dt, by="cell_barcode")
  DT$lane = paste0("Lane ", L)
  DT
}))
p7 <- ggplot(plot_dt, aes(x=factor(num_features), y=norm_RNA)) +
  geom_violin(scale = "width", fill="grey90") +
  geom_boxplot(outlier.colour = NA, width = 0.2) +
  facet_grid(lane ~ class) +
  labs(y="Normalized RNA expression", x="Number guides per cell")
p7
ggsave("plots/20220714_GDO_dCas9_RNAexprVsNumGuides.pdf", p7, height = 5, width = 9)


## Save the QC seurat objects with the cDNA and GDOs containing all cells passing cDNA QC
for (L in 1:2) {
  saveRDS(d2n_L[[L]], file = paste0("processed_data/GDOs_cDNA_dCas9_L", L, "_postQC.RDS"))
}

### (3) QC HTOs 
## Load data
# HTO mapping
mapping_hto <- merge.data.table(fread("../006-fastq2features/processed_data/hto_mapping.tsv",  col.names = c("HTO", "id"), header = F),
                                fread("../006-fastq2features/processed_data/d2n_HTOs.txt")[, .(HTO, cell_line)],
                                by = "HTO")

# cDNA and GDO data
d2n_L1 <- readRDS("processed_data/GDOs_cDNA_dCas9_L1_postQC.RDS")
d2n_L2 <- readRDS("processed_data/GDOs_cDNA_dCas9_L2_postQC.RDS")
d2n_L <- list(d2n_L1, d2n_L2)

# Guides metadata
guides_features <- fread("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/001-library_design/processed_data/guides/D2N_ALLguides_features_sortID.csv")
guides_features[, id := gsub("_", ".", ID), ID]


### Distribution of HTO UMIs per unique HTO
for (L in 1:2){
  mtx <- tximport(files = paste0("../006-fastq2features/salmon/d2n_HTO_L", L, "/alevin/quants_mat.gz"), type="alevin")[["counts"]]
  
  S <- d2n_L[[L]]
  cells <- colnames(S)
  
  S[["HTO"]] <- CreateAssayObject(counts = mtx[,cells])
  S <- NormalizeData(S, assay = "HTO", normalization.method = "CLR")
  
  DT <- data.table(S@meta.data)
  
  p1 <- ggplot(DT, aes(x=factor(num_features), y=nCount_HTO)) +
    geom_violin(scale = "width", fill="grey90") +
    geom_boxplot(outlier.colour = NA, width = 0.2, notch = T) +
    labs(y="Total HTO UMI counts", x="Number guides per cell") +
    scale_y_continuous(limits = c(100, 2500))
  p1
  suppressWarnings(ggsave(paste0("plots/20220718_HTO_dCas9_L", L, "_nCountsHTOperNumGuides.pdf"), p1, width = 6, height = 4))
  
  plot_dt <- melt(data.table(t(S[["HTO"]]@data)), measure.vars = gsub("_", "-", mapping_hto$id), variable.name = "HTO", value.name = "HTO_norm_expr")
  
  p2 <- ggplot(plot_dt, aes(x = HTO_norm_expr)) + 
    geom_histogram(bins = 50) +
    facet_grid(HTO ~ .) +
    scale_x_continuous(name = "HTO normalized expression", limits = c(-0.1,6))
  p2
  suppressWarnings(ggsave(paste0("plots/20220718_HTO_dCas9_L", L, "_DistHTONormExpr.pdf"), p2, width = 5, height = 6))
  
}

# # Use only single guides cells, and run HTODemux to select the best Thrshold that gives biggest num of singlets.
# selected_thrs_L <- list()
# for (L in 1:2) {
# 
#   # cDNA and GDOs
#   S <- subset(d2n_L[[L]], subset = num_guides == "1")
#   cells <- colnames(S)
#   
#   # Load alevin data 
#   mtx <- tximport(files = paste0("../006-fastq2features/salmon/d2n_HTO_L", L, "/alevin/quants_mat.gz"), type="alevin")[["counts"]]
#   mtx_subset <- mtx[, cells]
#   
#   # Create new HTO assay
#   S[["HTO"]] <- CreateAssayObject(counts = mtx[,cells])
#   DefaultAssay(S) <- "HTO"
#   S <- NormalizeData(S, assay = "HTO", normalization.method = "CLR")
#   S <- ScaleData(S, assay = "HTO", verbose = F)
#   
#   quant_thrs_seq = seq(0.20,0.995,0.01)
#   
#   plot_dt <- do.call("rbind", mclapply(quant_thrs_seq, function(x){
#     tmp <- HTODemux(S, assay = "HTO", positive.quantile = x, nsamples = 1000, verbose = F)
#     DT_1 <- data.table(tmp@meta.data)
#     DT_2 <- DT_1[, .N, HTO_classification.global]
#     DT_2[, quant_thrs := x]
#     return(DT_2)
#   }, mc.cores = ncores))
#   
#   p <- ggplot(plot_dt, aes(x=quant_thrs, y=N)) + 
#     geom_point(aes(color = HTO_classification.global)) +
#     labs(y = "# of cells in HTO category", x = "HTODemux Positive quantile threshold") +
#     scale_color_brewer("HTO global class", palette = "Set1") +
#     theme(legend.position = "bottom", legend.direction = "horizontal")
#   p
#   ggsave(paste0("plots/20220718_HTO_dCas9_HTODemux_L", L, "_GuidesSinglets_NumHTOGlobalClass.pdf"), p, width = 9, height = 6)
#   
#   selected_thrs <- plot_dt[N == max(plot_dt[HTO_classification.global == "Singlet", N]), quant_thrs]
#   selected_thrs_L[[L]] <- selected_thrs
# }

selected_thrs_L <- list(0.98, 0.98)

## Manually correct the demultiplexing 
col_HTO_id_all <- c(rev(c("#90171f", "#bc1e28", "#dd303b", "#e55c64", "#134b13", "#1e741e", "#289c28", "#33c533")), "grey70", "grey50")
names(col_HTO_id_all) <- c(gsub("_", "-", mapping_hto[, id]), "Doublet", "Negative")


d2n_curated_L <- list()
d2n_mix_L <- list()
for (L in 1:2) {
  # cDNA and GDOs
  S <- subset(d2n_L[[L]], subset = num_guides == "1")
  cells <- colnames(S)
  
  # Load alevin data 
  mtx <- tximport(files = paste0("../006-fastq2features/salmon/d2n_HTO_L", L, "/alevin/quants_mat.gz"), type="alevin")[["counts"]]
  mtx_subset <- mtx[, cells]
  
  # Create new HTO assay
  S[["HTO"]] <- CreateAssayObject(counts = mtx[,cells])
  DefaultAssay(S) <- "HTO"
  S <- NormalizeData(S, assay = "HTO", normalization.method = "CLR", margin = 2)
  S <- ScaleData(S, assay = "HTO", verbose = F)
  
  d2n_curated_L[[L]] <- S
  
  # HTODemux
  S <- HTODemux(S, assay = "HTO", positive.quantile = selected_thrs_L[[L]], nsamples = 1000, verbose = F)
  
  # UMAP and clustering
  VariableFeatures(S) <- rownames(S[["HTO"]]@counts)
  S <- RunPCA(S, reduction.name = "hto.pca", reduction.key = "HPC_", verbose = F)
  S <- RunUMAP(S, reduction = "hto.pca", dims = 1:7, reduction.name = "hto.umap", reduction.key = "HUMAP_", verbose = F)
  
  p <- DimPlot(S,reduction = "hto.umap", group.by = "hash.ID") + 
    ggtitle(paste0("Single guide cells, HTO UMAP (L", L, "), HTODmux ", selected_thrs_L[[L]])) +
    scale_color_manual("HTO class", values = col_HTO_id_all)
  ggsave(paste0("plots/20220718_HTO_dCas9_L", L, "_", selected_thrs_L[[L]], "_UMAP.pdf"), p, width = 6, height = 5)
  
  DefaultAssay(S) <- "RNA"
  p2 <- FeaturePlot(S, reduction = "hto.umap", features = c("CRISPRi", "CRISPRa"))
  ggsave(paste0("plots/20220718_HTO_dCas9_L", L, "_", selected_thrs_L[[L]], "_UMAP_ColorCRISPRiaRNA.pdf"), p2, width = 10, height = 5)
  
  # Create a new assay combining HTOs and CRISPRi expression to normalize, scale and cluster data based on both RNA and HTO quantitative data
  HTOs_dCas9_mtx <- rbind(S[["HTO"]]@counts, S[["RNA"]]@counts["CRISPRi", ])
  rownames(HTOs_dCas9_mtx) <- c(rownames(S[["HTO"]]), "CRISPRi")
  S[["MIX"]] <- CreateAssayObject(counts = HTOs_dCas9_mtx)
  S <- NormalizeData(S, assay = "MIX", normalization.method = "CLR", margin = 2)
  S <- ScaleData(S, assay = "MIX", verbose = F)
  VariableFeatures(S, assay = "MIX") <- rownames(S[["MIX"]])
  S <- RunPCA(S, assay = "MIX", reduction.name = "mix.pca", reduction.key = "MIX_", verbose = F)
  S <- RunUMAP(S, assay = "MIX", reduction = "mix.pca", dims = 1:8, reduction.name = "mix.umap", reduction.key = "MUMAP_", verbose = F)
  
  DefaultAssay(S) <- "MIX"
  p3 <- DimPlot(S, reduction = "mix.umap", group.by = "hash.ID") + 
    scale_color_manual("HTO class", values = col_HTO_id_all) +
    ggtitle("UMAP on HTO + CRISPRi counts")
  ggsave(paste0("plots/20220718_MIX_dCas9_L", L, "_UMAP_HashID.pdf"), p3, width = 6, height = 5)
  
  DefaultAssay(S) <- "RNA"
  p4i <- FeaturePlot(S, reduction = "mix.umap", features = "CRISPRi")
  ggsave(paste0("plots/20220718_MIX_dCas9_L", L, "_UMAP_FeatureCRISPRi.pdf"), p4i, width = 5.5, height = 5)
  p4a <- FeaturePlot(S, reduction = "mix.umap", features = "CRISPRa")
  ggsave(paste0("plots/20220718_MIX_dCas9_L", L, "_UMAP_FeatureCRISPRa.pdf"), p4a, width = 5.5, height = 5)
  
  DefaultAssay(S) <- "HTO"
  p5 <- FeaturePlot(S, reduction = "mix.umap", features = rownames(S[["HTO"]]))
  ggsave(paste0("plots/20220718_MIX_dCas9_L", L, "_UMAP_FeatureHTOs.pdf"), p5, width = 10, height = 10)
  
  Temp <- S
  DefaultAssay(Temp) <- "MIX"
  pL <- lapply(seq(0.05, 0.45, 0.05), function(x){
    Temp <- FindNeighbors(object = Temp, reduction = 'mix.umap', dims = 1:2)
    Temp <- FindClusters(object = Temp, verbose = FALSE, algorithm = 3, resolution = x)
    DimPlot(Temp, reduction = "mix.umap", group.by = "seurat_clusters", label = TRUE) + ggtitle(paste0("Resolultion = ", x))
  })
  p6 <- do.call("grid.arrange", c(pL, ncol=3))
  ggsave(paste0("plots/20220718_MIX_dCas9_L", L, "_UMAP_VariableResolution.pdf"), p6, width = 10, height = 10)
  
  d2n_mix_L[[L]] <- Temp
}

### Manually annotate the cell-lines based on the curated clustering & remove putative missclassification using CIRSPRi expression
mix_res_thrs = list(0.05, 0.05)
clust2crispr = list(L1 = list(a = c(2,3,4,5),
                              i = c(0,1)),
                    L2 = list(a = c(2,3,4,5),
                              i = c(0,1)))

d2n_annotated_L <- list()
for (L in 1:2) {
  
  Temp <- d2n_mix_L[[L]]
  
  # Find clusters using the mixed counts from HTOs and CRISPRi RNA
  Temp <- FindNeighbors(object = Temp, reduction = 'mix.umap', dims = 1:2)
  Temp <- FindClusters(object = Temp, verbose = FALSE, algorithm = 3, resolution = mix_res_thrs[[L]])
  
  # Assign clusters to cells
  Temp_MD <- data.table(Temp@meta.data)
  cells_a = Temp_MD[ seurat_clusters %in% clust2crispr[[paste0("L", L)]][["a"]], cell_barcode]
  cells_i = Temp_MD[ seurat_clusters %in% clust2crispr[[paste0("L", L)]][["i"]], cell_barcode]
  
  # Add cell line in metadata
  S <- d2n_curated_L[[L]]
  S$cell_line = "Undefined"
  S$cell_line[S$cell_barcode %in% cells_a] = "CRISPRa"
  S$cell_line[S$cell_barcode %in% cells_i] = "CRISPRi"
  
  d2n_annotated_L[[L]] <- S
}


### Remove top % CRISPRa cells expressing CRISPRi gene, and remove bottom % CRISPRi cells expressing CRISPRi gene
CRISPRi_expr_dt <- do.call("rbind", lapply(1:2, function(L){
  data.table(class_cell = d2n_annotated_L[[L]]@meta.data$cell_line,
             cell_barcode = d2n_annotated_L[[L]]@meta.data$cell_barcode,
             CRISPRi_norm_expr = d2n_annotated_L[[L]][["RNA"]]@data["CRISPRi", ],
             lane = paste0("L", L))
}))
thrs_dt <- rbind(cbind(CRISPRi_expr_dt[class_cell == "CRISPRi", quantile(CRISPRi_norm_expr, 0.05), lane], data.table(class_cell = "CRISPRi")), 
                 cbind(CRISPRi_expr_dt[class_cell == "CRISPRa", quantile(CRISPRi_norm_expr, 0.95), lane], data.table(class_cell = "CRISPRa")))


p <- ggplot(CRISPRi_expr_dt, aes(x=CRISPRi_norm_expr)) +
  geom_histogram(bins = 30) +
  facet_grid(class_cell ~ lane) +
  geom_vline(data = thrs_dt, aes(xintercept = V1), linetype = 2)
p
ggsave("plots/20220719_MIX_dCas9_CRISPRiExprDist_OutlierRemoval.pdf", p, width = 6, height = 3.5)


d2n_final_L <- list()
for (L in 1:2) {
  
  cells_i <- CRISPRi_expr_dt[lane == paste0("L", L) & class_cell == "CRISPRi" & CRISPRi_norm_expr >= thrs_dt[lane == paste0("L", L) & class_cell == "CRISPRi", V1], cell_barcode] 
  cells_a <- CRISPRi_expr_dt[lane == paste0("L", L) & class_cell == "CRISPRa" & CRISPRi_norm_expr <= thrs_dt[lane == paste0("L", L) & class_cell == "CRISPRa", V1], cell_barcode] 
  
  # Plot MIX UMAPs but 
  Temp <- subset(d2n_mix_L[[L]], subset = cell_barcode %in% c(cells_a, cells_i))
  DefaultAssay(Temp) <- "RNA"
  p <- FeaturePlot(Temp, reduction = "mix.umap", features = "CRISPRi")
  ggsave(paste0("plots/20220718_MIX_dCas9_L", L, "_UMAP_FeatureCRISPRi_OutliersCiRemoved.pdf"), p, width = 5.5, height = 5)
  
  d2n_final_L[[L]] <- subset(d2n_annotated_L[[L]], subset = cell_barcode %in% c(cells_a, cells_i))
}

### Merge and save final processed dataset
d2n_final <- merge(x = d2n_final_L[[1]],
                   y = d2n_final_L[[2]],
                   add.cell.ids = c("L1", "L2"), project = "d2n")
DefaultAssay(d2n_final) <- "RNA"
d2n_final <- NormalizeData(d2n_final, verbose = F)


### Final stats 
MD <- data.table(d2n_final@meta.data)
MD[, .N, cell_line]
# 1:   CRISPRa  9354
# 2:   CRISPRi 10647


## Save data ##
saveRDS(d2n_final, file = paste0("processed_data/HTOs_GDOs_cDNA_dCas9_merged_postQC.RDS"))
saveRDS(d2n_final_L[[1]], file = paste0("processed_data/HTOs_GDOs_cDNA_dCas9_merged_postQC_L1.RDS"))
saveRDS(d2n_final_L[[2]], file = paste0("processed_data/HTOs_GDOs_cDNA_dCas9_merged_postQC_L2.RDS"))



### (6) Plot distribution of UMIs per TSS sgRNAs

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

