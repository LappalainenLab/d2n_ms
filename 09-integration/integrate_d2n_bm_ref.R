library(here)
here::i_am("integrate_d2n_aml_ref.R")

#specify monocle lib path
#.libPaths( c( "/proj/lappalainen_lab1/users/marii/chip_seq_ann" , .libPaths() ) )


library(stringr)
library(monocle3)
library(Seurat)
library(corrplot)
library(cowplot)
library(ggplot2)
library(SeuratWrappers)
library(dplyr)
library(data.table)
library(Signac)
library(Matrix)
library(patchwork)
library(destiny)
library(viridis)
library(geosphere)
library(lsa)
library(foreach)
library(doParallel)	
library(tidyr)
library(gtable)
library(gridExtra)

registerDoParallel(cores=20)

##------------------------------------------Setting paths---------------------------------------
folder_name = basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)

##------------------------------------------Load reference dataset------------------------------
reference =  readRDS(paste0(data_dir, "/bm_ref.rds"))
reference <- subset(reference, subset = celltype.l1 %in% c('DC', 'HSPC', 'Mono', 'NK'))

##------------------------------------------Load and preprocess D2N dataset---------------------
d2n_full <- readRDS(file = paste0(data_dir, "/d2n_SeuratObj_PostQC_dCas9_NewMD.RDS"))
counts <- GetAssayData(d2n_full, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('CRISPRi', 'CRISPRa'))),]
d2n <- subset(d2n_full, features = rownames(counts))
d2n <- SCTransform(
  object = d2n,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference),
  reference.SCT.model = reference[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 20000,
  n_genes = 92,
  do.correct.umi = F,
  do.scale = FALSE,
  do.center = TRUE
)
VariableFeatures(d2n) = rownames(d2n[["RNA"]])

##------------------------------------------Integration------------------------------------------
anchors <- FindTransferAnchors(
  reference = reference, 
  query = d2n, 
  k.filter = NA, 
  #reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference), VariableFeatures(object = d2n)),
  dims = 1:50, 
  n.trees = 20,
  mapping.score.k = 100
)

d2n <- MapQuery(
  anchorset = anchors,
  query = d2n,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2"),
  reference.reduction = "refDR", 
  reduction.model = "refUMAP")

###------------------------------------------UMAP for shared embeddings---------------------------
d2n$celltype.l2 <- d2n$predicted.celltype.l2
d2n$group <- 'query'
reference$group <- 'reference'

emb.merge <- merge( reference[['refDR']], d2n[['ref.refDR']])
emb.merge@assay.used <- 'refAssay'
DefaultAssay(reference) <- DefaultAssay(d2n) <- 'refAssay'


obj.merge <- merge(DietSeurat(reference, assays = 'refAssay'), DietSeurat(d2n, assays = 'refAssay') )
obj.merge[['refDR']] <- emb.merge
obj.merge <- RunUMAP(obj.merge, dims = 1:30, reduction="refDR")
d2n_types = unique(d2n$celltype.l2)
ref_types = setdiff(unique(reference$celltype.l2), d2n_types)
obj.merge$celltype.l2 = factor(obj.merge$celltype.l2, levels=c(d2n_types, ref_types))

cols = DiscretePalette(24)

#------------------------------------------Panel A---------------------------------------------------
p1 <- DimPlot(obj.merge, cols = alpha(cols, 0.33), group.by = 'celltype.l2',  pt.size = 0.01, raster = F) + 
	NoLegend() + 
	ggtitle("All") +
	xlab("") + 
        ylab("UMAP_2")
p2 <- DimPlot(obj.merge,  cols = alpha(cols, 0.33), group.by = 'celltype.l2', cells = Cells(d2n),  pt.size = 0.01, raster = F) + 
	NoLegend() + 
	ggtitle("Perturbed") + 
	xlab("UMAP_1") +
	ylab("")
p3 <- DimPlot(obj.merge,   cols = alpha(cols, 0.33), group.by = 'celltype.l2', cells = Cells(reference),  pt.size = 0.01, raster = F) + 
	NoLegend() + 
	ggtitle("Reference") +
	xlab("") + 
        ylab("")
p4 =  DimPlot(obj.merge,  cols = alpha(cols, 0.33), group.by = 'celltype.l2', cells = Cells(reference), raster = F) + guides(colour = guide_legend(override.aes = list(size=2), nrow=3)) +
														      theme(legend.position= "bottom",
                                                                                                                                legend.text = element_text(size = 4),
                                                                                                                                legend.key.size = unit(2, 'mm'),
																legend.spacing.x = unit(1, 'mm'),
																legend.spacing.y = unit(1, 'mm'))

design = "
 111111112222222233333333
 111111112222222233333333
 111111112222222233333333
 111111112222222233333333
 111111112222222233333333
 111111112222222233333333
 111111112222222233333333
 111111112222222233333333
 111111112222222233333333
 #44444444444444444444444
"
leg = gtable_filter(ggplot_gtable(ggplot_build(p4)), "guide-box")
svg(paste0(plots_dir, "/Panel_A_BM_dataset.svg"), width = 118.8, height = 59.4,  res=720, units= "mm")
p1 + p2 + p3 + leg + plot_layout(design = design)&theme(axis.text = element_text(size = 4),
                                                                axis.title = element_text(size = 5),
                                                                plot.title = element_text(size = 7))
dev.off()


#------------------------------------------Diffusion map analysis-----------------------------------
df.map <- DiffusionMap(data = obj.merge[['refDR']]@cell.embeddings[,1:30])
obj.merge[['DC']] <- CreateDimReducObject(embeddings = df.map@eigenvectors*1000, key = 'DC_', assay = 'refAssay')

obj.merge.filt_cells = subset(reference, subset = celltype.l2 %in% d2n@meta.data$`predicted.celltype.l2`)

obj.merge@meta.data$cell_name = rownames(obj.merge@meta.data)

obj.merge.query <- subset(obj.merge, cells = Cells(d2n))
d2n = subset(obj.merge, cells = Cells(obj.merge.query))
d2n[['DC']] <- obj.merge.query[['DC']]


##------------------------------------------Pseudotime analysis---------------------------------------

###-----------------------------------------Erythroid Cells-------------------------------------------

eryth.monocle = subset(d2n, subset = (celltype.l2 == "Early Eryth") | (celltype.l2 == "Late Eryth"))
eryth.monocle[['umap']] = NULL
eryth.monocle[['umap']] <- eryth.monocle[['DC']]
eryth.monocle[['umap']]@cell.embeddings <- eryth.monocle[['umap']]@cell.embeddings[,1:2]

##------------------------------------------Select root cells------------------------------------------
early_eryth = subset(eryth.monocle, subset = (celltype.l2 == "Early Eryth"))[['umap']]@cell.embeddings[,1:2]
lim_x = quantile(early_eryth[, 1], probs = c(0.05, 0.95))
lim_y = quantile(early_eryth[, 2], probs = c(0.25, 0.65))

as.data.table(early_eryth) %>% filter(umap_1 > lim_x[1] & umap_1 < lim_x[2])  %>% filter(umap_2 > lim_y[1] & umap_2 < lim_y[2])-> early_eryth


eryth.centr = centroid(early_eryth)

start = Sys.time()
dist_cent = rbindlist(foreach(cell = rownames(subset(eryth.monocle, subset = (celltype.l2 == "Early Eryth"))[['umap']]@cell.embeddings[,1:2])) %dopar% {
	data.table(cell, as.numeric(dist(cbind(as.vector(eryth.centr), as.vector(subset(eryth.monocle, subset = (celltype.l2 == "Early Eryth"))[['umap']]@cell.embeddings[cell,1:2])))))
})
print(Sys.time() - start)
colnames(dist_cent) = c("cell", "dist")
dist_cent %>% arrange(desc(abs(dist))) %>% head(40) %>% select(cell) -> init_cells

##------------------------------------------Trajectory inference------------------------------------------
cds <- as.cell_data_set(x = eryth.monocle, assay = "refAssay", reductions = "umap", )
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")#, k=5)
cds <- learn_graph(cds)

##------------------------------------------Pseudotime Calculation-----------------------------------------
cds <- order_cells(cds, root_cells = init_cells$cell)

cds@principal_graph_aux@listData$UMAP$pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime/max(cds@principal_graph_aux@listData$UMAP$pseudotime)
d2n$DC.time <-   cds@principal_graph_aux@listData$UMAP$pseudotime
saveRDS(d2n, paste0(processed_data_dir, "/d2n_SeuratObj_PostQC_dCas9_NewMD_PT_BM.RDS"))


##------------------------------------------Panel C--------------------------------------------------------
p2 = plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
	   label_leaves=FALSE,
	   label_roots = FALSE,
           label_branch_points=FALSE,
           graph_label_size=0.0) + 
	theme(#legend.position = c(0.25, 0.15), legend.direction = "horizontal",
		legend.text = element_text(size = 3),
		legend.key.size = unit(2.3, 'mm'),
		legend.background = element_blank(),
		legend.title = element_text(size=5)) + 
	ggtitle("Pseudotime") +
        xlab("") +
        ylab("")

p3 = plot_cells(cds,
           color_cells_by = "celltype.l2",
	   show_trajectory_graph = F,
           label_cell_groups=FALSE, 
           label_leaves=TRUE, 
           label_branch_points=TRUE, 
           graph_label_size=1) + 
	guides(color = guide_legend(override.aes = list(size=1))) +
	theme(#legend.position = c(0.5, 0.1), legend.direction = "horizontal",
		legend.title = element_blank(),
		legend.background = element_blank(),
		legend.text = element_text(size = 4),
                legend.key.size = unit(1, 'mm')) +
	ggtitle("Cell type") +
        xlab("") +
        ylab("DC_2")

p4 <- DimPlot(d2n, sizes.highlight = .01, raster = F, pt.size = 0.5,  reduction = 'DC', cells.highlight = init_cells$cell) + 
	ggtitle("Initial Cells") +  
	theme(plot.title = element_text(size = 7, face = "plain"),
		axis.line=element_line(size=0.25)) + 
	NoLegend() +
	xlab("DC_1") +
        ylab("")



svg(paste0(plots_dir, "/Panel_C_BM_dataset.svg"))
p3 + p4  + p2 + plot_layout(ncol = 3, guides='collect')&theme(axis.text = element_text(size = 4),
                                                                axis.title = element_text(size = 5),
                                                                plot.title = element_text(size = 7, face = "bold"))
dev.off()





obj.merge.query_norm = NormalizeData(obj.merge.query)

#------------------------------------------Panel B-----------------------------------------------
p2 <- FeaturePlot(obj.merge.query_norm, features = "GFI1B", reduction = 'DC', cells = WhichCells(obj.merge.query_norm, expression = gene == "GFI1B"),  pt.size = 0.01, raster = F) + 
		xlim(c(-1.2, 1.0)) +
		ylim(c(0.75, 1.5)) + xlab("") +
                ylab("DC_2") + NoLegend()
p3 <- FeaturePlot(obj.merge.query_norm, features = "NFE2", reduction = 'DC', cells = WhichCells(obj.merge.query_norm, expression = gene == "NFE2"), pt.size = 0.01, raster = F) +
                xlim(c(-1.2, 1.0)) +
                ylim(c(0.75, 1.5)) + theme(legend.position = c(0.2, 0.95), legend.direction="horizontal", legend.text = element_text(size = 3), legend.key.size = unit(2.5, "mm"), legend.spacing.y = unit(0.1, 'mm')) +
		xlab("DC_1") +
                ylab("")
p4 <- FeaturePlot(obj.merge.query_norm,  features = "MYB", reduction = 'DC', cells = WhichCells(obj.merge.query_norm, expression = gene == "MYB"),  pt.size = 0.01, raster = F) +
                xlim(c(-1.2, 1.0)) +
                ylim(c(0.75, 1.5)) + 
		xlab("") +
                ylab("") + NoLegend()
 

svg(paste0(plots_dir, "/Panel_B_BM_dataset.svg"))
(p2|p3|p4) + plot_layout(ncol = 3)&theme(axis.text = element_text(size = 4),
                                                                axis.title = element_text(size = 5),
                                                                plot.title = element_text(size = 7))
dev.off() 




##------------------------------------------Time dist by dosage gene-----------------------------------
d2n_logFC = readRDS(paste0(data_dir, "/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS"))

gene.set = c("GFI1B", "NFE2", "MYB")

d2n_norm = NormalizeData(d2n)
time.list <- rbindlist(lapply(gene.set, function(x) {
  return(data.table(gene = x, time = d2n_norm$DC.time[ which(d2n_norm$gene == x ) ], expr = d2n_norm@assays$refAssay@data[, which(d2n_norm$gene == x ) ][x,]))
})) 


time.list_pseudo <- rbindlist(lapply(unique(d2n_norm$guide_crispr), function(x) {
  guide = unlist(str_split(x, pattern = "-"))[1]
  g = unlist(str_split(guide, pattern = "_"))[1]  
  c_l = unlist(str_split(x, pattern = "-"))[2]
  out = data.table(gene = g, time = mean(d2n_norm$DC.time[ which(d2n_norm$guide_crispr == x) ], na.rm = T), expr = d2n_logFC %>% 
															filter((guide_1 == guide) & (cell_line == c_l)) %>% 
															select(dosage_gene_log2FC) %>% 
															distinct())
  colnames(out) = c("gene", "time", "expr")
  out
}))  


time.list$gene = factor(time.list$gene, levels = c("GFI1B", "NFE2", "MYB"))


##------------------------------------------Panel D-------------------------------------------------
p = ggplot(time.list, aes(x = expr, y = time)) +
	geom_point(alpha=0.5, size=1) +
	geom_smooth(data = time.list %>% filter(!expr == 0.0),
			method=lm, alpha=0.5) +
	facet_wrap(~ gene, scale = "free_x") +
	theme_bw() + 
	xlab("Log Normalized UMI Expression") +
	ylab("Pseudotime") +
	theme(axis.text = element_text(size = 4),
        	axis.title = element_text(size = 5),
		strip.text = element_text(size = 5),
                plot.title = element_text(size = 7, face = "bold"))


svg(paste0(plots_dir, "/Panel_D_BM_dataset.svg"))
p
dev.off()

##------------------------------------------Panel E------------------------------------------
time.list_pseudo %>% filter(gene %in% gene.set) -> time.list_pseudo
p = ggplot(time.list_pseudo, aes(y = time, x = expr)) +
        geom_point(size = 1, alpha = 0.75) +
        geom_smooth(method=loess, alpha = 0.25, size = 0.75) +
        facet_wrap(~ gene, scale = "free_x") +
        theme_bw() +
	xlab("Log2 FC") +
        ylab("Pseudotime") +
        theme(axis.text = element_text(size = 4),
                axis.title = element_text(size = 5),
		strip.text = element_text(size = 5),
                plot.title = element_text(size = 7, face = "bold"))


svg(paste0(plots_dir, "/Panel_E_BM_dataset.svg"))
p
dev.off()

##------------------------------------Gene PT correlation---------------------------------------
start = Sys.time()
d2n_corr_df = rbindlist(foreach(gene = rownames(d2n)) %dopar% {
	c = cor.test(d2n$DC.time, as.vector(d2n[gene, ]@assays$refAssay@data), method = "spearman")$estimate
	p = cor.test(d2n$DC.time, as.vector(d2n[gene, ]@assays$refAssay@data), method = "spearman")$p.value
	data.table(gene = gene, corr = c, pval = p)
})
Sys.time() - start

d2n_corr_df %>% 
	drop_na() %>% 
	filter(pval < 0.01) %>% 
	arrange(desc(corr)) -> d2n_corr_df

d2n_corr_df$gene = factor(d2n_corr_df$gene, levels = d2n_corr_df$gene)
d2n_corr_df$pval = p.adjust(d2n_corr_df$pval, method = "bonferroni")

p = ggplot(d2n_corr_df, aes(x=as.factor(gene), y=corr)) +
  geom_segment( aes(x=as.factor(gene), xend=as.factor(gene), y=0, yend=corr, alpha = -log10(pval) > 2), color="grey") +
  geom_point( color="orange", aes(alpha=-log10(pval) > 2)) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=90, size=6)
  ) +
  guides(size = "none") +
  ylab("Spearman Correlation") +
  xlab("Gene")

svg(paste0(plots_dir, "/BM_dataset_integr_gene_corr_with_PT.svg"), width = 118.8, height = 59.4,  res=720, units= "mm")
p
dev.off()

# Top 10 corr genes with PT
d2n_corr_df %>% arrange(desc(abs(corr))) %>% head(10) %>% select(gene) -> corr_genes 


time.list_pseudo <- rbindlist(lapply(unique(d2n_norm$guide_crispr), function(x) {
  guide = unlist(str_split(x, pattern = "-"))[1]
  c_l = unlist(str_split(x, pattern = "-"))[2]
  out = data.table(d2n_logFC %>%
                              filter((guide_1 == guide) & (cell_line == c_l) & (gene %in% corr_genes$gene)) %>%
                              select(gene, avg_log2FC) %>%
                              distinct(),
		    time = mean(d2n_norm$DC.time[ which(d2n_norm$guide_crispr == x) ], na.rm = T))
  colnames(out) = c("gene", "expr", "time")
  out
}))

p = ggplot(time.list_pseudo, aes(x = time, y = expr)) +
        geom_point() +
        #geom_smooth(method=glm) +
        facet_wrap(~ gene) +
        theme_bw()


svg(paste0(plots_dir, "/BM_dataset_integr_trans_gene_vs_time_guide_agg.svg"), width = 118.8, height = 59.4,  res=720, units= "mm")
p
dev.off()

