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

registerDoParallel(cores=20)


##------------------------------------------Setting paths---------------------------------------
folder_name = basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)




##------------------------------------------Preparation of reference-----------------------------
#reference = readRDS(paste0(data_dir, "/seurat_main_cohort.rds"))
#cells = unique(reference@meta.data$ct_simple[!str_detect(reference@meta.data$ct_simple, "B ") & !str_detect(reference@meta.data$ct_simple, "T ")])
#reference <- subset(reference, subset = ct_simple %in% cells) 
#alldata.list <- SplitObject(reference, split.by = "cohort")
#for (i in 1:length(alldata.list)) {
#    alldata.list[[i]] <- NormalizeData(alldata.list[[i]], verbose = FALSE)
#    alldata.list[[i]] <- FindVariableFeatures(alldata.list[[i]], selection.method = "vst",
#        nfeatures = 5000, verbose = FALSE)
#}
#alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:30,
#    reduction = "cca", anchor.features = 5000)
#reference <- IntegrateData(anchorset = alldata.anchors, dims = 1:30, new.assay.name = "CCA")
#reference = ScaleData(reference, verbose = FALSE)
#reference <- RunPCA(reference, npcs = 30, verbose = FALSE)
#reference <- RunUMAP(reference, dims = 1:30)
#saveRDS(reference, paste0(processed_data_dir, "/aml_reference_integrated_cca.rds"))



##------------------------------------------AML dataset-----------------------------------------

# Load reference dataset
reference =  readRDS(paste0(processed_data_dir, "/aml_reference_integrated_cca.rds"))
cells = unique(reference@meta.data$ct_simple[!str_detect(reference@meta.data$ct_simple, "B ") & !str_detect(reference@meta.data$ct_simple, "T ")])
reference <- subset(reference, subset = ct_simple %in% cells)
DefaultAssay(reference) = "CCA"


##------------------------------------------Load and preprocess D2N dataset---------------------
d2n_full <- readRDS(file = paste0(data_dir, "/d2n_SeuratObj_PostQC_dCas9_NewMD.RDS"))
counts <- GetAssayData(d2n_full, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('CRISPRi', 'CRISPRa'))),]
d2n <- subset(d2n_full, features = rownames(counts))

d2n = NormalizeData(d2n)
VariableFeatures(d2n) = rownames(d2n[["RNA"]])
d2n[["CCA"]] = d2n[["RNA"]]
reference = RunPCA(reference)
reference <- RunUMAP(reference, dims = 1:30, return.model = T)

##------------------------------------------Integration------------------------------------------
anchors <- FindTransferAnchors(
  reference = reference,
  query = d2n,
  k.filter = NA,
  reference.reduction = "pca",
  reference.assay = "CCA",
  query.assay = "CCA",
  features = intersect(rownames(x = reference), rownames(x = d2n)),
  dims = 1:50)

d2n <- MapQuery(
  anchorset = anchors,
  query = d2n,
  reference = reference,
  refdata = list(
    ct = "ct"),
  reduction.model = "umap")

###------------------------------------------UMAP for shared embeddings---------------------------
d2n$ct <- d2n$predicted.ct
d2n$group <- 'query'
reference$group <- 'reference'

emb.merge <- merge( reference[['pca']], d2n[['ref.pca']])
emb.merge@assay.used <- 'CCA'
DefaultAssay(reference) <- "CCA"
DefaultAssay(d2n) <- 'CCA'


obj.merge <- merge(DietSeurat(reference, assays = 'CCA'), DietSeurat(d2n, assays = 'CCA') )
obj.merge[['pca']] <- emb.merge
obj.merge <- RunUMAP(obj.merge, dims = 1:50, reduction="pca")



#------------------------------------------Diffusion map analysis-----------------------------------
df.map <- DiffusionMap(data = obj.merge[['pca']]@cell.embeddings[,1:30])
obj.merge[['DC']] <- CreateDimReducObject(embeddings = df.map@eigenvectors*1000, key = 'DC_', assay = 'CCA')

d2n_types = unique(d2n$ct)
ref_types = setdiff(unique(reference$ct), d2n_types)
obj.merge$ct = factor(obj.merge$ct, levels=c(d2n_types, ref_types))
cols = DiscretePalette(25)
 
#------------------------------------------Panel A---------------------------------------------------
p1 <- DimPlot(obj.merge, cols = alpha(cols, 0.33), group.by = 'ct',  pt.size = 0.01, raster = F) + NoLegend() + ggtitle("All")
p2 <- DimPlot(obj.merge,  cols = alpha(cols, 0.33), group.by = 'ct', cells = Cells(d2n),  pt.size = 0.01, raster = F) + NoLegend() + ggtitle("Perturbed")
p3 <- DimPlot(obj.merge,   cols = alpha(cols, 0.33), group.by = 'ct', cells = Cells(reference),  pt.size = 0.01, raster = F) + NoLegend() + ggtitle("Reference")
p4 =  DimPlot(obj.merge,  cols = alpha(cols, 0.33), group.by = 'ct', cells = Cells(d2n), raster = F) + theme(legend.position= "bottom",
																legend.text = element_text(size = 9),
																legend.key.size = unit(2, 'cm'))

leg = gtable_filter(ggplot_gtable(ggplot_build(p4)), "guide-box")

png(paste0(plots_dir, "/Panel_A_AML_dataset.png"), width = 9, height = 4.5,  res=720, units= "in")
(p1|p2|p3)/leg + plot_layout(nrow = 2, height = c(4, 1))&theme(axis.text = element_text(size = 7),
                                                                axis.title = element_text(size = 9),
                                                                plot.title = element_text(size = 12))
dev.off()

 
obj.merge@meta.data$cell_name = rownames(obj.merge@meta.data)
 
obj.merge.query <- subset(obj.merge, cells = Cells(d2n))
d2n = subset(obj.merge, cells = Cells(obj.merge.query))
d2n[['DC']] <- obj.merge.query[['DC']]
obj.merge.query_norm = NormalizeData(obj.merge.query)

##------------------------------------------Pseudotime analysis---------------------------------------

###-----------------------------------------Erythroid Cells-------------------------------------------

eryth.monocle = subset(d2n, subset = (ct == "Early erythroid progenitor") | (ct == "Late erythroid progenitor"))
eryth.monocle[['umap']] = NULL
eryth.monocle[['umap']] <- eryth.monocle[['DC']]

eryth.monocle[['umap']]@cell.embeddings <- eryth.monocle[['umap']]@cell.embeddings[,1:2]


##------------------------------------------Select root cells------------------------------------------
early_eryth = as.data.table(subset(eryth.monocle, subset = (ct== "Early erythroid progenitor"))[['umap']]@cell.embeddings[,1:2])
early_eryth$cell = colnames(subset(eryth.monocle, subset = (ct== "Early erythroid progenitor")))
lim_y = quantile(early_eryth$umap_2, probs = c(0.9, 1.0))
lim_x =  quantile(early_eryth$umap_1, probs = c(0.3, 0.6))

early_eryth %>% filter(umap_1 > lim_x[1] & umap_1 < lim_x[2])  %>% filter(umap_2 > lim_y[1] & umap_2 < lim_y[2])-> early_eryth


eryth.centr = centroid(early_eryth[, 1:2])

start = Sys.time()
dist_cent = rbindlist(foreach(cell = early_eryth$cell) %dopar% {
        data.table(cell, 
			as.numeric(dist(cbind(as.vector(eryth.centr), 
			as.vector(subset(eryth.monocle, 
						subset = (ct == "Early erythroid progenitor"))[['umap']]@cell.embeddings[cell,1:2])))))
})
print(Sys.time() - start)
colnames(dist_cent) = c("cell", "dist")
dist_cent %>% arrange(desc(abs(dist))) %>% head(40) %>% select(cell) -> init_cells

##------------------------------------------Trajectory inference------------------------------------------
cds <- as.cell_data_set(x = eryth.monocle, assay = "CCA", reductions = "umap")
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")#, k=5)
cds <- learn_graph(cds)


##------------------------------------------Pseudotime Calculation-----------------------------------------
cds <- order_cells(cds, root_cells = init_cells$cell)
cds@principal_graph_aux@listData$UMAP$pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime/max(cds@principal_graph_aux@listData$UMAP$pseudotime)
d2n$DC.time <-   cds@principal_graph_aux@listData$UMAP$pseudotime
saveRDS(d2n, paste0(processed_data_dir, "/d2n_SeuratObj_PostQC_dCas9_NewMD_PT_AML.RDS"))

##------------------------------------------Panel C--------------------------------------------------------
 
p2 = plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) +
        theme(legend.position = "bottom") +
	ggtitle("Pseudotime") + 
	xlab("DC_1") +
	ylab("DC_2")

p3 = plot_cells(cds,
           color_cells_by = "ct",
           show_trajectory_graph = F,
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5) +
        theme(legend.position = "bottom",
                legend.title = element_blank()) +
        ggtitle("Cell type") +
        xlab("DC_1") +
        ylab("DC_2")

p4 <- DimPlot(d2n, sizes.highlight = .01, raster = F, pt.size = 0.5,  reduction = 'DC', cells.highlight = init_cells$cell) +
        ggtitle("Initial Cells") +
        theme(plot.title = element_text(size = 12, face = "plain"),
                axis.line=element_line(size=0.25)) +
        NoLegend()



png(paste0(plots_dir, "/Panel_C_AML_dataset.png"), width = 9, height = 4,  res=600, units= "in")
p3 + p4  + p2 + plot_layout(ncol = 3)&theme(axis.text = element_text(size = 7),
                                                                axis.title = element_text(size = 9),
                                                                plot.title = element_text(size = 12, face = "bold"))
dev.off()

#------------------------------------------Panel B-----------------------------------------------
p2 <- FeaturePlot(obj.merge.query_norm, 
			features = "GFI1B", 
			reduction = 'DC', 
			cells = WhichCells(obj.merge.query_norm, expression = gene == "GFI1B"),  
			pt.size = 0.01, 
			raster = F) +
	NoLegend()
p3 <- FeaturePlot(obj.merge.query_norm, 
			features = "NFE2", 
			reduction = 'DC', 
			cells = WhichCells(obj.merge.query_norm, expression = gene == "NFE2"), 
			pt.size = 0.01, 
			raster = F) +
	theme(legend.text = element_text(size = 7))

p4 <- FeaturePlot(obj.merge.query_norm,  
			features = "MYB", 
			reduction = 'DC', 
			cells = WhichCells(obj.merge.query_norm, expression = gene == "MYB"),  
			pt.size = 0.01, 
			raster = F) +
	NoLegend()
 
png(paste0(plots_dir, "/Panel_b_AML_dataset.png"), width = 9, height = 4,  res=600, units= "in")
p2+p3+p4 + plot_layout(guides = "collect", ncol=3)&theme(axis.text = element_text(size = 7),
                                                                axis.title = element_text(size = 9),
                                                                plot.title = element_text(size = 12))
dev.off()



##------------------------------------------Time dist by dosage gene-----------------------------------
d2n_logFC = readRDS(paste0(data_dir, "/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS"))

d2n_true_erythro = subset(d2n, subset = DC_1 < eryth.centr[1])
gene.set = c("GFI1B", "NFE2", "MYB")

d2n_true_erythro_norm = NormalizeData(d2n_true_erythro)
time.list <- rbindlist(lapply(gene.set, function(x) {
  return(data.table(gene = x, time = d2n_true_erythro_norm$DC.time[ which(d2n_true_erythro_norm$gene == x ) ], expr = d2n_true_erythro_norm@assays$CCA@data[, which(d2n_true_erythro_norm$gene == x ) ][x,]))
}))

time.list_pseudo <- rbindlist(lapply(unique(d2n_true_erythro_norm$guide_crispr), function(x) {
  guide = unlist(str_split(x, pattern = "-"))[1]
  g = unlist(str_split(guide, pattern = "_"))[1]
  c_l = unlist(str_split(x, pattern = "-"))[2]
  out = data.table(gene = g, time = mean(d2n_true_erythro_norm$DC.time[ which(d2n_true_erythro_norm$guide_crispr == x) ], na.rm = T), expr = d2n_logFC %>%
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
        theme(axis.text = element_text(size = 7),
                axis.title = element_text(size = 9),
                plot.title = element_text(size = 12, face = "bold"))
 
 
png(paste0(plots_dir, "/Panel_D_AML_dataset.png"), width = 9, height = 3.33,  res=600, units= "in")
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
        theme(axis.text = element_text(size = 7),
                axis.title = element_text(size = 9),
                plot.title = element_text(size = 12, face = "bold"))
 

png(paste0(plots_dir, "/Panel_E_AML_dataset.png"), width = 9, height = 3.33,  res=600, units= "in")
p
dev.off()

##------------------------------------Gene PT correlation---------------------------------------
start = Sys.time()
d2n_true_erythro_corr_df = rbindlist(foreach(gene = rownames(d2n_true_erythro)) %dopar% {
        c = cor.test(d2n_true_erythro$DC.time, as.vector(d2n_true_erythro[gene, ]@assays$CCA@data), method = "spearman")$estimate
        p = cor.test(d2n_true_erythro$DC.time, as.vector(d2n_true_erythro[gene, ]@assays$CCA@data), method = "spearman")$p.value
        data.table(gene = gene, corr = c, pval = p)
})
Sys.time() - start

d2n_true_erythro_corr_df %>% 
	drop_na() %>% 
	filter(pval < 0.01) %>% 
	arrange(desc(corr)) -> d2n_true_erythro_corr_df

d2n_true_erythro_corr_df$gene = factor(d2n_true_erythro_corr_df$gene, levels = d2n_true_erythro_corr_df$gene)
d2n_true_erythro_corr_df$pval = p.adjust(d2n_true_erythro_corr_df$pval, method = "bonferroni")

p = ggplot(d2n_true_erythro_corr_df, aes(x=as.factor(gene), y=corr)) +
  geom_segment( aes(x=as.factor(gene), xend=as.factor(gene), y=0, yend=corr), color="grey") +
  geom_point( color="orange") +
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

png(paste0(plots_dir, "/AML_dataset_integr_gene_corr_with_PT.png"), width = 9, height = 6,  res=600, units= "in")
p
dev.off()
 
# Top 10 corr genes with PT
d2n_true_erythro_corr_df %>% arrange(desc(abs(corr))) %>% head(10) %>% select(gene) -> corr_genes

time.list_pseudo <- rbindlist(lapply(unique(d2n_true_erythro_norm$guide_crispr), function(x) {
  guide = unlist(str_split(x, pattern = "-"))[1]
  c_l = unlist(str_split(x, pattern = "-"))[2]
  out = data.table(d2n_logFC %>%
                              filter((guide_1 == guide) & (cell_line == c_l) & (gene %in% corr_genes$gene)) %>%
                              select(gene, avg_log2FC) %>%
                              distinct(),
                    time = mean(d2n_true_erythro_norm$DC.time[ which(d2n_true_erythro_norm$guide_crispr == x) ], na.rm = T))
  colnames(out) = c("gene", "expr", "time")
  out
}))
 
p = ggplot(time.list_pseudo, aes(x = time, y = expr)) +
        geom_point() +
        #geom_smooth(method=glm) +
        facet_wrap(~ gene) +
        theme_bw()


png(paste0(plots_dir, "/AML_dataset_integr_trans_gene_vs_time_guide_agg.png"), width = 9, height = 6,  res=600, units= "in")
p
dev.off()

