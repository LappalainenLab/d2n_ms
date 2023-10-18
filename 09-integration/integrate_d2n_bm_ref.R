#specify monocle lib path
.libPaths( c( "/proj/lappalainen_lab1/users/marii/chip_seq_ann" , .libPaths() ) )


library(stringr)
library(monocle3)
library(Seurat)
library(corrplot)
library(cowplot)
library(grid)
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
if (T){
reference =  readRDS(paste0(data_dir, "/bm_ref.rds"))
reference <- subset(reference, subset = celltype.l1 %in% c('DC', 'HSPC', 'Mono', 'NK'))
VariableFeatures(reference) <- rownames(reference@reductions$refDR@feature.loadings)

##------------------------------------------Load and preprocess D2N dataset---------------------
d2n_full <- readRDS(file = paste0(data_dir, "/d2n_SeuratObj_PostQC_dCas9_NewMD.RDS"))
counts <- GetAssayData(d2n_full, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('CRISPRi', 'CRISPRa'))),]#, 'MYB', "GFI1B", 'NFE2'))),]
d2n <- subset(d2n_full, features = rownames(counts))

VariableFeatures(d2n) = rownames(d2n[["RNA"]])

d2n@meta.data["MYB"] = t(d2n@assays$RNA["MYB"])
d2n@meta.data["NFE2"] = t(d2n@assays$RNA["NFE2"])
d2n@meta.data["GFI1B"] = t(d2n@assays$RNA["GFI1B"])


##------------------------------------------Integration------------------------------------------
genes = intersect(rownames(x = reference), VariableFeatures(object = d2n))
genes = genes[!genes %in% c("MYB", "NFE2", "GFI1B")]

## Difussion map
df.map <- DiffusionMap(data = reference[["refDR"]]@cell.embeddings)
reference[['DC']] <- CreateDimReducObject(embeddings = df.map@eigenvectors*1000, key = 'DC_', assay = 'refAssay')
reference@meta.data$cell_name = rownames(reference@meta.data)

eryth.monocle = subset(reference, subset = celltype.l2 %in% c("Early Eryth", "Late Eryth", "HSC", "EMP"), features = rownames(reference))
#eryth.monocle[['umap']] = NULL
eryth.monocle[['umap']] <- eryth.monocle[['DC']]
eryth.monocle[['umap']]@cell.embeddings <- eryth.monocle[['umap']]@cell.embeddings[,1:2]

early_eryth = subset(eryth.monocle, subset = celltype.l2 == "HSC")[['umap']]@cell.embeddings[,1:2]
 
 
eryth.centr = centroid(early_eryth) 
 
start = Sys.time() 
dist_cent = rbindlist(foreach(cell = rownames(subset(eryth.monocle, subset = (celltype.l2 == "HSC"))[['umap']]@cell.embeddings[,1:2])) %dopar% {
        data.table(cell, as.numeric(dist(cbind(as.vector(eryth.centr), as.vector(subset(eryth.monocle, subset = (celltype.l2 == "HSC"))[['umap']]@cell.embeddings[cell,1:2])))))
})

colnames(dist_cent) = c("cell", "dist") 
dist_cent %>% arrange(abs(dist)) %>% head(40) %>% select(cell) -> init_cells

##------------------------------------------Trajectory inference------------------------------------------
cds <- as.cell_data_set(x = eryth.monocle, assay = "refAssay", reductions = "umap", )
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")#, cluster_method = 'louvain') #, k=5)
cds <- learn_graph(cds)

##------------------------------------------Pseudotime Calculation-----------------------------------------
cds <- order_cells(cds, root_cells = init_cells$cell)

cds@principal_graph_aux@listData$UMAP$pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime/max(cds@principal_graph_aux@listData$UMAP$pseudotime)
eryth.monocle@meta.data$DC.time <-   cds@principal_graph_aux@listData$UMAP$pseudotime
eryth.monocle@meta.data$label = paste0(eryth.monocle@meta.data$celltype.l2, "_", eryth.monocle@meta.data$group)
 
left_join(reference@meta.data, eryth.monocle@meta.data) -> reference@meta.data
rownames(reference@meta.data) <- colnames(reference@assays$refAssay)

anchors <- FindTransferAnchors(
  reference = reference, 
  query = d2n, 
  k.filter = 10, 
  reduction='pcaproject',
  reference.assay = "refAssay",
  query.assay = "RNA",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = genes,
  dims = 1:50)

d2n <- MapQuery(
  anchorset = anchors,
  query = d2n,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    DC.time = "DC.time"),
  reference.reduction = "refDR", 
  reduction.model = "refUMAP")


###------------------------------------------UMAP for shared embeddings---------------------------
d2n$celltype.l2 <- d2n$predicted.celltype.l2
d2n$celltype.l1 <- d2n$predicted.celltype.l1
d2n$DC.time = d2n$predicted.DC.time
 
d2n$group <- 'query'
reference$group <- 'reference'

emb.merge <- merge( reference[['refDR']], d2n[['ref.refDR']])
emb.merge@assay.used <- 'refAssay'
DefaultAssay(reference) <- 'refAssay'
DefaultAssay(d2n) = "RNA"

obj.merge <- merge(DietSeurat(reference, assays = 'refAssay'), DietSeurat(d2n, assays = 'RNA') )
obj.merge[['refDR']] <- emb.merge
obj.merge <- RunUMAP(obj.merge, dims = 1:30, reduction="refDR")
d2n_types = unique(d2n$celltype.l2)
ref_types = setdiff(unique(reference$celltype.l2), d2n_types)
obj.merge$celltype.l2 = factor(obj.merge$celltype.l2, levels=c(d2n_types, ref_types))

scores = MappingScore(
  anchors=anchors,
  verbose = TRUE,
  ndim=50
)

saveRDS(obj.merge, paste0(processed_data_dir, "/integrated_SeuratObject_BM_d2n.RDS"))
}


cols = DiscretePalette(24)
Idents(obj.merge) <- obj.merge@meta.data$group
to_plot_obj.merge = subset(obj.merge, downsample=5000)

#------------------------------------------Panel A---------------------------------------------------
p1 <- DimPlot(to_plot_obj.merge, cols = alpha(cols, 0.33), group.by = 'celltype.l2',  pt.size = 0.01, raster = F) + 
	NoLegend() + 
	ggtitle("All") +
	xlab("") + 
        ylab("UMAP_2")
p2 <- DimPlot(to_plot_obj.merge,  cols = alpha(cols, 0.33), group.by = 'celltype.l2', cells = WhichCells(to_plot_obj.merge, expression =  group == "query"),  pt.size = 0.01, raster = F) + 
	NoLegend() + 
	ggtitle("Perturbed") + 
	xlab("UMAP_1") +
	ylab("")

p3 <- DimPlot(to_plot_obj.merge,   cols = alpha(cols, 0.33), group.by = 'celltype.l2', cells = WhichCells(to_plot_obj.merge, expression =  group == "reference"),  pt.size = 0.01, raster = F) + 
	NoLegend() + 
	ggtitle("Reference") +
	xlab("") + 
        ylab("")

p4 =  DimPlot(to_plot_obj.merge,  cols = alpha(cols, 0.33), group.by = 'celltype.l2', raster = F) + 
	guides(colour = guide_legend(override.aes = list(size=2), nrow=3)) +
	theme(legend.position= "bottom",
        	legend.text = element_text(size = 9),
		legend.key.size = unit(7, 'pt'),
		legend.spacing.x = unit(6, 'pt'),
		legend.spacing.y = unit(3, 'pt'))

leg = cowplot::get_legend(p4)
svg(paste0(plots_dir, "/Panel_A_BM_dataset_new.svg"), width=8.3, height=3.1)
p1 + p2 + p3 + plot_layout(ncol = 3)&theme(axis.text = element_text(size = 9),
                                           axis.title = element_text(size = 11),
                                           plot.title = element_text(size = 12))
dev.off()

svg(paste0(plots_dir, "/Panel_A_BM_dataset_legend.svg"))
grid.newpage()
grid.draw(leg)
dev.off()

d2n = subset(obj.merge, subset = group == "query")

to_plot_cds = subset(cds, downsample=5000)
to_plot_eryth.monocle = subset(eryth.monocle, downsample=5000)
##------------------------------------------Panel B--------------------------------------------------------
p2 = plot_cells(to_plot_cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) +
        xlab("") +
        ylab("DC_2") +
        labs(color = "Pseudotime") +
        theme(legend.position = c(0.15, 0.9),
               legend.title = element_text(size = 9),
               legend.direction="horizontal",
               legend.text = element_text(size = 9),
               legend.key.size = unit(7, "pt"),
               legend.spacing.y = unit(3, 'pt'))
 

p3 = DimPlot(to_plot_eryth.monocle, sizes.highlight = .01, label=F, raster = F, pt.size = 0.5,  reduction = 'DC', group.by = "celltype.l2") +
        theme(plot.title = element_text(size = 12, face = "plain"),
                axis.line=element_line(size=0.25),
                legend.position = c(0.15, 0.9),
                #legend.direction = "horizontal",
                legend.text = element_text(size = 9),
                legend.key.size = unit(7, 'pt'),
                legend.background = element_blank(),
                legend.title = element_text(size=9)) +
        xlab("") +
        ggtitle("") +
        ylab("DC_2")

p4 <- DimPlot(to_plot_eryth.monocle, sizes.highlight = .01, raster = F, pt.size = 0.5,  reduction = 'DC', cells.highlight = init_cells$cell) +
        theme(plot.title = element_text(size = 12, face = "plain"),
                axis.line=element_line(size=0.25),
                legend.position = c(0.15, 0.9),
                #legend.direction = "horizontal",
                legend.text = element_text(size = 9),
                legend.key.size = unit(7, 'pt'),
                legend.background = element_blank(),
                legend.title = element_text(size=9)) +
        scale_color_manual(name = "", labels = c("other", "centroid"), values = c("gray", 'red')) +
        xlab("DC_1") +
        ylab("")


svg(paste0(plots_dir, "/Panel_B_BM_dataset_new.svg"), width = 8.3, height=3.1)
p3 + p4  + p2 + plot_layout(ncol = 3)&theme(axis.text = element_text(size = 9),
                                                                axis.title = element_text(size = 11),
                                                                plot.title = element_text(size = 12, face = "bold"))
dev.off()

to_plot_d2n = subset(d2n, downsample=5000)

#------------------------------------------Panel C-----------------------------------------------
p2 <- FeaturePlot(to_plot_d2n, reduction = 'refDR', 
			cells = WhichCells(d2n, expression = gene == "GFI1B"),  
			pt.size = 0.01, raster = F,  features = "GFI1B") + 
		xlab("") +
                ylab("DC_2") +
		labs(color = "Norm UMI") + 
		theme(legend.position = c(0.4, 0.9), 
			legend.title = element_text(size = 9), 
			legend.direction="horizontal", 
			legend.text = element_text(size = 9), 
			legend.key.size = unit(7, "pt"), 
			legend.spacing.y = unit(3, 'pt'))
p3 <- FeaturePlot(to_plot_d2n, features = "NFE2", reduction = 'refDR', cells = WhichCells(d2n, expression = gene == "NFE2"), pt.size = 0.01, raster = F) +
		labs(color = "Norm UMI") +
		xlab("DC_1") +
                ylab("") +
                theme(legend.position = c(0.4, 0.9),
                        legend.title = element_text(size = 9),
                        legend.direction="horizontal",
                        legend.text = element_text(size = 9),
                        legend.key.size = unit(7, "pt"),
                        legend.spacing.y = unit(3, 'pt'))
p4 <- FeaturePlot(to_plot_d2n,  features = "MYB", reduction = 'refDR', cells = WhichCells(d2n, expression = gene == "MYB"),  pt.size = 0.01, raster = F) +
		labs(color = "Norm UMI") +
		xlab("") +
                ylab("") +
                theme(legend.position = c(0.4, 0.9),
                        legend.title = element_text(size = 9),
                        legend.direction="horizontal",
                        legend.text = element_text(size = 9),
                        legend.key.size = unit(7, "pt"),
                        legend.spacing.y = unit(3, 'pt'))
 

svg(paste0(plots_dir, "/Panel_C_BM_dataset_new.svg"), width = 8.3, height=3.06)
(p2|p3|p4) + plot_layout(ncol = 3)&theme(axis.text = element_text(size = 9),
                                                                axis.title = element_text(size = 11),
                                                                plot.title = element_text(size = 12))
dev.off() 





##------------------------------------------Time dist by dosage gene-----------------------------------
d2n_logFC = readRDS(paste0(data_dir, "/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS"))

gene.set = c("GFI1B", "NFE2", "MYB")

d2n_norm = NormalizeData(d2n)
d2n$DC.time = as.numeric(d2n$DC.time)
time.list <- rbindlist(lapply(gene.set, function(x) {
  cells = which(d2n_norm$gene == x )
  return(data.table(gene = x, time = d2n$DC.time[cells], expr = d2n_norm@assays$RNA@data[,cells][x,]))
})) 


time.list_pseudo <- rbindlist(lapply(unique(d2n_norm$guide_crispr), function(x) {
  guide = unlist(str_split(x, pattern = "-"))[1]
  g = unlist(str_split(guide, pattern = "_"))[1]  
  c_l = unlist(str_split(x, pattern = "-"))[2]
  out = data.table(gene = g, time = mean(d2n$DC.time[ which(d2n$guide_crispr == x) ], na.rm = T), expr = d2n_logFC %>% 
															filter((guide_1 == guide) & (cell_line == c_l)) %>% 
															select(dosage_gene_log2FC) %>% 
															distinct())
  colnames(out) = c("gene", "time", "expr")
  out
}))  


time.list$gene = factor(time.list$gene, levels = c("GFI1B", "NFE2", "MYB"))
time.list$time  = as.numeric(time.list$time)

##------------------------------------------Panel D-------------------------------------------------
theme_set(theme_bw())
p = ggplot(time.list, aes(x = expr, y = time)) + 
        geom_hex() +
	stat_bin2d(bins = 50, show.legend = F) +
	scale_fill_gradient(low = "grey85", high = "grey30", guide = guide_colorbar(direction = "horizontal", barheight = 0.5))  +
        geom_smooth(data = time.list %>% filter(!expr == 0.0),
                        method=lm, alpha=0.5, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free_x") + 
        theme_bw() +  
        xlab("Cis gene normalized UMI expression") +
        ylab("Pseudotime") + 
        theme(legend.key = element_blank(),
		legend.position = "bottom",
                strip.background = element_rect(colour="white", fill="white"),
                axis.text = element_text(size = 9),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 9),
                plot.title = element_text(size = 12, face = "bold"))
 
 
svg(paste0(plots_dir, "/Panel_D_BM_dataset_hex_new.svg"), width = 8.3, height=3.85)
p
dev.off() 
 


##------------------------------------------Panel E------------------------------------------
time.list_pseudo %>% filter(gene %in% gene.set) -> time.list_pseudo
time.list_pseudo$gene = factor(time.list_pseudo$gene, levels = c("GFI1B", "NFE2", "MYB"))
p = ggplot(time.list_pseudo, aes(y = time, x = expr)) +
        geom_point(size = 1, alpha = 0.75) +
        geom_smooth(method=loess, alpha = 0.25, size = 0.75, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free_x") +
        theme_bw() +
	xlab("Cis gene log2(FC)") +
        ylab("Mean pseudotime") +
        theme(axis.text = element_text(size = 9),
		legend.key = element_blank(), 
		strip.background = element_rect(colour="white", fill="white"),
                axis.title = element_text(size = 11),
		strip.text = element_text(size = 9),
                plot.title = element_text(size = 12, face = "bold"))


svg(paste0(plots_dir, "/Panel_E_BM_dataset_new.svg"), width = 8.3, height=3.24)
p
dev.off()






#------------------------------------------NOT_USED----------------------------------------------------   
#------------------------------------------Pseudotime_analysis_for_Ref_dataset-------------------------
if (F){
#------------------------------------------Reference---------------------------------------------------
##------------------------------------------Time dist by dosage gene-----------------------------------
reference_true_eryth = subset(reference, subset = (ct == "Early Eryth") | (ct == "Late Eryth"))

gene.set = c("GFI1B", "NFE2", "MYB")

time.list <- rbindlist(lapply(gene.set, function(x) {
  return(data.table(gene = x, time = reference_true_eryth$DC.time, expr = reference_true_eryth@assays$RNA@data[x,]))
}))


time.list$gene = factor(time.list$gene, levels = c("GFI1B", "NFE2", "MYB"))
 
 

#Hex version
p = ggplot(time.list, aes(x = expr, y = time)) +
        geom_hex() +
        stat_bin2d(bins = 50, show.legend = F) +
        scale_fill_gradient(low = "grey85", high = "grey30", guide = guide_colorbar(direction = "horizontal", barheight = 0.5))  +
        geom_smooth(data = time.list %>% filter(!expr == 0.0),
                        method=lm, alpha=0.5, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free_x") +
        theme_bw() +
        xlab("Cis gene normalized UMI expression") +
        ylab("Pseudotime") +
        theme(legend.key = element_blank(),
                legend.position = "bottom",
                strip.background = element_rect(colour="white", fill="white"),
                axis.text = element_text(size = 9),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))

svg(paste0(plots_dir, "/Panel_D_BM_dataset_ref_hex.svg"), width = 8.3, height=3.85)
p
dev.off()

reference@meta.data$GFI1B_bin = paste0("GFI1B_", t(round(reference@assays$RNA["GFI1B"] / 0.1 + 1, 0)))
reference@meta.data$NFE2_bin = paste0("NFE2_", t(round(reference@assays$RNA["NFE2"] / 0.1 + 1, 0)))
reference@meta.data$MYB_bin = paste0("MYB_", t(round(reference@assays$RNA["MYB"] / 0.1 + 1, 0)))


time.list_pseudo <- rbindlist(lapply(gene.set, function(y){
        rbindlist(lapply(unique(unlist(reference@meta.data[paste0(y, "_bin")])), function(x) {
          bin = unlist(str_split(x, pattern = "_"))[2]
          gene = unlist(str_split(x, pattern = "_"))[1]
          cl = rownames(reference@meta.data[reference@meta.data[paste0(y, "_bin")] == x,])
          out = data.table(gene = gene,
                                time = mean(reference@meta.data[cl,]$DC.time, na.rm = T),
                                expr = mean(reference@assays$RNA[gene, cl], na.rm = T),
                                bin = bin)
          colnames(out) = c("gene", "time", "expr", "bin")
          out
        }))
}))


##------------------------------------------Panel E------------------------------------------
time.list_pseudo %>% filter(gene %in% gene.set) -> time.list_pseudo
time.list_pseudo$gene = factor(time.list_pseudo$gene, levels = c("GFI1B", "NFE2", "MYB"))

p = ggplot(time.list_pseudo, aes(y = time, x = expr)) +
        geom_point(size = 1, alpha = 0.75) +
        geom_smooth(method=loess, alpha = 0.25, size = 0.75, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free_x") +
        theme_bw() +
        xlab("Cis gene log2(FC)") +
        ylab("Mean pseudotime") +
        theme(axis.text = element_text(size = 9),
                legend.key = element_blank(),
                strip.background = element_rect(colour="white", fill="white"),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))

svg(paste0(plots_dir, "/Panel_E_BM_dataset_ref.svg"), width = 8.3, height=3.24)
p
dev.off()
}
