library(here)
here::i_am("integrate_d2n_aml_ref.R")

#specify monocle lib path
.libPaths( c( "/proj/lappalainen_lab1/users/marii/chip_seq_ann" , .libPaths() ) )

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
library(grid)
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




##------------------------------------------AML dataset-----------------------------------------

##------------------------------------------Preparation of reference-----------------------------
if (F){
reference = readRDS(paste0(data_dir, "/seurat_main_cohort.rds"))
alldata.list <- SplitObject(reference, split.by = "cohort")
for (i in 1:length(alldata.list)) {
    alldata.list[[i]] <- FindVariableFeatures(alldata.list[[i]], selection.method = "vst",
        nfeatures = 5000, verbose = FALSE)
}
alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:30,
    reduction = "cca", anchor.features = 5000)
reference <- IntegrateData(anchorset = alldata.anchors, dims = 1:30, new.assay.name = "CCA")
reference = ScaleData(reference, verbose = FALSE)
reference <- RunPCA(reference, npcs = 30, verbose = FALSE)
reference <- RunUMAP(reference, dims = 1:30)
saveRDS(reference, paste0(processed_data_dir, "/aml_reference_integrated_cca.rds"))
}

# Load reference dataset
if (F){
reference = readRDS(paste0(processed_data_dir, "/aml_reference_integrated_cca.rds"))
reference = subset(reference, subset = (cohort == 'Cohort A') | (cohort == 'Triana et al. (healthy)'))
reference <- subset(reference, subset = ct_simple %in% c("Early myeloid", "Immature", "Erythroid")) 
reference = RunPCA(reference, npcs = 30, verbose = FALSE)
reference = RunUMAP(reference, reduction = "pca", return.model = TRUE, dims=1:30)


##------------------------------------------Load and preprocess D2N dataset---------------------
d2n_full <- readRDS(file = paste0(data_dir, "/d2n_SeuratObj_PostQC_dCas9_NewMD.RDS"))
counts <- GetAssayData(d2n_full, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('CRISPRi', 'CRISPRa'))),]
d2n <- subset(d2n_full, features = rownames(counts))

d2n = NormalizeData(d2n)
d2n = ScaleData(d2n)
d2n[["CCA"]] = d2n[["RNA"]]
DefaultAssay(d2n) = "CCA"
VariableFeatures(d2n) = rownames(d2n[["CCA"]])

##------------------------------------------Integration------------------------------------------
anchors <- FindTransferAnchors(
  reference = reference,
  query = d2n,
  k.filter = NA,
  reference.reduction = "pca",
  reference.assay = "CCA",
  query.assay = "CCA",
  features = intersect(rownames(x = reference), rownames(x = d2n)),
  dims = 1:30)

d2n <- MapQuery(
  anchorset = anchors,
  query = d2n,
  reference = reference,
  refdata = list(
    ct = "ct",
    status = "status",
    ct_simple = "ct_simple"),
  reference.reduction = "pca",
  reduction.model = "umap")

###------------------------------------------UMAP for shared embeddings---------------------------
d2n$ct <- d2n$predicted.ct
d2n$status = d2n$predicted.status
d2n$ct_simple = d2n$predicted.ct_simple
d2n$group <- 'query'
reference$group <- 'reference'

emb.merge <- merge( reference[['pca']], d2n[['ref.pca']])
emb.merge@assay.used <- 'CCA'
DefaultAssay(reference) <- "CCA"
DefaultAssay(d2n) <- 'CCA'


obj.merge <- merge(DietSeurat(reference, assays = 'CCA'), DietSeurat(d2n, assays = 'CCA') )
obj.merge[['pca']] <- emb.merge
obj.merge <- RunUMAP(obj.merge, dims = 1:30, reduction="pca")



#------------------------------------------Diffusion map analysis-----------------------------------
df.map <- DiffusionMap(data = obj.merge[['pca']]@cell.embeddings[,1:30])
obj.merge[['DC']] <- CreateDimReducObject(embeddings = df.map@eigenvectors*1000, key = 'DC_', assay = 'CCA')

d2n_types = unique(d2n$ct)
ref_types = setdiff(unique(reference$ct), d2n_types)
obj.merge$ct = factor(obj.merge$ct, levels=c(d2n_types, ref_types))
saveRDS(obj.merge, paste0(processed_data_dir, "/integrated_SeuratObject_AML_d2n.RDS"))
}

obj.merge = readRDS(paste0(processed_data_dir, "/integrated_SeuratObject_AML_d2n.RDS"))
DefaultAssay(obj.merge) <- "CCA"

cols = DiscretePalette(24)
Idents(obj.merge) <- obj.merge@meta.data$group
to_plot_obj.merge = subset(obj.merge, downsample=10000)


#------------------------------------------Panel A---------------------------------------------------
p1 <- DimPlot(to_plot_obj.merge, cols = alpha(cols, 0.33), group.by = 'ct',  pt.size = 0.01, raster = F) +
	NoLegend() +
        ggtitle("All") +
        xlab("") +
        ylab("UMAP_2")
p2 <- DimPlot(to_plot_obj.merge,  cols = alpha(cols, 0.33), group.by = 'ct', cells = WhichCells(to_plot_obj.merge, expression =  group == "query"),  pt.size = 0.01, raster = F) + 
	NoLegend() +
        ggtitle("Perturbed") +
        xlab("UMAP_1") +
        ylab("")
p3 <- DimPlot(to_plot_obj.merge,   cols = alpha(cols, 0.33), group.by = 'ct', cells = WhichCells(to_plot_obj.merge, expression =  group == "reference"),  pt.size = 0.01, raster = F) + 
	NoLegend() +
        ggtitle("Reference") +
        xlab("") +
        ylab("")
p4 =  DimPlot(obj.merge,  cols = alpha(cols, 0.33), group.by = 'ct', cells = WhichCells(to_plot_obj.merge, expression =  group == "query"), raster = F) + theme(legend.position= "bottom",
																legend.text = element_text(size = 9),
																legend.key.size = unit(7, 'pt'))

leg = cowplot::get_legend(p4)

svg(paste0(plots_dir, "/Panel_A_AML_dataset.svg"), width = 8.3, height=3.1)
p1 + p2 + p3 + plot_layout(ncol = 3)&theme(axis.text = element_text(size = 9),
                                           axis.title = element_text(size = 11),
                                           plot.title = element_text(size = 12))
dev.off()

svg(paste0(plots_dir, "/Panel_A_AML_dataset_legend.svg"), width=8.3, height=0.5)
grid.newpage()
grid.draw(leg)
dev.off()
reference = subset(obj.merge, subset = group == "reference")

cols = DiscretePalette(24)
p1 <- DimPlot(reference, cols = alpha(cols, 0.33), group.by = 'cohort',  pt.size = 0.01, raster = F) +
        NoLegend() +
        xlab("") +
        ylab("UMAP_2")
p2 <- DimPlot(reference, cols = alpha(cols, 0.33), group.by = 'patient',  pt.size = 0.01, raster = F) +
        NoLegend() +
        xlab("") +
        ylab("UMAP_2")
p3 <- DimPlot(reference, cols = alpha(cols, 0.33), group.by = 'ct_simple',  pt.size = 0.01, raster = F) +
        xlab("") +
        ylab("UMAP_2")
svg(paste0(plots_dir, "/new_aml_ref.svg"), width = 8.3, height=3.1)
p1 + p2 + p3 + plot_layout(ncol = 3)&theme(axis.text = element_text(size = 9),
                                           axis.title = element_text(size = 11),
                                           plot.title = element_text(size = 12))
dev.off()

 
##------------------------------------------Pseudotime analysis---------------------------------------
###-----------------------------------------Erythroid Cells-------------------------------------------

if (T){
eryth.monocle = subset(obj.merge, subset = (ct %in% c("Aberrant erythroid",
                                                        "Early erythroid progenitor",
                                                        "Late erythroid progenitor",
                                                        "Erythro-myeloid progenitors")))
eryth.monocle[['umap']] = NULL
eryth.monocle[['umap']] <- eryth.monocle[['DC']]

eryth.monocle[['umap']]@cell.embeddings <- eryth.monocle[['umap']]@cell.embeddings[,1:2]

##------------------------------------------Select root cells------------------------------------------
early_eryth = as.data.table(subset(eryth.monocle, subset = (ct== "Erythro-myeloid progenitors"))[['umap']]@cell.embeddings[,1:2])

lim_x = quantile(early_eryth$umap_1, probs=0.85)
lim_y =  quantile(early_eryth$umap_2, probs=0.85)

early_eryth$cell = colnames(subset(eryth.monocle, subset = (ct== "Erythro-myeloid progenitors")))

early_eryth %>%
        filter(umap_1 > lim_x & umap_2 < lim_y) -> init_cells
##------------------------------------------Trajectory inference------------------------------------------
cds <- as.cell_data_set(x = eryth.monocle, assay = "CCA", reductions = "umap")
cds <- cluster_cells(cds = cds, reduction_method = "UMAP", k=15)
cds <- learn_graph(cds)


##------------------------------------------Pseudotime Calculation-----------------------------------------
cds <- order_cells(cds, root_cells = init_cells$cell)

quantile(cds@principal_graph_aux@listData$UMAP$pseudotime)

cds@principal_graph_aux@listData$UMAP$pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime/max(cds@principal_graph_aux@listData$UMAP$pseudotime)
eryth.monocle@meta.data$DC.time <-   cds@principal_graph_aux@listData$UMAP$pseudotime
eryth.monocle@meta.data$label = paste0(eryth.monocle@meta.data$ct, "_", eryth.monocle@meta.data$group)
eryth.monocle@meta.data$m_cluster = cds@clusters$UMAP$clusters
saveRDS(eryth.monocle, paste0(processed_data_dir, "/integrated_erythro_SeuratObject_AML_full.RDS"))
}
eryth.monocle = readRDS(paste0(processed_data_dir, "/integrated_erythro_SeuratObject_AML_full.RDS"))
d2n = subset(eryth.monocle, subset = group == "query")
reference = subset(eryth.monocle, subset = group == "reference")
DefaultAssay(d2n) = "CCA"


##------------------------------------------Panel B--------------------------------------------------------
p2 = plot_cells(cds,
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



p3 = DimPlot(eryth.monocle, sizes.highlight = .01, label=F, raster = F, pt.size = 0.5,  reduction = 'DC', group.by = "ct") +
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

p4 <- DimPlot(eryth.monocle, sizes.highlight = .01, raster = F, pt.size = 0.5,  reduction = 'DC', cells.highlight = init_cells$cell) +
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



svg(paste0(plots_dir, "/Panel_B_AML_dataset_full.svg"), width = 8.3, height=3.1)
p3 + p4  + p2 + plot_layout(ncol = 3)&theme(axis.text = element_text(size = 9),
                                                                axis.title = element_text(size = 11),
                                                                plot.title = element_text(size = 12, face = "bold"))
dev.off()


##------------------------------------------Time dist by dosage gene-----------------------------------
d2n_logFC = readRDS(paste0(data_dir, "/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS"))

d2n_true_erythro = d2n
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
                                                                                                                        distinct(),
			std_time = sd(d2n_true_erythro_norm$DC.time[ which(d2n_true_erythro_norm$guide_crispr == x) ], na.rm = T))
  colnames(out) = c("gene", "time", "expr", "std_time")
  out
}))

time.list$gene = factor(time.list$gene, levels = c("GFI1B", "NFE2", "MYB"))

##------------------------------------------Panel D-------------------------------------------------

p = ggplot(time.list, aes(x = expr, y = time)) +
        geom_hex() +
        stat_bin2d(bins = 50, show.legend = F) +
        scale_fill_gradient(low = "grey85", high = "grey30", guide = guide_colorbar(direction = "horizontal", barheight = 0.5))  +
        geom_smooth(data = time.list %>% filter(!expr == 0.0),
                        method=lm, alpha=0.5, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free_x") +
        theme_bw() +
        ylim(c(0.0, 1.0)) +
        xlab("Cis gene normalized UMI expression") +
        ylab("Pseudotime") +
        theme(legend.key = element_blank(),
                legend.position = "bottom",
                strip.background = element_rect(colour="white", fill="white"),
                axis.text = element_text(size = 9),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))

svg(paste0(plots_dir, "/Panel_D_AML_dataset_d2n_full.svg"), width = 8.3, height=3.85)
p
dev.off()






##------------------------------------------Panel E------------------------------------------
time.list_pseudo %>% filter(gene %in% gene.set) -> time.list_pseudo
time.list_pseudo$gene = factor(time.list_pseudo$gene, levels = c("GFI1B", "NFE2", "MYB"))

p = ggplot(time.list_pseudo, aes(y = time, x = expr, ymin = time - 1 * std_time, ymax = time + 1 * std_time)) +
        geom_point(size = 1, alpha = 0.5) +
        geom_linerange(alpha=0.25, size=0.25) +        
	geom_smooth(method=loess, alpha = 0.25, size = 0.75, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free_x") +
        theme_bw() +
	coord_cartesian(ylim = c(0.0, 1.0), expand = FALSE) +
        xlab("Cis gene log2(FC)") +
        ylab("Mean pseudotime") +
        theme(axis.text = element_text(size = 9),
                legend.key = element_blank(),
                strip.background = element_rect(colour="white", fill="white"),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))

png(paste0(plots_dir, "/Panel_E_AML_dataset_d2n_full.png"), width = 8.3, height=3.24, res=720, units="in")
p
dev.off()


#------------------------------------------NOT_USED----------------------------------------------------
#----------------------Pseudotime_analysis_for_Ref_dataset_and_expression_averaged_d2n-----------------
if (F){
ncells = 20
d2n_true_erythro_norm_sub = subset(d2n_true_erythro_norm, subset = gene %in% gene.set)
d2n_true_erythro_norm_sub = ScaleData(d2n_true_erythro_norm_sub)
time.list_pseudo_exp <- rbindlist(lapply(gene.set, function(y){
        temp_d2n = subset(d2n_true_erythro_norm, subset = gene == y)
	rbindlist(lapply(seq(1, ncol(d2n_true_erythro_norm_sub@assays$CCA), ncells), function(x) {
          as.data.frame(t(temp_d2n@assays$CCA[y])) %>% arrange(pick(y)) %>% slice(x:(x + ncells)) %>% rownames -> temp_cells
          out = data.table(gene = y,
                                time = mean(temp_d2n@meta.data[temp_cells,]$DC.time, na.rm = T),
                                expr = mean(temp_d2n@assays$CCA[y, temp_cells], na.rm = T),
                                bin = as.integer(x / ncells) + 1,
                                std_time = sd(temp_d2n@meta.data[temp_cells,]$DC.time, na.rm = T),
                                std_expr = sd(temp_d2n@assays$CCA[y, temp_cells], na.rm = T))
          colnames(out) = c("gene", "time", "expr", "bin", "std_time", "std_expr")
          out 
        }))
}))


time.list_pseudo_super = c()
time.list_pseudo_exp$dataset = "D2N"
bind_rows(time.list_pseudo_super, time.list_pseudo_exp) -> time.list_pseudo_super


#------------------------------------------Reference AML---------------------------------------------------
reference_aml = subset(reference, subset = (cohort == "Cohort A") & (patient != "A.0"))
#------------------------------------------Time dist by dosage gene-----------------------------------
gene.set = c("GFI1B", "NFE2", "MYB")

time.list <- rbindlist(lapply(gene.set, function(x) {
  return(data.table(gene = x, time = reference_aml$DC.time, expr = reference_aml@assays$CCA@data[x,], ct = reference_aml@meta.data$ct, status = reference_aml@meta.data$status))
}))

time.list$gene = factor(time.list$gene, levels = c("GFI1B", "NFE2", "MYB"))

#Hex version
p = ggplot(time.list, aes(x = expr, y = time, color = status)) +
        geom_hex() +
        stat_bin2d(bins = 50, show.legend = F) +
        scale_fill_gradient(low = "grey85", high = "grey30", guide = guide_colorbar(direction = "horizontal", barheight = 0.5))  +
        geom_smooth(data = time.list %>% filter(!expr == 0.0),
                        method=lm, alpha=0.5, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free_x") +
        theme_bw() +
        ylim(c(0.0, 1.0)) +
        xlab("Cis gene normalized UMI expression") +
        ylab("Pseudotime") +
        theme(legend.key = element_blank(),
                legend.position = "bottom",
                strip.background = element_rect(colour="white", fill="white"),
                axis.text = element_text(size = 9),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))

svg(paste0(plots_dir, "/Panel_D_AML_dataset_ref_aml_full_status.svg"), width = 8.3, height=3.85)
p
dev.off()

#Hex version
p = ggplot(time.list, aes(x = expr, y = time, color = ct)) +
        geom_hex() +
        stat_bin2d(bins = 50, show.legend = F) +
        scale_fill_gradient(low = "grey85", high = "grey30", guide = guide_colorbar(direction = "horizontal", barheight = 0.5))  +
        geom_smooth(data = time.list %>% filter(!expr == 0.0),
                        method=lm, alpha=0.5, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free_x") +
        theme_bw() +
        ylim(c(0.0, 1.0)) +
        xlab("Cis gene normalized UMI expression") +
        ylab("Pseudotime") +
        theme(legend.key = element_blank(),
                legend.position = "bottom",
                strip.background = element_rect(colour="white", fill="white"),
                axis.text = element_text(size = 9),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))

svg(paste0(plots_dir, "/Panel_D_AML_dataset_ref_aml_full_ct.svg"), width = 8.3, height=3.85)
p
dev.off()

ncells = 20
time.list_pseudo <- rbindlist(lapply(gene.set, function(y){
        rbindlist(lapply(seq(1, ncol(reference_aml@assays$CCA), ncells), function(x) {
          as.data.frame(t(reference_aml@assays$CCA[y])) %>% arrange(pick(y)) %>% slice(x:(x + ncells)) %>% rownames -> temp_cells 
          out = data.table(gene = y,
                                time = mean(reference_aml@meta.data[temp_cells,]$DC.time, na.rm = T),
                                expr = mean(reference_aml@assays$CCA[y, temp_cells], na.rm = T),
                                bin = as.integer(x / ncells) + 1,
				std_time = sd(reference_aml@meta.data[temp_cells,]$DC.time, na.rm = T),
				std_expr = sd(reference_aml@assays$CCA[y, temp_cells], na.rm = T))
	  colnames(out) = c("gene", "time", "expr", "bin", "std_time", "std_expr")
          out
        }))
}))

##------------------------------------------Panel E------------------------------------------
time.list_pseudo %>% filter(gene %in% gene.set) -> time.list_pseudo
time.list_pseudo$gene = factor(time.list_pseudo$gene, levels = c("GFI1B", "NFE2", "MYB"))

p = ggplot(time.list_pseudo, aes(y = time, x = expr, ymin = time - 1 * std_time, ymax = time + 1 * std_time, xmin = expr - 1 * std_expr, xmax = expr + 1 * std_expr)) +
        geom_point(size = 1, alpha = 0.5) +
	geom_linerange(alpha=0.25, size=0.25) +
	ggstance::geom_linerangeh(alpha=0.25, size=0.25) +
        geom_smooth(method=loess, alpha = 0.25, size = 0.75, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free_x") +
        theme_bw() +
	coord_cartesian(ylim = c(0.0, 1.0), expand = FALSE) +
        xlab("Aggregated cis gene UMI expression") +
        ylab("Mean pseudotime") +
        theme(axis.text = element_text(size = 9),
                legend.key = element_blank(),
                strip.background = element_rect(colour="white", fill="white"),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))

png(paste0(plots_dir, "/Panel_E_AML_dataset_ref_aml_full.png"), width = 8.3, height=3.24, units="in", res=720)
p
dev.off()

time.list_pseudo$dataset = "AML"
bind_rows(time.list_pseudo_super, time.list_pseudo) -> time.list_pseudo_super

#------------------------------------------Reference Healthy---------------------------------------------------
reference_healthy = subset(reference, subset = (patient == "Reference") | (patient == "A.0"))
##------------------------------------------Time dist by dosage gene-----------------------------------

gene.set = c("GFI1B", "NFE2", "MYB")

time.list <- rbindlist(lapply(gene.set, function(x) {
  return(data.table(gene = x, time = reference_healthy$DC.time, expr = reference_healthy@assays$CCA@data[x,]))
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
        ylim(c(0.0, 1.0)) +
        xlab("Cis gene normalized UMI expression") +
        ylab("Pseudotime") +
        theme(legend.key = element_blank(),
                legend.position = "bottom",
                strip.background = element_rect(colour="white", fill="white"),
                axis.text = element_text(size = 9),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))

svg(paste0(plots_dir, "/Panel_D_AML_dataset_ref_healthy_full.svg"), width = 8.3, height=3.85)
p
dev.off()


time.list_pseudo <- rbindlist(lapply(gene.set, function(y){
        rbindlist(lapply(seq(1, ncol(reference_healthy@assays$CCA), ncells), function(x) {
          as.data.frame(t(reference_healthy@assays$CCA[y])) %>% arrange(pick(y)) %>% slice(x:(x + ncells)) %>% rownames -> temp_cells
          out = data.table(gene = y,
                                time = mean(reference_healthy@meta.data[temp_cells,]$DC.time, na.rm = T),
                                expr = mean(reference_healthy@assays$CCA[y, temp_cells], na.rm = T),
                                bin = as.integer(x / ncells) + 1,
                                std_time = sd(reference_healthy@meta.data[temp_cells,]$DC.time, na.rm = T),
                                std_expr = sd(reference_healthy@assays$CCA[y, temp_cells], na.rm = T))
          colnames(out) = c("gene", "time", "expr", "bin", "std_time", "std_expr")
          out
        }))
}))


##------------------------------------------Panel E------------------------------------------
time.list_pseudo %>% filter(gene %in% gene.set) -> time.list_pseudo
time.list_pseudo$gene = factor(time.list_pseudo$gene, levels = c("GFI1B", "NFE2", "MYB"))


p = ggplot(time.list_pseudo, aes(y = time, x = expr, ymin = time - 1 * std_time, ymax = time + 1 * std_time, xmin = expr - 1 * std_expr, xmax = expr + 1 * std_expr)) +
        geom_point(size = 1, alpha = 0.5) +
        geom_linerange(alpha=0.25, size=0.25) +
        ggstance::geom_linerangeh(alpha=0.25, size=0.25) +
        geom_smooth(method=loess, alpha = 0.25, size = 0.75, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free_x") +
        theme_bw() +
        coord_cartesian(ylim = c(0.0, 1.0), expand = FALSE) +
        xlab("Aggregated cis gene UMI expression") +
        ylab("Mean pseudotime") +
        theme(axis.text = element_text(size = 9),
                legend.key = element_blank(),
                strip.background = element_rect(colour="white", fill="white"),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))


png(paste0(plots_dir, "/Panel_E_AML_dataset_ref_healthy_full.png"), width = 8.3, height=3.24, units="in", res=720)
p
dev.off()

time.list_pseudo$dataset = "Healthy"
bind_rows(time.list_pseudo_super, time.list_pseudo) -> time.list_pseudo_super







##-------------------------------------------Super Panel E---------------------------------------------------------------
time.list_pseudo_super %>% filter(gene %in% gene.set) -> time.list_pseudo_super
time.list_pseudo_super$gene = factor(time.list_pseudo_super$gene, levels = c("GFI1B", "NFE2", "MYB"))
time.list_pseudo_super$dataset = factor(time.list_pseudo_super$dataset, levels = c("AML", "Healthy"))




p = ggplot(time.list_pseudo_super, aes(y = time, x = expr, ymin = time - 1 * std_time, ymax = time + 1 * std_time, xmin = expr - 1 * std_expr, xmax = expr + 1 * std_expr)) +
        geom_point(size = 1, alpha = 0.5) +
        geom_linerange(alpha=0.25, size=0.25) +
        #ggstance::geom_linerangeh(alpha=0.25, size=0.25) +
        geom_smooth(method=loess, alpha = 0.25, size = 0.75, color="#FF7F00") +
        facet_grid(dataset ~ gene) +
        theme_bw() +
        xlab("Mean cis gene normalized UMI expression") + 
        ylab("Mean Pseudotime") +
        coord_cartesian(ylim = c(0.0, 1.0), expand = T) +
        theme(axis.text = element_text(size = 9), 
                legend.key = element_blank(), 
                strip.background = element_rect(colour="white", fill="white"), 
                axis.title = element_text(size = 11), 
                strip.text = element_text(size = 11), 
                plot.title = element_text(size = 12, face = "bold")) 
 
 
pdf(paste0(plots_dir, "/Panel_E_AML_dataset_full_super_ref.pdf"), width = 8.3,height=5.6)#, units="in", res=720) 
p
dev.off() 
}
