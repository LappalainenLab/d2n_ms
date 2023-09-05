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
data_dir = file.path("../../data", "09-integration")
processed_data_dir = file.path("../../processed_data", "09-integration")
plots_dir = file.path("../../plots", folder_name)



#-------------------------------------------Dosage<->Disease_analysis-----------------------------------------
##------------------------------------------leukemic_vs_Progenitors-------------------------------------------
###-----------------------------------------leukemic--------------------------------------------------------
reference = readRDS(paste0(data_dir, "/seurat_main_cohort.rds"))
reference = subset(reference, subset = (cohort == "Cohort A") | (cohort == "Triana et al. (healthy)"))
reference = subset(reference, subset = ct_simple %in% c("Immature", "Early myeloid", "Monocytes"))

#------------------------------------------Diffusion map analysis-----------------------------------  
df.map <- DiffusionMap(data = reference[['scanorama']]@cell.embeddings[,1:30])
reference[['DC']] <- CreateDimReducObject(embeddings = df.map@eigenvectors*1000, key = 'DC_', assay = 'RNA')
 
cts = setdiff(unique(reference@meta.data$ct), c("Pre-pro-B cells", "Lymphomyeloid prog"))
myelo.monocle = subset(reference, subset = (ct %in% cts) & (id != 29))

myelo.monocle[['umap']] = NULL
myelo.monocle[['umap']] <- myelo.monocle[['DC']]

myelo.monocle[['umap']]@cell.embeddings <- myelo.monocle[['umap']]@cell.embeddings[,1:2]

##------------------------------------------Select root cells------------------------------------------

early_myelo = as.data.table(subset(myelo.monocle, subset = ct == "HSCs & MPPs")[['umap']]@cell.embeddings[,1:2])
early_myelo$cell = colnames(subset(myelo.monocle, subset = (ct == "HSCs & MPPs")))
lim_x = quantile(early_myelo$umap_1, probs=0.95)
lim_y =  quantile(early_myelo$umap_2, probs=0.05)

early_myelo$cell = colnames(subset(myelo.monocle, subset = (ct== "HSCs & MPPs")))

early_myelo %>%
        filter(umap_1 > lim_x & 
		umap_2 < lim_y) -> init_cells
##------------------------------------------Trajectory inference------------------------------------------
cds <- as.cell_data_set(x = myelo.monocle, assay = "RNA", reductions = "umap")
cds <- cluster_cells(cds = cds, reduction_method = "UMAP", k=15)
cds <- learn_graph(cds)

##------------------------------------------Pseudotime Calculation-----------------------------------------
cds <- order_cells(cds, root_cells = init_cells$cell)

quantile(cds@principal_graph_aux@listData$UMAP$pseudotime)

cds@principal_graph_aux@listData$UMAP$pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime/max(cds@principal_graph_aux@listData$UMAP$pseudotime)
myelo.monocle@meta.data$DC.time <-   cds@principal_graph_aux@listData$UMAP$pseudotime
myelo.monocle@meta.data$label = paste0(myelo.monocle@meta.data$ct, "_", myelo.monocle@meta.data$group)
myelo.monocle@meta.data$m_cluster = cds@clusters$UMAP$clusters

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
 
p3 = DimPlot(myelo.monocle, sizes.highlight = .01, label=F, raster = F, pt.size = 0.5,  reduction = 'DC', group.by = "ct") +
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


p4 <- DimPlot(myelo.monocle, sizes.highlight = .01, raster = F, pt.size = 0.5,  reduction = 'DC', cells.highlight = init_cells$cell) +
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
 

 
png(paste0(plots_dir, "/Panel_B_AML_dataset_full.png"), width = 8.3, height=3.24, res=720, units="in")
p3 + p4  + p2 + plot_layout(ncol = 3)&theme(axis.text = element_text(size = 9),
                                                                axis.title = element_text(size = 11),
                                                                plot.title = element_text(size = 12, face = "bold"))
dev.off()


reference_leukemic = subset(myelo.monocle, subset = status == "leukemic")
Idents(reference_leukemic) = reference_leukemic@meta.data$ct
reference_leukemic = subset(reference_leukemic, downsample = 300)
gene.set = c("NPM1", "FLT3", "IDH2", "CEBPA", "RUNX1", "IDH1")
time.list_pseudo_super = c()


time.list <- rbindlist(lapply(gene.set, function(x) {
  return(data.table(gene = x, time = reference_leukemic$DC.time, expr = reference_leukemic@assays$RNA@data[x,], ct = reference_leukemic@meta.data$ct, status = reference_leukemic@meta.data$status))
}))
 

time.list$gene = factor(time.list$gene, levels = gene.set)



ncells = 25
time.list_pseudo_aml <- rbindlist(lapply(gene.set, function(y){
        rbindlist(lapply(seq(1, ncol(reference_leukemic@assays$RNA), ncells), function(x) {
          # Arrange by expr
	  as.data.frame(t(reference_leukemic@assays$RNA[y])) %>% arrange(pick(y)) %>% slice(x:(x + ncells)) %>% rownames -> temp_cells
          
	  # Arrange by time
	  #reference_leukemic@meta.data %>% arrange(DC.time) %>% slice(x:(x + ncells)) %>% rownames -> temp_cells
	  out = data.table(gene = y,
                                time = mean(reference_leukemic@meta.data[temp_cells,]$DC.time, na.rm = T),
                                expr = mean(reference_leukemic@assays$RNA[y, temp_cells], na.rm = T),
                                bin = as.integer(x / ncells) + 1,
				ct = reference_leukemic@meta.data[temp_cells,]$ct %>% table %>% as.data.frame %>% arrange(desc(Freq)) %>% slice(1) %>% select(".") %>% as.vector,
                                std_time = sd(reference_leukemic@meta.data[temp_cells,]$DC.time, na.rm = T),
                                std_expr = sd(reference_leukemic@assays$RNA[y, temp_cells], na.rm = T))
          colnames(out) = c("gene", "time", "expr", "bin", "ct", "std_time", "std_expr")
          out
        }))
}))

##------------------------------------------Panel E------------------------------------------
time.list_pseudo_aml %>% filter(gene %in% gene.set) -> time.list_pseudo_aml
time.list_pseudo_aml$gene = factor(time.list_pseudo_aml$gene, levels = gene.set)
 

p = ggplot(time.list_pseudo_aml, aes(y = time, x = expr, color=as.factor(unlist(ct)), label=as.factor(unlist(ct)))) +
        geom_point(size = 1, alpha = 0.75) +
        geom_smooth(method=loess, alpha = 0.25, size = 0.75, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free") +
        theme_bw() +
	#ggrepel::geom_text_repel() +
	#ylim(c(0.0, 1.0)) +
        xlab("Aggregated cis gene UMI expression") +
        ylab("Mean pseudotime") +
        labs(color="ct") +
	theme(axis.text = element_text(size = 9),
                #legend.position='none',
                strip.background = element_rect(colour="white", fill="white"),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))
 
png(paste0(plots_dir, "/Panel_E_AML_dataset_ref_leukemic_expr.png"), width = 8.3, height=3.24, res=720, units="in")
p 
dev.off()

time.list_pseudo_aml$dataset = "AML" 
bind_rows(time.list_pseudo_super, time.list_pseudo_aml) -> time.list_pseudo_super 


###-----------------------------------------Healthy--------------------------------------------------------
reference_prog = subset(myelo.monocle, subset = patient == "Reference")

time.list <- rbindlist(lapply(gene.set, function(x) {
  return(data.table(gene = x, time = reference_prog$DC.time, expr = reference_prog@assays$RNA@data[x,], ct = reference_prog@meta.data$ct, status = reference_prog@meta.data$status))
}))


time.list$gene = factor(time.list$gene, levels = gene.set)

ncells = 25
time.list_pseudo <- rbindlist(lapply(gene.set, function(y){
        rbindlist(lapply(seq(1, ncol(reference_prog@assays$RNA), ncells), function(x) {
	  # Arrange by expr 
          as.data.frame(t(reference_prog@assays$RNA[y])) %>% arrange(pick(y)) %>% slice(x:(x + ncells)) %>% rownames -> temp_cells
           
          # Arrange by time
          #reference_prog@meta.data %>% arrange(DC.time) %>% slice(x:(x + ncells)) %>% rownames -> temp_cells
	  out = data.table(gene = y,
                                time = mean(reference_prog@meta.data[temp_cells,]$DC.time, na.rm = T),
                                expr = mean(reference_prog@assays$RNA[y, temp_cells], na.rm = T),
                                bin = as.integer(x / ncells) + 1,
				ct = reference_prog@meta.data[temp_cells,]$ct %>% table %>% as.data.frame %>% arrange(desc(Freq)) %>% slice(1) %>% select(".") %>% as.vector,
                                std_time = sd(reference_prog@meta.data[temp_cells,]$DC.time, na.rm = T),
                                std_expr = sd(reference_prog@assays$RNA[y, temp_cells], na.rm = T))
          colnames(out) = c("gene", "time", "expr", "bin", "ct", "std_time", "std_expr")
          out
        }))
}))


##------------------------------------------Panel E------------------------------------------
time.list_pseudo %>% filter(gene %in% gene.set) -> time.list_pseudo
time.list_pseudo$gene = factor(time.list_pseudo$gene, levels = gene.set)

p = ggplot(time.list_pseudo, aes(y = time, x = expr, color=as.factor(unlist(ct)), label=as.factor(unlist(ct)))) + 
        geom_point(size = 1, alpha = 0.75) +
        geom_smooth(method=loess, alpha = 0.25, size = 0.75, color="#FF7F00") +
        facet_wrap(~ gene, scale = "free") +
        theme_bw() +
        #ylim(c(0.6, 3.15)) +
        xlab("Aggregated cis gene UMI expression") +
        ylab("Mean pseudotime") +
        theme(axis.text = element_text(size = 9),
                legend.key = element_blank(),
                strip.background = element_rect(colour="white", fill="white"),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))
 
png(paste0(plots_dir, "/Panel_E_AML_dataset_ref_prog_expr.png"), width = 8.3, height=3.24, res=720, units="in")
p
dev.off()


time.list_pseudo$dataset = "Healthy"
bind_rows(time.list_pseudo_super, time.list_pseudo) -> time.list_pseudo_super


##-------------------------------------------Super Panel E---------------------------------------------------------------
time.list_pseudo_super %>% filter(gene %in% gene.set) -> time.list_pseudo_super
time.list_pseudo_super$gene = factor(time.list_pseudo_super$gene, levels = gene.set)
time.list_pseudo_super$dataset = factor(time.list_pseudo_super$dataset, levels = c("AML", "Healthy"))

p = ggplot(time.list_pseudo_super, aes(y = time, x = expr, ymin = time - 1 * std_time, ymax = time + 1 * std_time, xmin = expr - 1 * std_expr, xmax = expr + 1 * std_expr)) +
        geom_point(size = 1, alpha = 0.5) +
        geom_linerange(alpha=0.25, size=0.25) +
        #ggstance::geom_linerangeh(alpha=0.25, size=0.25) +
        geom_smooth(method="loess", alpha = 0.25, size = 0.75, color="#FF7F00", span=0.8) +
        facet_grid(dataset ~ gene, scales="free_x") +
        theme_bw() +
        xlab(" cis gene expression measure") +
        ylab("Mean Pseudotime") +
	#coord_cartesian(ylim = c(-0.5, 1.0), expand = FALSE) +
        #coord_cartesian(xlim = c(0.0, 0.5), expand = FALSE) +
	theme(axis.text = element_text(size = 9),
                legend.key = element_blank(),
                strip.background = element_rect(colour="white", fill="white"),
                axis.title = element_text(size = 11),
                strip.text = element_text(size = 11),
                plot.title = element_text(size = 12, face = "bold"))


png(paste0(plots_dir, "/Panel_E_AML_dataset_full_super_myelo_ref_expr.png"), width = 10, height=8.3, units="in", res=720)
p
dev.off()
