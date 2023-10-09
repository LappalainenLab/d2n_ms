##### Plot clustered heat map updated SSv1 data for selected target genes in d2n project#####
# JDE, Nov 2021

library(data.table)
library(ggplot2)
library(ggdendro)
library(RColorBrewer)
library(gridExtra)
library(WGCNA)


# Directories
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/001-library_design/")
processed_data_dir = "/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/001-GWAS_CRISPRi/202108_coexpression_SSv1/processed_data/"

# source functions
functions_list <- list.files("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/001-GWAS_CRISPRi/202108_coexpression_SSv1/functions/")
invisible(sapply(paste0("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/001-GWAS_CRISPRi/202108_coexpression_SSv1/functions/", functions_list), source, .GlobalEnv))


### Load data
# Updated Expression matrices
load(file = file.path(processed_data_dir, "20210907_Exp_l.RData"))

# Selected target genes
SG <- fread(file = "processed_data/20211111_list_genes_targeted_seq.txt")

### Get dendogram and heatmap dataframes
snp_vec <- c("snp63", "snp83")


CM_SNP_l <- lapply(snp_vec, function(s){
  list(MAGIC = ClustCorrMatrixFromTransfExpr(exp_matrix = Exp_l[[s]][["magic"]][["sub"]],
                                             hclust_method = "ward.D2"))
})
names(CM_SNP_l) <- snp_vec

# Cluster information
Clust_info <- list(
  snp63 = data.table::setnames(fread(file = file.path(processed_data_dir, "SSv1_ClustIDs_SNP63.txt"))[, .(nclust_3, gene)], c("cluster", "gene")),
  snp83 = data.table::setnames(fread(file = file.path(processed_data_dir, "SSv1_ClustIDs_SNP83.txt"))[, .(nclust_4, gene)], c("cluster", "gene"))
)


# Color palettes
pal_clust_snp63_num <-  brewer.pal(4, "Spectral")
names(pal_clust_snp63_num) <- c("1","4","2","3")


pal_clust_snp83_num <-  brewer.pal(4, "Spectral")
names(pal_clust_snp83_num) <- c("1","2","4","3") 

pal_clust_num_l <- list(snp63 = pal_clust_snp63_num,
                    snp83 = pal_clust_snp83_num)

s="snp63"
m="MAGIC"

gene_name <- list(
  snp63 = "GFI1B",
  snp83 = "NFE2"
)

# Iterate over the different SNPs
for (s in snp_vec) {
  
  CM_l <- CM_SNP_l[[s]]
  clust_info <- Clust_info[[s]]
  clust_info$row = 1
  
  # Generate plots with old cluster labels for each method
  for (m in names(CM_l)){
    
    # Dendogram-heatmap plot elements
    cor_df_temp = CM_l[[m]]$Cor_df
    tile_fill = "bicor"
    plot_title = gene_name[[s]]
    
    # new clust object for only selected genes
    Cor <- CM_l[[m]]$Cor
    sc <- colnames(Cor) %in% SG$external_gene_name
    sr <- rownames(Cor) %in% SG$external_gene_name
    Cor_sel <- Cor[sr, sc]
    set.seed(1234)
    clust_object <- stats::hclust(as.dist(1 - Cor_sel), method = "ward.D2")
    
    # Palette
    pal_clust <- pal_clust_num_l[[s]]
    
    # Create the data frame for plotting heatmap
    cor_df <- cor_df_temp[gene_A %in% SG$external_gene_name & gene_B %in% SG$external_gene_name,]
    
    DF <- data.frame(cor_df)
    df <- setNames(DF[, c(tile_fill, "gene_A", "gene_B")], c("fill_value", "gene_A", "gene_B"))
    
    df$gene_A <- factor(df$gene_A, levels = clust_object$labels[clust_object$order])
    df$gene_B <- factor(df$gene_B, levels = clust_object$labels[clust_object$order])
    
    
    #dendogram including information about which gene belongs to which cluster
    dend_df = ggdendro::dendro_data(stats::as.dendrogram(clust_object), type = "rectangle")
    
    dend_data = dend_df
    dend_data$labels <- merge(dend_data$labels, clust_info, by.x="label", by.y="gene", all.x = T)
    dend_segments_temp <- merge(dend_data$segments[dend_data$segments$yend == 0,], 
                                dend_data$labels[, c("x","label", "cluster")], by = "x", all.x = T)
    dend_data$segments <- rbind(dend_segments_temp,
                                cbind(dend_data$segments[!(dend_data$segments$yend == 0),],
                                      data.frame(cluster=NA),
                                      data.frame(label=NA)))
    
    
    p_dend <- ggplot(dend_data$segments) + 
      geom_segment(aes(x = x, 
                       y = y, 
                       xend = xend, 
                       yend = yend)) +
      labs(x="", y="") +
      theme(panel.background = element_rect(color="white", fill="white"),
            axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            axis.title = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) 
    
    
    range_cor <- range(df$fill_value[(df$gene_A != df$gene_B)])
    p_heatmap <- ggplot(df, aes(x=gene_B, y=gene_A)) +
      geom_raster(aes(fill=fill_value)) +
      scale_fill_gradient2(tile_fill, 
                           midpoint = 0, 
                           low = "#9161A8", 
                           high = "#F7941E", 
                           mid = "white", 
                           na.value = "grey90", 
                           limits=range_cor) +
      theme(legend.position = "bottom", legend.direction = "horizontal") +
      guides(fill = guide_colourbar(barwidth = 10, barheight = 0.5)) +
      theme(panel.background = element_rect(color="white", fill="white"),
            axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            axis.title = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_discrete(expand = expansion(mult = c(0.01, 0.01)))
    
    
    clust_df <- clust_info[gene %in% SG$external_gene_name,]
    clust_df$gene <- factor(clust_df$gene, levels = clust_object$labels[clust_object$order])
    
    p_clust <- ggplot(clust_df, aes(x=gene, y=row)) +
      geom_raster(aes(fill=factor(cluster)), show.legend = F) +
      scale_fill_manual("Cluster",
                        values = pal_clust, 
                        na.value = "white") +
      theme(panel.background = element_rect(color="white", fill="white"),
            axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            axis.title = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
      scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
      theme(legend.position = "bottom", legend.direction = "horizontal")
    
    lay <- rbind(c(1,1,1,1,1,1,1,1,1,1),
                 c(1,1,1,1,1,1,1,1,1,1),
                 c(1,1,1,1,1,1,1,1,1,1),
                 c(1,1,1,1,1,1,1,1,1,1),
                 c(1,1,1,1,1,1,1,1,1,1),
                 c(1,1,1,1,1,1,1,1,1,1),
                 c(3,3,3,3,3,3,3,3,3,3),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2),
                 c(2,2,2,2,2,2,2,2,2,2))
    p <- grid.arrange(grobs = list(p_dend, p_heatmap, p_clust), layout_matrix = lay, top=plot_title)
    #ggsave(p, filename = paste0("plots/20210909_DendroHeatmap_oldClusts_", m, "_", s, ".pdf"), width = 8, height = 9.5)
    ggsave(p, filename = paste0("plots/20211129_SelectedGenes_DendroHeatmap_", s, ".png"), width = 8, height = 9.5)
  }
}


### Write tables with cluster category when using different number of clusters
