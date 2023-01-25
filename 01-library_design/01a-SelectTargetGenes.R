##### D2N - Targeted-seq gene selection #####
# (1) Get all genes from the NEF2 and GFI1B trans network, with cluster identity
# (2) Find which genes have 5'ALT isoforms
# (3) Get average absolute correlation with the rest of the genes in the cluster
# (4) Get average expression in wt K562 scRNA-seq
# (5) Add each gene's cumulative sum PPI score with the rest of the genes

# JDE, September 2021
# Last modified: January 2023


# Libraries
library(biomaRt)
library(data.table)
library(stringr)
library(WGCNA)
library(parallel)
library(plyr)
library(ggplot2)



# Dir

setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/001-library_design/")
processed_data = "processed_data"
data_dir = "../../../data"


# Functions
remove_trailing_digit <- function(x) stringr::str_remove(x, "\\..+")

# Variables
num_tg_total = 86
control_genes = c("GAPDH", "LHX3")
titration_genes = c("GFI1B", "NFE2", "MYB", "TET2")

### (0) Load data

# Trans associations with each gene cluster id
GFI1B_trans <- cbind(fread(input = file.path(data_dir, "GWAS-CRISPRi_phaseI/210709_UpdatedAssociations/210813_SSv1-Expression-Trans-SNP-63-Results.txt"))[sig_trans == 1,],
                     data.table::setnames(fread(input = "../../001-GWAS_CRISPRi/202108_coexpression_SSv1/processed_data/SSv1_ClustIDs_SNP63.txt")[, .(nclust_4)], c("cluster")))
NFE2_trans <- cbind(fread(input = file.path(data_dir, "GWAS-CRISPRi_phaseI/210709_UpdatedAssociations/210813_SSv1-Expression-Trans-SNP-83-Results.txt"))[sig_trans == 1,],
                    data.table::setnames(fread(input = "../../001-GWAS_CRISPRi/202108_coexpression_SSv1/processed_data/SSv1_ClustIDs_SNP83.txt")[, .(nclust_4)], c("cluster")))

# Trans networks expression matrix
GFI1B_Exp <- t(as.matrix(fread(file = paste0("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/001-GWAS_CRISPRi/202108_coexpression_SSv1/processed_data/output_magic_snp63.csv"))))
NFE2_Exp <- t(as.matrix(fread(file = paste0("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/001-GWAS_CRISPRi/202108_coexpression_SSv1/processed_data/output_magic_snp83.csv"))))

# Correlation matrices
GFI1B_Cor <- WGCNA::bicor(x = t(GFI1B_Exp), nThreads = detectCores()/2)
NFE2_Cor <- WGCNA::bicor(x = t(NFE2_Exp), nThreads = detectCores()/2)


# Features datasets
GFI1B_features <- merge.data.table(GFI1B_trans[, .(gene, log2fc_norm, qvalue, cells_with_umis, umis_per_cell_norm, cluster)], 
                                   data.table(gene = rownames(GFI1B_Cor), mean_abs_cor = rowMeans(abs(GFI1B_Cor))),
                                   by = "gene")
NFE2_features <- merge.data.table(NFE2_trans[, .(gene, log2fc_norm, qvalue, cells_with_umis, umis_per_cell_norm, cluster)], 
                                  data.table(gene = rownames(NFE2_Cor), mean_abs_cor = rowMeans(abs(NFE2_Cor))),
                                  by = "gene")

Shared_features <- merge.data.table(GFI1B_features, NFE2_features, by="gene", suffixes = c("_gfi1b", "_nfe2"))
Uniq_GFI1B_features <- GFI1B_features[! (gene %in% Shared_features$gene), ]
Uniq_NFE2_features <- NFE2_features[! (gene %in% Shared_features$gene), ]

# Unique list of genes
universe_target_gene_symbols <- unique(c(GFI1B_features$gene, NFE2_features$gene))
# rename a pair of genes not found by BioMart
universe_target_gene_symbols <- gsub("H3F3B", "H3-3B", universe_target_gene_symbols) 

# 5'Alternative structure genes list
alt5utr_all <- fread(input = file.path(data_dir, "ONT_Glinos_2021/K562_expressed_transcripts_alt5utrs.txt"))
alt5utr_k562 <- fread(input = file.path(data_dir, "ONT_Glinos_2021/K562_expressed_transcripts_alt5utrs_within.txt"))
alt5utr_all[, ensembl_gene_id := remove_trailing_digit(gene_id)]
alt5utr_k562[, ensembl_gene_id := remove_trailing_digit(gene_id)]


# K562 scRNA-seq data
ECCITE_K562 <- fread(input = file.path(data_dir, "K562_RNAseq/K562_scRNAseq_bulkRNAseq.txt"))


# Get ENSEMBL ids
# Genereate Ensembl Mart
ensembl <- biomaRt::useEnsembl(
  biomart = "ensembl", 
  dataset = "hsapiens_gene_ensembl",
  host = "http://www.ensembl.org"
)
# Use 'attributes' to tell biomaRt what values to return; in this case, return a data frame with the Ensembl gene IDs and the gene symbols
# Use 'filters' to tell biomaRt what to search by; in this case, search by Gene Symbol ('external_gene_name')
# Use 'values' to input the search query; in this case, a vector of Ensembl gene IDs from a data frame I have in R
symbols <- as.data.table(biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = universe_target_gene_symbols,
  mart = ensembl
)) # some gene symbols have more than one Ensembl id

# Merge Ensmbl ids with 5'Alt structure
symbols[, alt5utr_all := ensembl_gene_id %in% alt5utr_all$ensembl_gene_id]


# Find genes with multiple Ensembl IDs 
gene_ensembl_count <- data.table::setnames(data.table(table(symbols[, external_gene_name])), c("gene", "ensembl_count"))
genes_multiple_ensembl <- gene_ensembl_count[ensembl_count > 1, gene]

symbols_unique <- symbols[, list(any_alt5utr = any(alt5utr_all), 
                                 ensembl_id_list = paste0(ensembl_gene_id, collapse = ";"),
                                 multiple_ensembl = length(ensembl_gene_id)>1), by=list(gene = external_gene_name)]

universe_target_gene_symbols_noAlt5 <- symbols_unique[any_alt5utr == F & multiple_ensembl == F, gene]

## Add extra info related to TFs and PPI

# Load in-house ADTs
ADT <- fread(file = file.path(data_dir, "Antibodies_TotalSeq/TSC_BAP_v2.txt"))

# Load TF data
TFS <- fread(file = file.path(data_dir, "TFs/GarciaAlonso_GenRes_2019/merged_interactions.txt"), header = T) %>%
  dplyr::select(TF, target) %>%
  dplyr::group_by(TF) %>%
  dplyr::distinct(target, .keep_all = TRUE) %>%
  dplyr::mutate(
    N_TfTargets = length(target),
    N_TfTargets_in_transnet = sum(target %in% universe_target_gene_symbols_noAlt5),
    TfTargets_in_transnet = paste(target[target %in% universe_target_gene_symbols_noAlt5], collapse = "__"),
    N_ADTsTargets = sum(target %in% ADT$Gene)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(TF, N_TfTargets, N_TfTargets_in_transnet, TfTargets_in_transnet, N_ADTsTargets) %>%
  dplyr::distinct(TF, .keep_all = TRUE) %>%
  dplyr::rename(TfGene = TF) %>%
  as.data.table()

# Load PPI data
PPI_exp <- fread(file = file.path(data_dir, "STRING/processed_data/STRING_PPI_physical_all.txt"))
FindPPIinTransnet <- function(x) {
  sum(PPI_exp[gene1 == x, gene2] %in% universe_target_gene_symbols_noAlt5)
}
symbols_unique[, PPI_in_transnet := FindPPIinTransnet(gene), gene]

# Merge TFs data with rest of symbol information
symbols_unique_features <- merge.data.table(symbols_unique, TFS, all.x = T, by.x = "gene", by.y = "TfGene")


Shared_features_symbols <- merge.data.table(Shared_features, symbols_unique_features, by="gene", all.x = T)
Uniq_GFI1B_features_symbols <- merge.data.table(Uniq_GFI1B_features, symbols_unique_features, by="gene", all.x = T)
Uniq_NFE2_features_symbols <- merge.data.table(Uniq_NFE2_features, symbols_unique_features, by="gene", all.x = T)


# Select all shared transnetwork genes that are shared but don't have Alt 5'
Selected_shared_genes <- Shared_features_symbols[any_alt5utr == F & multiple_ensembl == F & !(gene %in% titration_genes), ]


# Find number of targets per cluster
N_clust <- setnames(rbind(cbind(as.data.table(table(Uniq_GFI1B_features$cluster)), rep("GFI1B", "4")),
                          cbind(as.data.table(table(Uniq_NFE2_features$cluster)), rep("NFE2", "4"))), 
                    c("cluster", "n", "network")) 
n_shared <- nrow(Selected_shared_genes)
n_rand <- num_tg_total - n_shared

N_clust[, pct := n/sum(N_clust$n)]
N_clust[cluster %in% 1:3 & network == "GFI1B", n_target := floor(pct*n_rand)]
N_clust[!(cluster %in% 1:3 & network == "GFI1B"), n_target := ceiling(pct*n_rand)]
N_clust[network == "GFI1B", adj_n_target := n_target - 1]
N_clust[network == "NFE2", adj_n_target := n_target + 1]

sum(N_clust$n_target)

# Find number of TFs to select per cluster
Net_L <- list(GFI1B = Uniq_GFI1B_features_symbols[ multiple_ensembl == F & any_alt5utr == F & !(gene %in% titration_genes),],
              NFE2 = Uniq_NFE2_features_symbols[multiple_ensembl == F & any_alt5utr == F & !(gene %in% titration_genes), ])

plot_dt <- rbind(cbind(Net_L[["GFI1B"]][ !is.na(N_TfTargets_in_transnet),], data.frame(network = "GFI1B")), 
                 cbind(Net_L[["NFE2"]][ !is.na(N_TfTargets_in_transnet),], data.frame(network = "NFE2")))
p <- ggplot(plot_dt, aes(x=N_TfTargets_in_transnet, y=umis_per_cell_norm*10000)) +
  geom_point(aes(color=network)) +
  scale_y_log10() +
  theme_bw() +
  geom_vline(lty=2, xintercept = 10) +
  geom_hline(lty=2, yintercept = 0.1) +
  labs(y = "Normalized UMIs per cell", x = "# TF targets in trans-network")
p
ggsave("../001-library_design/plots/20211111_TransTFs_UMIsCell_NTargets.pdf", p, width = 6, height = 4)


# A third of the target genes to be TFs
N_clust[, n_TFs := ceiling(n_target/3)]

# Make sure there are enough TFs in represented in each cluster
Uniq_GFI1B_features_symbols[multiple_ensembl == F & any_alt5utr == F & N_TfTargets_in_transnet >10 & umis_per_cell_norm*10000 >0.1 & !(gene %in% titration_genes), .N, cluster]
Uniq_NFE2_features_symbols[multiple_ensembl == F & any_alt5utr == F & N_TfTargets_in_transnet >10 & umis_per_cell_norm*10000 >0.1 & !(gene %in% titration_genes), .N, cluster]

# Modify numbers depending on absence of TFs to be selected
N_clust[network == "NFE2" & cluster %in% 3:2, n_TFs := 0]
N_clust[network == "GFI1B" & cluster == 2, n_TFs := 3]

# The rest of the selected genes should be target or not of these TFs (50% representation of each category)
N_clust[, n_TFtargets := ceiling( (adj_n_target - n_TFs)/2 )]
N_clust[, n_others := floor( (adj_n_target - n_TFs)/2 )]

# Select TF target genes
set.seed(1234)
Selected_TFs <- lapply(c("GFI1B", "NFE2"), function(net){
  do.call("c", lapply(1:4, function(c){
    tfs = Net_L[[net]][N_TfTargets_in_transnet >10 & umis_per_cell_norm*10000 >0.1 & !(gene %in% titration_genes) & cluster == c , gene]
    sample(tfs, size = N_clust[ network == net & cluster == c, n_TFs])
  }))
})
names(Selected_TFs) <- c("GFI1B", "NFE2")

# Once selected the TFs, find which genes in the network are TF targets of those
Targets_selected_TFs <- lapply(c("GFI1B", "NFE2"), function(net){
  unlist(lapply(TFS[TfGene %in% Selected_TFs[[net]], TfTargets_in_transnet ], strsplit, "__"))
})
names(Targets_selected_TFs) <- c("GFI1B", "NFE2")
Targets_shared_selected_TFs <- unlist(lapply(Selected_shared_genes[ N_TfTargets_in_transnet > 0, TfTargets_in_transnet], strsplit, "__"))

# Count the number of times those genes are targets of any TF found in both transnetworks
Net_L[["GFI1B"]][, Target_of_selected_TF := sum(c(Targets_selected_TFs[["GFI1B"]], Targets_selected_TFs[["NFE2"]], Targets_shared_selected_TFs) %in% gene), gene]
Net_L[["NFE2"]][, Target_of_selected_TF := sum(c(Targets_selected_TFs[["GFI1B"]], Targets_selected_TFs[["NFE2"]], Targets_shared_selected_TFs) %in% gene), gene]


# Make sure in each cluster there are enough genes from both categories (target or not of a TF)
# TF targets
Net_L[["GFI1B"]][ is.na(N_TfTargets) & mean_abs_cor > 0.2 & umis_per_cell_norm*10000 >0.1 & !(gene %in% titration_genes) & Target_of_selected_TF > 0, .N, cluster]
Net_L[["NFE2"]][ is.na(N_TfTargets) & mean_abs_cor > 0.2 & umis_per_cell_norm*10000 >0.1 & !(gene %in% titration_genes) & Target_of_selected_TF > 0, .N, cluster]
# Others
Net_L[["GFI1B"]][ is.na(N_TfTargets) & mean_abs_cor > 0.2 & umis_per_cell_norm*10000 >0.1 & !(gene %in% titration_genes) & Target_of_selected_TF == 0, .N, cluster]
Net_L[["NFE2"]][ is.na(N_TfTargets) & mean_abs_cor > 0.2 & umis_per_cell_norm*10000 >0.1 & !(gene %in% titration_genes) & Target_of_selected_TF == 0, .N, cluster]

# Change to have enough genes to chose from
N_clust[network == "NFE2" & cluster == 1, n_TFtargets := 2]
N_clust[network == "NFE2" & cluster == 1, n_others := 0 ]


# Select targets of the selected TFs
set.seed(3412)
Selected_targets <- lapply(c("GFI1B", "NFE2"), function(net){
  do.call("c", lapply(1:4, function(c){
    targets = Net_L[[net]][ is.na(N_TfTargets) & mean_abs_cor > 0.2 & umis_per_cell_norm*10000 >0.1 & !(gene %in% titration_genes) & cluster == c & Target_of_selected_TF > 0, gene]
    targets = targets[!(grepl("^(MT-)|(LINC)|(MRP)", targets))]
    sample(targets, size = N_clust[ network == net & cluster == c, n_TFtargets])
  }))
})
names(Selected_targets) <- c("GFI1B", "NFE2") 

hist(Net_L[["GFI1B"]][gene %in% Selected_targets[["GFI1B"]], Target_of_selected_TF], breaks = 7)
hist(Net_L[["NFE2"]][gene %in% Selected_targets[["NFE2"]], Target_of_selected_TF], breaks = 7)


# Select the rest of the target genes
set.seed(2341)
Selected_others <- lapply(c("GFI1B", "NFE2"), function(net){
  do.call("c", lapply(1:4, function(c){
    others = Net_L[[net]][ is.na(N_TfTargets) & mean_abs_cor > 0.2 & umis_per_cell_norm*10000 >0.1 & !(gene %in% c("GFI1B", "NFE2")) & cluster == c & Target_of_selected_TF == 0, gene]
    others = others[!(grepl("^(MT-)|(LINC)|(MRP)|(RP)", others))]
    sample(others, size = N_clust[ network == net & cluster == c, n_others])
  }))
})
names(Selected_others) <- c("GFI1B", "NFE2") 

GFI1B_unique_selected_genes <- c(Selected_TFs[["GFI1B"]], Selected_targets[["GFI1B"]], Selected_others[["GFI1B"]])
Selected_GFI1B_genes <- Net_L[["GFI1B"]][gene %in% GFI1B_unique_selected_genes]

NFE2_unique_selected_genes <- c(Selected_TFs[["NFE2"]], Selected_targets[["NFE2"]], Selected_others[["NFE2"]])
Selected_NFE2_genes <- Net_L[["NFE2"]][gene %in% NFE2_unique_selected_genes]

Selected_shared_genes[, Target_of_selected_TF := sum(c(Targets_selected_TFs[["GFI1B"]], Targets_selected_TFs[["NFE2"]], Targets_shared_selected_TFs) %in% gene), gene]
Shared_unique_selected_genes <- Selected_shared_genes$gene


Selected_genes_total <- unique(c(GFI1B_unique_selected_genes, NFE2_unique_selected_genes, Shared_unique_selected_genes, control_genes, titration_genes))


X <- as.data.table(biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = Selected_genes_total,
  mart = ensembl
))
#fwrite(X, file = "processed_data/20211111_list_genes_targeted_seq_test.txt", quote = F, row.names = F)
