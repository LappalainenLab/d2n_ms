##### D2N - OMIM and GWAS genes enrichments for non-linear dosage responses #####
# (1) OMIM disease genes
# (2) Blood traits GWAS genes


# JDE, February 2023
# Last modified: August 2023


## Libs
library(data.table)
library(ggplot2); theme_set(theme_bw());
library(foreach)
library(doParallel); registerDoParallel(detectCores()/2)
library(cowplot)


## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/10-complex_traits_and_disease")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


## Gene annotations
GA.1 <- readRDS("../../processed_data/08-gene_dosage_properties/GA.1_GWASgenesUKB.RDS")
GA.2 <- readRDS("../../processed_data/08-gene_dosage_properties/GA.2_DiseaseGenesOMIM.RDS")

# Get model parameters and fitting stats
RMSE_dt <- fread(file = "../../processed_data/07-nonlinear_fit/D2N_AIC_LmSig4Loess.txt")
S4Param_dt <- fread(file = "../../processed_data/07-nonlinear_fit/D2N_S4Param_10fCV.txt")


RMSE_S4P_dt <- merge.data.table(RMSE_dt[, .(dosage_gene, gene, loess_rmse, lm_rmse, delta_RMSE, delta_AIC, delta_RMSE_sig_loess)], 
                                S4Param_dt, 
                                by=c("gene", "dosage_gene"))
RMSE_S4P_dt[, non_linear := ifelse(delta_AIC >= 0, T, F)]


### (1) Are OMIM genes enriched in non-linear dosage-responses
omim_genes <- GA.2[value == T, gene]

RMSE_S4P_dt[, omim_gene := ifelse(gene %in% omim_genes, T, F)]

Odds_plot_dt <- foreach(dg = c("GFI1B", "NFE2", "MYB"), .combine = rbind) %do% {
  test = fisher.test(table(RMSE_S4P_dt[dosage_gene == dg & unresponsive == F, .(omim_gene, non_linear)]))
  data.table(odds = log2(test$estimate), 
             ci_down = log2(test$conf.int[1]), 
             ci_up = log2(test$conf.int[2]), 
             pval = test$p.value,
             signif = ifelse(test$p.value < 0.05, T, F),
             dosage_gene = dg)
}
pA <- ggplot(Odds_plot_dt, aes(x=odds, y=dosage_gene, color = signif)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 3, show.legend = F) +
  geom_errorbarh(aes(xmin = ci_down, xmax = ci_up), height = 0.2, show.legend = F) +
  scale_color_manual("OMIM gene?", values = c(`TRUE` = "#56B1F7", `FALSE`= "black")) +
  ylab("Trans network") + 
  theme(axis.title.x = element_blank()) +
  ggtitle("OMIM genes") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 


### (2) Are GWAS genes enriched in non-linear dosage-responses
gwas_genes <- GA.1[value == T, gene]

RMSE_S4P_dt[, gwas_gene := ifelse(gene %in% gwas_genes, T, F)]
RMSE_S4P_dt[, Platelets := gene %in% GA.1[gene_set == "Platelets(UKB GWAS)" & value == T, gene]]
RMSE_S4P_dt[, RBCs := gene %in% GA.1[gene_set == "RBCs(UKB GWAS)" & value == T, gene]]
RMSE_S4P_dt[, WBCs := gene %in% GA.1[gene_set == "WBCs(UKB GWAS)" & value == T, gene]]
RMSE_S4P_dt[, Reticulocytes := gene %in% GA.1[gene_set == "Reticulocytes(UKB GWAS)" & value == T, gene]]

gwas_traits_v <- c("Platelets",  "RBCs",  "WBCs", "Reticulocytes")

RMSE_S4P_GWASmelt_dt <- melt.data.table(RMSE_S4P_dt[, .SD , .SDcols = c("delta_AIC", "dosage_gene", "unresponsive", "gene", "non_linear", gwas_traits_v)], 
                                        id.vars = c("delta_AIC", "dosage_gene", "unresponsive", "gene", "non_linear"),
                                        variable.name = "GWAS_trait", value.name = "GAWAS_gene")


Odds_plot_GWAS_dt <- foreach(dg = c("GFI1B", "NFE2", "MYB"), .combine = rbind) %do% {
  foreach(gt = gwas_traits_v, .combine = rbind) %do% {
    test = fisher.test(table(RMSE_S4P_GWASmelt_dt[dosage_gene == dg & unresponsive == F & GWAS_trait == gt, .(GAWAS_gene, non_linear)]))
    data.table(odds = log2(test$estimate), 
               ci_down = log2(test$conf.int[1]), 
               ci_up = log2(test$conf.int[2]), 
               pval = test$p.value,
               signif = ifelse(test$p.value < 0.05, T, F),
               dosage_gene = dg,
               GWAS_trait = gt)
  }
}
Odds_plot_GWAS_dt[, pval_fdr := p.adjust(pval, method = "fdr")]


pB <- ggplot(Odds_plot_GWAS_dt, aes(x=odds, y=dosage_gene, color = signif)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 3, show.legend = F) +
  geom_errorbarh(aes(xmin = ci_down, xmax = ci_up), height = 0.2, show.legend = F) +
  scale_color_manual("OMIM gene?", values = c(`TRUE` = "#56B1F7", `FALSE`= "black")) +
  xlab("log2(Odds) for non-linearities") + ylab("Trans network") +
  ggtitle("GWAS genes") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  facet_grid(GWAS_trait ~ .)


p <- cowplot::plot_grid(pA, pB, align = "v", axis = "lr", nrow = 2, rel_heights = c(1/4, 3/4))
p
ggsave(file.path(plots_dir, "10a_OMIM_GWAS_Nonlinear_Odds.pdf"), width = 2.5 , height = 5)



## (3) Examples dosage-response curves of different disease gwas genes across transnetwork

# Expression data
DE_dt <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS")

Fit_curves_dt <- fread("../../processed_data/07-nonlinear_fit/D2N_S4LoessModelFits.txt")

union_genes <- omim_genes[omim_genes %in% gwas_genes]

dt_unresp <- S4Param_dt[, sum(unresponsive), gene]

plot_genes <- dt_unresp[gene %in% union_genes & V1 < 2, gene]

p <- ggplot(DE_dt[gene %in% plot_genes & dosage_gene != "TET2" & dosage_gene != "NTC", ], aes(x = dosage_gene_log2FC, y=avg_log2FC)) +
  facet_grid(gene ~ dosage_gene, scales = "free") +
  geom_point(color = "grey20", alpha=0.7) +
  geom_line(data = Fit_curves_dt[gene %in% plot_genes & dosage_gene != "TET2", ], color="#FF7F00", linewidth = 0.75) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  labs(x="Cis gene log2(FC)", y="Trans gene log2(FC)")
p
ggsave(filename = file.path(plots_dir, "10b_OMIM_GWAS_DosageResponseCurves.pdf"), p, width = 5, height = 13)

