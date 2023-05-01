##### D2N - OMIM and GWAS genes enrichments for non-linear dosage responses #####
# (1) OMIM disease genes
# (2) Blood traits GWAS genes


# JDE, February 2023
# Last modified: April 2023


## Libs
library(data.table)
library(ggplot2); theme_set(theme_bw());
library(foreach)
library(doParallel); registerDoParallel(detectCores()/2)
library(cowplot)


## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/08-gene_dosage_properties/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


## Gene annotations
GA.9 <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/016-gene_annotations/processed_data/GA.9_DiseaseGenesOMIM.RDS")
GA.10 <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/016-gene_annotations/processed_data/GA.10_GWASgenesUKB.RDS")

RMSE_dt <- fread(file = "../../processed_data/07-nonlinear_fit/D2N_AIC_LmSig4Loess.txt")
S4Param_dt <- fread(file = "../../processed_data/07-nonlinear_fit/D2N_S4Param_10fCV.txt")

### (1) Are OMIM genes enriched in non-linear dosage-responses
omim_genes <- GA.9[value == T, gene]



RMSE_S4P_dt <- merge.data.table(RMSE_dt[, .(dosage_gene, gene, loess_rmse, lm_rmse, delta_RMSE, delta_AIC, delta_RMSE_sig_loess)], 
                                S4Param_dt, 
                                by=c("gene", "dosage_gene"))

RMSE_S4P_dt[, omim_gene := ifelse(gene %in% omim_genes, T, F)]
RMSE_S4P_dt[, non_linear := ifelse(delta_AIC >= 0, T, F)]

Odds_plot_dt <- foreach(dg = c("GFI1B", "NFE2", "MYB"), .combine = rbind) %do% {
  test = fisher.test(table(RMSE_S4P_dt[dosage_gene == dg & unresponsive == F, .(omim_gene, non_linear)]))
  data.table(odds = log2(test$estimate), 
             ci_down = log2(test$conf.int[1]), 
             ci_up = log2(test$conf.int[2]), 
             pval = test$p.value,
             signif = ifelse(test$p.value < 0.05, T, F),
             dosage_gene = dg)
}
ggplot(Odds_plot_dt, aes(x=odds, y=dosage_gene, color = signif)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 3, show.legend = F) +
  geom_errorbarh(aes(xmin = ci_down, xmax = ci_up), height = 0.2, show.legend = F) +
  scale_color_manual("OMIM gene?", values = c(`TRUE` = "#56B1F7", `FALSE`= "black")) +
  xlab("log2(Odds)") + ylab("") +
  ggtitle("OMIM genes having non-linear dosage responses")

