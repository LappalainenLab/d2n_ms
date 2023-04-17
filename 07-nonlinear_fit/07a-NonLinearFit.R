##### D2N - Non-linear fit #####
# (1) Fit loess, linear and sigmoidal model  and measure RMSE and AIC from all models per unique dosage response
#     - Plot Lm vs. Simoid AIC & delta AIC
#     - Plot Lm vs. Sigmoid AIC & delta AIC in CNV range (Supplementary plot)
#     - Plot Loess vs. Sigmoid RMSE (Supplmentary plot)
# (2) Fit linear and sigmoidal model in 10f-CV scheme
#     - Plot Pred vs. expected with R^2 (Supplmentary plot)
# (3) Identify:
#     - Non-responsive genes (Supplementary plot)
#     - Non-monotonic genes (Supplementary plot)
# (4) Plot all dosage reponses with corresponding fit (loess for non-monotonic and sigmoid for rest)
# (5) Plot subsets of genes with all cis genes in a single plot
#     - GWAS and OMIM genes


# JDE, December 2022
# Last modified: April 2023


## Libs
library(ggplot2); theme_set(theme_bw())
library(data.table)
library(foreach)
library(doParallel); registerDoParallel(detectCores()/2)
library(gridExtra)
library(vegan)
library(drc)
library(ggrepel)


## Dir
setwd("/gpfs/commons/groups/lappalainen_lab/jdomingo/manuscripts/d2n/d2n_ms/07-nonlinear_fit/")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


## Color palette
dosage_genes = c("GFI1B", "MYB", "NFE2", "TET2")
col_genes = RColorBrewer::brewer.pal(5, "Set1")
names(col_genes) <- c("MYB", "TET2", "NFE2", "GFI1B", "NTC")
col_crispr = list(CRISPRa = "#FF7F00", CRISPRi = "#377EB8")
rect_trans_GFI1B <- c("FOXP1", "GATA2", "PITX1", "RHD", "MYB", "KLK1", "TUBB1", "SLC22A4", "LHX3")


## Functions
RMSE <- function(residuals) {
  sqrt(mean((residuals)^2))
}


## Data
DE_dt_raw <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS") 
DE_dt_raw[, guide_crispr := paste0(guide_1, "-", cell_line)][]
DE_dt <- DE_dt_raw[!(guide_crispr == "NFE2_9-CRISPRa" | dosage_gene %in% c("NTC") | grepl("CRISPR", gene))]



### (1) Fit loess, linear and sigmoidal model 

# Genes vector
d2n_genes = unique(DE_dt$gene)
trans_genes = d2n_genes[!(d2n_genes %in% dosage_genes)]

# Fit sigmoid, linear and Loess models and calculate AIC for Lm and Sig, or RMSE for each
RMSE_dt <- foreach(dg = dosage_genes, .combine = rbind) %do% {
  trans_genes_dg = d2n_genes[d2n_genes != dg]
  foreach(tg = trans_genes_dg, .combine = rbind) %dopar% {
    m.loess <- loess(data = DE_dt[dosage_gene == dg & gene == tg, .(avg_log2FC, dosage_gene_log2FC)], formula = 'avg_log2FC ~ dosage_gene_log2FC')
    m.lm <- lm(data = DE_dt[dosage_gene == dg & gene == tg, .(avg_log2FC, dosage_gene_log2FC)], formula = 'avg_log2FC ~ dosage_gene_log2FC')
    m.sig <- drm(avg_log2FC ~ dosage_gene_log2FC, fct = L.4(), data = DE_dt[dosage_gene == dg & gene == tg, .(dosage_gene_log2FC, avg_log2FC)])
    data.table(gene = tg,
               loess_rmse = RMSE(m.loess$residuals),
               lm_rmse = RMSE(m.lm$residuals),
               sig_rmse = RMSE(residuals(m.sig)),
               lm_aic = stats::AIC(m.lm),
               sig_aic = stats::AIC(m.sig),
               dosage_gene = dg,
               lm_slope = m.lm$coefficients[2],
               sig_slope = (-1)*m.sig$coefficients[1],
               min_max_ast = m.sig$coefficients[3]-m.sig$coefficients[2])
  }
}
RMSE_dt[, delta_RMSE := lm_rmse - sig_rmse]
RMSE_dt[, delta_AIC := lm_aic - sig_aic]
RMSE_dt[, delta_RMSE_sig_loess := sig_rmse - loess_rmse]

mean_AIC <- RMSE_dt[!(gene %in% c("LHX3", "GAPDH")), .(mean_delta_AIC = mean(delta_AIC)), dosage_gene]
pA <- ggplot(RMSE_dt[!(gene %in% c("LHX3", "GAPDH")),], aes(x=sig_aic, y=lm_aic)) +
  facet_grid(. ~ dosage_gene) +
  geom_abline(color="grey50", linetype=2) +
  geom_point(size=2, alpha=0.5) +
  labs(x = "AIC (Sigmoidal fit)", y = "AIC (Linear fit)") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
pB <- ggplot(RMSE_dt, aes(x = delta_AIC)) +
  geom_histogram(bins = 40) +
  facet_grid(. ~ dosage_gene) +
  labs(x=paste0(expression(delta), " AIC (lm - sigmoid)"))  +
  geom_vline(xintercept = 0, color="grey50") +
  geom_vline(data = mean_AIC, aes(xintercept = mean_delta_AIC), color = "#E41A1C", linetype = 2) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(panel.grid.minor = element_blank())
p <- cowplot::plot_grid(pB, pA, align = "v", axis = "tb", nrow = 2, rel_heights = c(0.3, 0.7))
p

ggsave(file.path(plots_dir, "07a_01_SigVsLM_AIC.pdf"), p, width = 10, height = 4.5)

## !!! Finish Just for GFI1B ####






### (3) Fit linear and sigmoidal model in 10f-CV scheme
seed_dt <- data.table(DG = do.call("c", lapply(dosage_genes, rep, length(d2n_genes))),
                      TG = rep(d2n_genes, length(dosage_genes)))
seed_dt[, s_num := 1:.N ]
n_cv = 10

S4Param_10fCV_dt <- foreach(dg = dosage_genes, .combine = rbind) %do% {
  foreach(tg = d2n_genes, .combine = rbind) %dopar% {
    dt <- DE_dt[dosage_gene == dg & gene == tg, .(avg_log2FC, dosage_gene_log2FC, guide_crispr)]
    set.seed(seed_dt[DG == dg & TG == tg, s_num])
    test_dt_l <- split(dt, cut_number(sample(1:nrow(dt)), n_cv))
    foreach(cv = 1:length(test_dt_l), .combine = rbind) %do% {
      test_perturb = test_dt_l[[cv]][, guide_crispr]
      test_data = dt[guide_crispr %in% test_perturb, .(dosage_gene_log2FC, avg_log2FC, guide_crispr)]
      train_pertubr = dt$guide_crispr[!(dt$guide_crispr) %in% test_perturb]
      length(train_pertubr)
      m.sig <- drm(avg_log2FC ~ dosage_gene_log2FC, fct = L.4(), data = dt[guide_crispr %in% train_pertubr, .(dosage_gene_log2FC, avg_log2FC)])
      data.table(slope_IF = (-1)*m.sig$coefficients[1],
                 min_assmp = m.sig$coefficients[2],
                 max_assmp = m.sig$coefficients[3],
                 x_IF = m.sig$coefficients[4],
                 gene = tg,
                 dosage_gene = dg,
                 cv_fold = cv)
    }
  }
}
S4Param_10fCV_dt[, min_max_range := max_assmp - min_assmp]

# Collapse and calc mean and sd
S4Param_dt <- S4Param_10fCV_dt[gene != dosage_gene, .(slope_IF = mean(slope_IF), 
                                                      slope_IF_sd = sd(slope_IF),
                                                      min_assmp = mean(min_assmp),
                                                      min_assmp_sd = sd(min_assmp),
                                                      max_assmp = mean(max_assmp),
                                                      max_assmp_sd = sd(max_assmp),
                                                      x_IF = mean(x_IF),
                                                      x_IF_sd = sd(x_IF),
                                                      min_max_range = mean(min_max_range),
                                                      min_max_range_sd = sd(min_max_range)), 
                               .(gene, dosage_gene)]

# Compare it to parameters when using all data to fit models
PT_singles_S4Param[, min_max_range := max_assmp - min_assmp]








#### (4)

trans_genes_ord <- trans_genes[order(trans_genes)]