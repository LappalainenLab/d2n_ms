##### D2N - Non-linear fit #####
# (1) Fit loess, linear and sigmoidal model  and measure RMSE and AIC from all models per unique dosage response
#     - Plot Lm vs. Simoid AIC & delta AIC
#     - Plot Lm vs. Sigmoid AIC & delta AIC in CNV range (Supplementary plot)
# (2) Fit linear and sigmoidal model in 10f-CV scheme
#     - Plot Pred vs. expected with R^2 (Supplmentary plot)
# (3) Identify:
#     - Non-responsive genes (Supplementary plot)
#     - Non-monotonic genes (Supplementary plot)
#     - Plot Loess vs. Sigmoid RMSE (Supplmentary plot)
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
rect_trans_GFI1B <- c("TUBB1", "RHD", "KLK1", "SLC22A4", "PITX1", "LHX3", "GATA2", "FOXP1", "MYB")
control_genes = c("LHX3", "GAPDH")

## Functions
RMSE <- function(residuals) {
  sqrt(mean((residuals)^2))
}

ztest <- function(av, se, mu, alternative = "two.sided") {
  zscore <- (av - mu)/se
  if (alternative == "two.sided") {
    pval <- 2*pnorm(abs(zscore), lower.tail = FALSE)
  } else if (alternative == "less") {
    pval <- pnorm(zscore, lower.tail = TRUE)
  } else if (alternative == "greater") {
    pval <- pnorm(zscore, lower.tail = FALSE)
  }
  return(pval)
}

SigmoidL4Predict <- function(x, b, c, d, e) {
  #x Input value
  #Slope parameter
  #Lower asymptote
  #Upper asymptote
  #Inflection point
  y = c + ((d-c)/(1+exp(-b*(x-e))))
  return(y)
}


## Data
DE_dt_raw <- readRDS("/gpfs/commons/groups/lappalainen_lab/jdomingo/projects/004-dosage_network/009-DEanalyses/processed_data/after_seurat4.3/d2n_DemuxHTOdCas9_Wilcoxon_AllGenesDE.RDS") 
DE_dt_raw[, guide_crispr := paste0(guide_1, "-", cell_line)][]
# Remove TET2 perturbations, the NFE2 outlier and the CRISPR genes
DE_dt <- merge.data.table(DE_dt_raw[!(guide_crispr == "NFE2_9-CRISPRa" | dosage_gene %in% c("NTC") | grepl("CRISPR", gene))], 
                          DE_dt_raw[gene == dosage_gene, .(avg_log2FC_dg = avg_log2FC), guide_crispr], 
                          by = "guide_crispr")


### (1) Fit loess, linear and sigmoidal model 

# Genes vector
d2n_genes = unique(DE_dt$gene)
trans_genes = d2n_genes[!(d2n_genes %in% dosage_genes)]

# Fit sigmoid, linear and Loess models and calculate AIC for Lm and Sig, or RMSE for each (Using full cis gene FC range)
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
               dosage_gene = dg)
  }
}
RMSE_dt[, delta_RMSE := lm_rmse - sig_rmse]
RMSE_dt[, delta_AIC := lm_aic - sig_aic]
RMSE_dt[, delta_RMSE_sig_loess := sig_rmse - loess_rmse]
RMSE_dt[, cis_trans := paste0(dosage_gene, "_", gene)]

# Plot AIC of Lm vs. Sigmoid model
mean_AIC <- RMSE_dt[!(gene %in% control_genes), .(mean_delta_AIC = mean(delta_AIC)), dosage_gene]
pA <- ggplot(RMSE_dt[!(gene %in% control_genes),], aes(x=sig_aic, y=lm_aic)) +
  facet_grid(. ~ dosage_gene) +
  geom_abline(color="grey50", linetype=2) +
  geom_point(size=2, alpha=0.5) +
  labs(x = "AIC (Sigmoidal fit)", y = "AIC (Linear fit)") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
pB <- ggplot(RMSE_dt[!(gene %in% control_genes),], aes(x = delta_AIC)) +
  geom_histogram(bins = 40) +
  facet_grid(. ~ dosage_gene) +
  labs(x=paste0(expression(delta), " AIC (lm - sigmoid)"))  +
  geom_vline(xintercept = 0, color="grey50") +
  geom_vline(data = mean_AIC, aes(xintercept = mean_delta_AIC), color = "#E41A1C", linetype = 2) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(panel.grid.minor = element_blank())
p <- cowplot::plot_grid(pB, pA, align = "v", axis = "tb", nrow = 2, rel_heights = c(0.3, 0.7))
p
ggsave(file.path(plots_dir, "07a_01_SigVsLM_AIC.pdf"), p, width = 9, height = 4)

# Plot AIC of Lm vs. Sigmoid model (only for GFI1B)
pA <- ggplot(RMSE_dt[!(gene %in% control_genes) & dosage_gene == "GFI1B",], aes(x=sig_aic, y=lm_aic)) +
  facet_grid(. ~ dosage_gene) +
  geom_abline(color="grey50", linetype=2) +
  geom_point(size=2, alpha=0.5) +
  labs(x = "AIC (Sigmoidal fit)", y = "AIC (Linear fit)") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
pB <- ggplot(RMSE_dt[!(gene %in% control_genes) & dosage_gene == "GFI1B",], aes(x = delta_AIC)) +
  geom_histogram(bins = 40) +
  facet_grid(. ~ dosage_gene) +
  labs(x=paste0(expression(delta), " AIC (linear - sigmoid)"))  +
  geom_vline(xintercept = 0, color="grey50") +
  geom_vline(data = mean_AIC[dosage_gene == "GFI1B", ], aes(xintercept = mean_delta_AIC), color = "#E41A1C", linetype = 2) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(panel.grid.minor = element_blank())
p <- cowplot::plot_grid(pB, pA, align = "v", axis = "tb", nrow = 2, rel_heights = c(0.3, 0.7))
p
ggsave(file.path(plots_dir, "07a_02_SigVsLM_AIC_GFI1B.pdf"), p, width = 2.75, height = 3.75)


# Limit cis gene FC range to the one gene copy gain or loss
CN_range <- c(log2(0.5), log2(3/2))

# Fit again the Lm, Loess and Sigmoid models and calculate AIC 
RMSE_CNrange <- foreach(dg = dosage_genes, .combine = rbind) %do% {
  
  # Get points from the cis gene within the FC range defined above
  CIS_dt <- DE_dt[dosage_gene == dg & gene == dg & (dosage_gene_log2FC >= CN_range[1] & dosage_gene_log2FC <= CN_range[2]), ]
  order_dosage <- CIS_dt[order(avg_log2FC), guide_crispr]
  CIS_dt[, guide_crispr := factor(guide_crispr, levels = order_dosage)]
  
  # Only get perturbations within the range
  G_DG <- DE_dt[dosage_gene == dg & guide_crispr %in% order_dosage, ]
  
  
  # Define min and max FC in for each dosage gene, and create a range of new cis FC values to get predicted trans FCs
  min_d = min(CIS_dt$avg_log2FC)
  max_d = max(CIS_dt$avg_log2FC)
  
  # For each trans gene, predict new values given the exact cis FC or equally binned FC within the min-max range, generate new data.table with smoothen values
  DG_DE_all <- foreach(tg = trans_genes, .combine = rbind) %dopar% {
    
    # Loess fitting
    m.sig <- drm(avg_log2FC ~ dosage_gene_log2FC, fct = L.4(), data = G_DG[gene == tg, .(dosage_gene_log2FC, avg_log2FC)])
    m.lm <- lm(data = G_DG[gene == tg, .(avg_log2FC, dosage_gene_log2FC)], formula = 'avg_log2FC ~ dosage_gene_log2FC')
    
    data.table(gene = tg,
               lm_aic = stats::AIC(m.lm),
               sig_aic = stats::AIC(m.sig),
               dosage_gene = dg)
  }
}
RMSE_CNrange[, delta_AIC := lm_aic - sig_aic]

# Repeat same scatter AIC plot and distribution of Delta AIC
mean_AIC <- RMSE_CNrange[!(gene %in% control_genes), .(mean_delta_AIC = mean(delta_AIC)), dosage_gene]
pA <- ggplot(RMSE_CNrange[!(gene %in% control_genes),], aes(x=sig_aic, y=lm_aic)) +
  facet_grid(. ~ dosage_gene) +
  geom_abline(color="grey50", linetype=2) +
  geom_point(size=2, alpha=0.5) +
  labs(x = "AIC (Sigmoidal fit)", y = "AIC (Linear fit)") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
pB <- ggplot(RMSE_CNrange[!(gene %in% control_genes),], aes(x = delta_AIC)) +
  geom_histogram(bins = 40) +
  facet_grid(. ~ dosage_gene) +
  labs(x=paste0(expression(delta), " AIC (lm - sigmoid)"))  +
  geom_vline(xintercept = 0, color="grey50") +
  geom_vline(data = mean_AIC, aes(xintercept = mean_delta_AIC), color = "#E41A1C", linetype = 2) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(panel.grid.minor = element_blank())
p <- cowplot::plot_grid(pB, pA, align = "v", axis = "tb", nrow = 2, rel_heights = c(0.3, 0.7))
p
ggsave(file.path(plots_dir, "07a_03_SigVsLM_AIC_CNVrange.pdf"), p, width = 9, height = 4)




### (3) Fit linear and sigmoidal model in 10f-CV scheme
seed_dt <- data.table(DG = do.call("c", lapply(dosage_genes, rep, length(d2n_genes))),
                      TG = rep(d2n_genes, length(dosage_genes)))
seed_dt[, s_num := 1:.N ]
n_cv = 10

# Get Predicted data when using 10-fold CV
SigPred10CV_dt <- foreach(dg = dosage_genes, .combine = rbind) %do% {
  foreach(tg = d2n_genes, .combine = rbind) %dopar% {
    dt <- DE_dt[dosage_gene == dg & gene == tg, .(avg_log2FC, dosage_gene_log2FC, guide_crispr)]
    set.seed(seed_dt[DG == dg & TG == tg, s_num])
    test_dt_l <- split(dt, cut_number(sample(1:nrow(dt)), n_cv))
    foreach(cv = 1:length(test_dt_l), .combine = rbind) %do% {
      test_perturb = test_dt_l[[cv]][, guide_crispr]
      test_data = dt[guide_crispr %in% test_perturb, .(dosage_gene_log2FC, avg_log2FC, guide_crispr)]
      train_pertubr = dt$guide_crispr[!(dt$guide_crispr) %in% test_perturb]
      m.sig <- drm(avg_log2FC ~ dosage_gene_log2FC, fct = L.4(), data = dt[guide_crispr %in% train_pertubr, .(dosage_gene_log2FC, avg_log2FC)])
      data.table( gene = tg,
                  dosage_gene = dg,
                  guide_crispr = test_data$guide_crispr,
                  pred_avg_log2 = predict(m.sig, newdata = test_data),
                  cv_fold = cv)
    }
  }
}

# Plot predicted vs. observed
ObsPred10CV_dt <- merge.data.table(SigPred10CV_dt, DE_dt[, .(gene, dosage_gene, guide_crispr, avg_log2FC)], by = c("gene", "dosage_gene", "guide_crispr"))

cor_plot <- unique(ObsPred10CV_dt[, r := cor(pred_avg_log2, avg_log2FC), dosage_gene][, .(dosage_gene, r)])
p <- ggplot(ObsPred10CV_dt, aes(y = avg_log2FC, x = pred_avg_log2)) +
  stat_bin2d(bins = 50) +
  facet_grid(. ~ dosage_gene) +
  scale_fill_gradient(low = "grey85", high = "grey30", guide = guide_colorbar(direction = "vertical", barwidth = 0.5))  +
  geom_abline(linetype = 2)  +
  labs(x = "Predicted log2(FC)", y = "Observed log2(FC)") +
  geom_text(aes(x= -1, y=1.25, label = paste0("r = ", round(r, 2))), data = cor_plot) +
  coord_cartesian(xlim = c(-2, 1.2)) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(legend.position = "bottom")
p
ggsave(file.path(plots_dir, "07a_04_10fCV_PredVsObs.pdf"), p, width = 9, height = 4.5)

# Get sigmoid parameters for each CV fold
S4Param_10fCV_dt <- foreach(dg = dosage_genes, .combine = rbind) %do% {
  foreach(tg = d2n_genes, .combine = rbind) %dopar% {
    dt <- DE_dt[dosage_gene == dg & gene == tg, .(avg_log2FC, dosage_gene_log2FC, guide_crispr)]
    set.seed(seed_dt[DG == dg & TG == tg, s_num])
    test_dt_l <- split(dt, cut_number(sample(1:nrow(dt)), n_cv))
    foreach(cv = 1:length(test_dt_l), .combine = rbind) %do% {
      test_perturb = test_dt_l[[cv]][, guide_crispr]
      test_data = dt[guide_crispr %in% test_perturb, .(dosage_gene_log2FC, avg_log2FC, guide_crispr)]
      train_pertubr = dt$guide_crispr[!(dt$guide_crispr) %in% test_perturb]
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

# Collapse and calculate mean and s.d.
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
S4Param_dt[, cis_trans := paste0(dosage_gene, "_", gene)]


#### (3) Classify genes into responsive, unresponsive and non-monotonic

# Unresponsive
S4Param_dt[, fdr_slope := p.adjust(S4Param_dt[, ztest(slope_IF, slope_IF_sd/sqrt(n_cv), 0)])]
S4Param_dt[, fdr_range := p.adjust(S4Param_dt[, ztest(min_max_range, min_max_range_sd/sqrt(n_cv), 0.05, "greater")])]
S4Param_dt[, unresponsive :=  ifelse(fdr_slope > 0.01 | fdr_range > 0.01, T, F)]
S4Param_dt[, .N, .(dosage_gene, unresponsive)]
#    dosage_gene unresponsive  N
# 1:       GFI1B         TRUE 18
# 2:       GFI1B        FALSE 73
# 3:         MYB        FALSE 57
# 4:         MYB         TRUE 34
# 5:        NFE2        FALSE 49
# 6:        NFE2         TRUE 42
# 7:        TET2        FALSE 21
# 8:        TET2         TRUE 70


# Non-monotonic (for GFI1B since has the largest range)
RMSE_G_dt <- RMSE_dt[dosage_gene == "GFI1B" ,]
nonmono_cistrans <- RMSE_G_dt[delta_RMSE_sig_loess >= quantile(delta_RMSE_sig_loess, 0.95), cis_trans]
mono_cistrans <- RMSE_dt$cis_trans[!(RMSE_dt$cis_trans %in% nonmono_cistrans)]
S4Param_dt[, non_monotonic := ifelse(cis_trans %in% nonmono_cistrans, T, F)]




RMSE_dt[, non_monotonic := ifelse(cis_trans %in% nonmono_cistrans, T, F)]

pA <- ggplot(RMSE_dt[!(gene %in% control_genes) & dosage_gene != "TET2", ], aes(x = delta_RMSE_sig_loess, color = non_monotonic, fill=non_monotonic)) +
  facet_grid(. ~ dosage_gene) +
  geom_histogram(bins = 40, alpha = 0.5) +
  scale_color_manual("Non-monotonic", values = c(`TRUE` = "#377EB8", `FALSE`= "black")) +
  scale_fill_manual("Non-monotonic", values = c(`TRUE` = "#377EB8", `FALSE`= "black")) +
  labs(x = "delta RMSE (Sigmoid - Loess)") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  theme(panel.grid.minor = element_blank())

pB <- ggplot(RMSE_dt[!(gene %in% control_genes) & dosage_gene != "TET2", ], aes(x=sig_rmse, y=loess_rmse)) +
  facet_grid(. ~ dosage_gene) +
  geom_abline(color="grey50", linetype=2) +
  geom_point(size=2, alpha=0.5, aes(color = non_monotonic)) +
  scale_color_manual("Non-monotonic", values = list(`TRUE` = "#377EB8", `FALSE` = "black")) +
  labs(x = "RMSE (Sigmoid)", y = "RMSE (Loess)") +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) 
p <- cowplot::plot_grid(pA, pB, nrow=2, align = "v", axis = "tb", rel_heights = c(0.3, 0.7))
p

ggsave(file.path(plots_dir, "07a_05_SigVsLoess_RMSE.pdf"), p, width = 8, height = 4)


#### (4) Plot all dosage reponses with corresponding fit (loess for non-monotonic and sigmoid for rest)

trans_genes_ord <- trans_genes[order(trans_genes)]

Fit_curves_dt <- foreach(dg = dosage_genes, .combine = rbind) %do% {
  
  cistrans_ord <- paste0(dg, "_", sort(d2n_genes[d2n_genes != dg]))
  dose_bins = seq(from = DE_dt[dosage_gene == dg, min(dosage_gene_log2FC)], to = DE_dt[dosage_gene == dg, max(dosage_gene_log2FC)], length.out = 100)
  
  foreach(ct = cistrans_ord, .combine = rbind) %dopar% {
    
    tg = unlist(strsplit(ct, "_"))[2]
    
    if (S4Param_dt[cis_trans == ct, non_monotonic]) {
      
      m.loess <- loess(data = DE_dt[gene == tg & dosage_gene == dg, .(avg_log2FC, dosage_gene_log2FC)], formula = 'avg_log2FC ~ dosage_gene_log2FC')
      
      newdata =  as.data.table(expand.grid(dosage_gene_log2FC = dose_bins))
      newdata$avg_log2FC <- predict(object = m.loess, newdata = dose_bins)
      newdata[, gene := tg]
      newdata[, dosage_gene := dg]
      
    } else {
      
      dt_params <- S4Param_dt[cis_trans == ct, ]
      
      newdata = as.data.table(expand.grid(dosage_gene_log2FC = dose_bins))
      newdata$avg_log2FC <- SigmoidL4Predict(x = newdata$dosage_gene_log2FC, dt_params$slope_IF, dt_params$min_assmp, dt_params$max_assmp, dt_params$x_IF)
      newdata[, gene := unlist(strsplit(ct, "_"))[2]]
      newdata[, dosage_gene := dg]
      
    }
  }
}


# Examples in heatmap GFI1B
plot_dt <- DE_dt[gene %in% rect_trans_GFI1B & dosage_gene == "GFI1B", ]
plot_dt[, gene_ord := factor(gene, levels = rect_trans_GFI1B)]

fit_dt <- Fit_curves_dt[gene %in% rect_trans_GFI1B & dosage_gene == "GFI1B", ]
fit_dt[, gene_ord := factor(gene, levels = rect_trans_GFI1B)]

p <- ggplot(plot_dt, aes(x = dosage_gene_log2FC, y=avg_log2FC)) +
  facet_grid(gene_ord ~ dosage_gene, ) +
  geom_point(color = "grey20", alpha=0.7) +
  geom_line(data = fit_dt, color="#FF7F00", linewidth = 0.75) +
  theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  labs(x="GFI1B log2(FC)", y="Trans gene log2(FC)")
p
ggsave(file.path(plots_dir, "07a_06_DosageRespCurves_ExamplesGFI1B.pdf"), p, width = 2.1, height = 12)


# Plot for each cis gene, all genes dosage responses observations and fits
foreach(dg = dosage_genes) %do% {
  ord_genes <- d2n_genes[d2n_genes != dg][order(d2n_genes[d2n_genes != dg])]
  plot_dt <- DE_dt[dosage_gene == dg & gene != dg, ]
  plot_dt[, gene := factor(gene, levels = ord_genes)]
  fit_dt <- Fit_curves_dt[dosage_gene == dg & gene != dg,  ]
  fit_dt[, gene := factor(gene, levels = ord_genes)]
  
  p <- ggplot(plot_dt, aes(x = dosage_gene_log2FC, y=avg_log2FC)) +
    facet_wrap(gene ~ ., ncol = 7) +
    geom_point(color = "grey20", alpha=0.7) +
    geom_line(data = fit_dt, color="#FF7F00", linewidth = 0.75) +
    theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
    labs(x=paste0(dg, " log2(FC)"), y="Trans gene log2(FC)") +
    coord_cartesian(ylim = c(min(plot_dt$avg_log2FC)-0.1, max(plot_dt$avg_log2FC+0.1)))
  
  ggsave(file.path(plots_dir, paste0("07a_07_AllDosageRespCurves_", dg, ".pdf")), p, width = 5, height = 12)
}




#### Save processed data
# RMSE/AIC data
fwrite(file = file.path(processed_data_dir, "D2N_AIC_LmSig4Loess.txt"), RMSE_dt, quote = F, row.names = F)

# Parameters per gene and cis gene perturbed
fwrite(file = file.path(processed_data_dir, "D2N_S4Param_10fCV.txt"), S4Param_dt, quote = F, row.names = F)

# Fitted dosage response curves (sigmoid or loess)
fwrite(file = file.path(processed_data_dir, "D2N_S4LoessModelFits.txt"), Fit_curves_dt, quote = F, row.names = F)



