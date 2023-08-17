library(foreach)
library(ggplot2); theme_set(theme_bw())
library(data.table)
library(dplyr)



## Dir
setwd("./")
folder_name <- basename(getwd())
data_dir = file.path("../../data", folder_name)
processed_data_dir = file.path("../../processed_data", folder_name)
plots_dir = file.path("../../plots", folder_name)


sigm_params = fread(paste0(processed_data_dir, "/D2N_S4Param_10fCV.txt"))
net_params = fread(paste0(processed_data_dir, "/topological_metrics_TRN_remap2022.tsv"))

net_params %>% select(-c("n_components")) -> net_params
sigm_params %>% mutate(abs_slope_IF = abs(slope_IF)) -> sigm_params
colnames(net_params)[1] = "gene"

net_params_melt = melt(net_params, id = c("gene"))
colnames(net_params_melt) = c("gene", "metric", "value")

## (4) Correlation between centrality metrics and all sigmoid features

sig_feat <-  c("slope_IF", "abs_slope_IF","min_assmp","max_assmp", "min_max_range", "x_IF")
Cor_dt <- foreach(dg = c("GFI1B", "MYB", "NFE2"), .combine = rbind) %do% {
  dt1 <- merge.data.table(net_params_melt, sigm_params[dosage_gene == dg & unresponsive == F & non_monotonic == F, ], by = "gene")
  dt2 <- dt1[!is.na(value), {
    lapply(.SD[, sig_feat, with = FALSE], function(x) {
      cor.test(as.numeric(value), x, method = "spearman")$p.value
    })
  }, by = .(metric)]
  dt3 <- dt1[!is.na(value), {
    lapply(.SD[, sig_feat, with = FALSE], function(x) {
      cor.test(as.numeric(value), x, method = "spearman")$estimate
    })
  }, by = .(metric)]
  
  dt4 <- merge.data.table(melt(dt2, id.vars = c("metric"), variable.name = "sigmoid_param", value.name = "pval"),
                          melt(dt3, id.vars = c("metric"), variable.name = "sigmoid_param", value.name = "r"),
                          by = c("metric", "sigmoid_param"))
  dt4[, dosage_gene := dg]
}
Cor_dt[, pval := p.adjust(pval, "fdr")]
Cor_dt[, log10_pval := -log10(pval)]


p <- ggplot(Cor_dt, aes(x= sigmoid_param, y = metric)) +
  geom_point(aes(color=r, size=log10_pval, alpha=log10_pval > 1.3)) +
  facet_grid(~ dosage_gene, scales = "free_y", space = "free_y") +
  scale_color_gradient2(low = "#377EB8", high = "#FF7F00", midpoint = 0) +
  guides(alpha = "none") + 
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = p, filename = paste0(plots_dir, "/08c_NetworkTopologicalProperties.png"), device = "png", dpi = 600, width=8, height=5)
