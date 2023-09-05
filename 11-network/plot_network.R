library(ggplot2)
library(tidygraph)
library(ggnetwork)
library(igraph)
library(data.table)
library(dplyr)
library(ggraph)

data = fread("../../../chip_seq_ann/data/remap2022/data/TF_target_mapping_filtered_merged_K562_with_motifs_with_ppi_with_coexpr_with_dnase_with_atac_with_dist_score.tsv", nThread=10)
data %>% filter(is_atac) -> data
d2n_genes = fread("../../../chip_seq_ann/data/d2n/d2n_genes.txt")
data %>% filter((tf %in% d2n_genes$external_gene_name) & (gene_symbol %in% d2n_genes$external_gene_name)) %>% select(tf, gene_symbol) %>% distinct() -> d2n_net

net <- igraph::graph_from_data_frame(d2n_net, directed=T) 
net <- igraph::simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
deg <- igraph::degree(net, mode = "all", normalized = TRUE) 
bet = igraph::betweenness(net, mode = "all", normalized = TRUE)
tg <- tidygraph::as_tbl_graph(net) %>% 
  tidygraph::activate(nodes) %>% 
  dplyr::mutate(label=name,
		size = deg)


plot <- tg %>%
   activate(edges) %>%
   mutate(bc = centrality_edge_betweenness()) %>%
   ggraph(layout = "centrality", centrality=deg) +
   geom_edge_arc(aes(colour= bc, width= bc),
                  lineend = "round",
                 strength = .1,
		arrow = arrow(length = unit(1, 'mm')),
		end_cap = circle(1, 'mm')) +
   geom_node_point(aes(color=bet), size= log10(deg * 10)) +
   geom_node_text(aes(label = name), 
                  repel = TRUE, 
                  point.padding = unit(0.5, "lines"), 
                  size= 2) +
  scale_edge_width(range = c(0.02, 0.5)) +
  scale_color_viridis(option="rocket", direction=-1) +
  scale_edge_colour_viridis(option = "rocket", direction=-1, guide = guide_colorbar(available_aes = "edge_colour")) +
  theme_graph(background = "white") +
  guides(edge_width = FALSE,
         edge_alpha = FALSE,
	color=F)

ggsave(plot = plot, filename = "../../plots/10-network/try.png", dpi=720, width=210, height=210, units="mm",device = "png")

