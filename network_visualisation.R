library(igraph)

source('clustering.R')

# function to get adjacency matrix from transcript counts
get_adjacency_matrix <- function(counts, method = c('pearson', 'spearman', 'kendall'), threshold = 0.85){
   
  corr_mat <- 1 - get_dist(counts, method = method) %>% as.matrix
  
  adj_mat <- ifelse(corr_mat > threshold, 1, 0)
  
  return(adj_mat)
  
}

retro_adj_mat <- get_adjacency_matrix(sig_retro_genes_infection_counts, method = 'pearson')
retro_graph <- graph_from_adjacency_matrix(retro_adj_mat)
layout <- layout_with_fr(retro_graph)

simplified_retro_graph <- igraph::simplify(retro_graph, remove.multiple=T, remove.loops = T, edge.attr.comb=list(Weight="sum","ignore"))
plot(simplified_retro_graph, edge.arrow.size = 0.4, layout = layout)

gene_sample <- dplyr::sample_n(as.data.frame(sig_human_genes_infection_counts), 100)
gene_adj_mat <- get_adjacency_matrix(gene_sample, method = 'pearson')
gene_graph <- graph_from_adjacency_matrix(gene_adj_mat)
layout <- layout_with_fr(gene_graph)

simplified_gene_graph <- igraph::simplify(gene_graph, remove.multiple=T, remove.loops = T, edge.attr.comb=list(Weight="sum","ignore"))
plot(simplified_gene_graph, edge.arrow.size = 0.4, layout = layout)
