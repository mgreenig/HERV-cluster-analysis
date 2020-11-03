library(dynamicTreeCut)
library(factoextra)
library(pheatmap)

source('DE_testing.R')

# set ggplot theme for increased text size
text_theme <- theme(plot.title = element_text(size = 18, hjust = 0.5), text = element_text(size = 14))

# combine significant HERVs with significant human genes
all_sig_genes_infection_expr <- rbind(sig_human_genes_infection_expr, sig_retro_genes_infection_expr) %>%
  as.data.frame

# get network adjacency matrix as 1 - correlation matrix squared
all_sig_genes_dist <- factoextra::get_dist(all_sig_genes_infection_expr, method = 'pearson')

# cluster and use dynamic tree cut
all_sig_genes_dendro <- hclust(all_sig_genes_dist, method = 'average')
all_sig_genes_clusters <- dynamicTreeCut::cutreeDynamicTree(all_sig_genes_dendro, maxTreeHeight = 1.95, deepSplit = TRUE, minModuleSize = 200)

# add cluster to the dataframe
all_sig_genes_infection$cluster <- all_sig_genes_clusters + 1
all_sig_genes_infection_expr$cluster <- all_sig_genes_clusters + 1

# get size of each cluster
cluster_sizes <- all_sig_genes_infection_expr %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::summarise(n = n()) 
cluster_sizes$cluster <- factor(cluster_sizes$cluster, levels = 1:length(unique(all_sig_genes_clusters)))

# plot size of each cluster
cluster_size_plot <- ggplot(cluster_sizes, aes(x = cluster, y = n)) + 
  geom_col(width = 0.75, fill = 'deepskyblue4', colour = 'black') + 
  scale_x_discrete(limits = rev(levels(cluster_sizes$cluster))) + coord_flip() + 
  labs(x = 'Cluster', y = 'Number of genes', title = 'Number of genes per cluster') +
  theme_minimal() + text_theme

# function for getting eigengene for a cluster
get_eigengene <- function(expr, cluster_number){
  # Function that takes in a count matrix and a cluster number and returns the cluster eigengene
  # 
  # Arguments
  # ------------
  # 
  # counts : transcript count matrix with samples as rows and genes as columns, including a column named cluster
  # cluster_number : which cluster to get the eigengene for 
  # 
  # 
  # Returns
  # ------------
  # 
  # eigengene : first principal component of the cluster expression matrix
  
  
  # get genes in the cluster
  in_cluster <- expr[expr$cluster == cluster_number,]
  cluster_gene_expr <- in_cluster[,colnames(in_cluster) != 'cluster']
  
  # run PCA on cluster genes
  PCA_res <- prcomp(cluster_gene_expr, scale. = FALSE, center = FALSE)
  
  # get eigengene as the 1st principal component
  eigengene <- PCA_res$rotation[,c('PC1')]
  
  return(eigengene)
  
}

# get eigengenes for all clusters
cluster_eigengenes <- lapply(unique(all_sig_genes_infection_expr$cluster), get_eigengene, expr = all_sig_genes_infection_expr) %>%
  do.call(rbind, .) %>%
  as.data.frame
rownames(cluster_eigengenes) <- unique(all_sig_genes_infection_expr$cluster) 

# get distance matrix for the eigengenes and cluster
eigengene_dists <- factoextra::get_dist(cluster_eigengenes, method = 'pearson')
eigengene_clustering <- hclust(eigengene_dists, method = 'average')

if(!interactive()){
  png('figures/dendrogram_gene_clusters.png', width = 1280, height = 720, units = 'px')
  plot(eigengene_clustering)
  dev.off()
}

# function for getting the hub gene(s) for a cluster
get_hub_genes <- function(expr, distance_matrix, cluster_number, n = 1){
  # Function that takes in a count matrix and a cluster number and returns the cluster hub genes,
  # defined as the genes with the highest average connectivity (lowest distance) to other genes in the cluster
  # 
  # Arguments
  # ------------
  # 
  # counts : transcript count matrix with samples as rows and genes as columns, including a column named cluster
  # distance_matrix : distance matrix used for the clustering
  # cluster_number : which cluster to get the eigengene for 
  # n : number of hub genes to return
  # 
  # 
  # Returns
  # ------------
  # 
  # top_n_genes : the top n hub genes for the cluster
  
  # convert to matrix if input is a dist object
  if(is(distance_matrix, 'dist')){
    distance_matrix <- as.matrix(distance_matrix)
  }
  
  # get genes in the cluster
  cluster_genes <- rownames(expr[expr$cluster == cluster_number,])
  
  # subset distance matrix
  cluster_dists <- distance_matrix[cluster_genes, cluster_genes]
  
  # iterate over rows and calculate the mean value for each gene
  avg_dists <- apply(cluster_dists, 1, mean)
  avg_dists_sorted <- sort(avg_dists)

  # get top n genes
  top_n_genes <- names(avg_dists_sorted)[1:n]
  
  return(top_n_genes)
  
}

# get hub genes for all clusters
cluster_hub_genes <- lapply(unique(all_sig_genes_infection_expr$cluster), get_hub_genes, 
                            expr = all_sig_genes_infection_expr, distance_matrix = all_sig_genes_dist, n = 1)
names(cluster_hub_genes) <- unique(all_sig_genes_infection_expr$cluster)

# get expression vectors for cluster hub genes
cluster_hub_gene_expr <- lapply(names(cluster_hub_genes), function(cluster){
  
  # get counts for genes in the cluster
  in_cluster <- all_sig_genes_infection_expr[all_sig_genes_infection_expr$cluster == cluster,]
  cluster_counts <- scale(in_cluster[,colnames(in_cluster) != 'cluster'])
  
  # get the count vector for the hub gene
  gene_name <- cluster_hub_genes[[cluster]]
  gene_vector <- cluster_counts[gene_name,]
    
  return(gene_vector)
  
}) %>% 
  do.call(rbind, .) %>%
  as.data.frame
rownames(cluster_hub_gene_expr) <- cluster_hub_genes

# cluster hub genes
hub_gene_dists <- factoextra::get_dist(cluster_hub_gene_expr, method = 'pearson')
hub_gene_clustering <- hclust(hub_gene_dists, method = 'average')

# order eigengene clusters 1 to k
ordered_rows <- rownames(cluster_eigengenes) %>% 
  as.numeric %>% 
  order

# get adjacency between hub and eigengenes (correlation squared)
hub_and_eigengene_dists <- cor(t(cluster_eigengenes[ordered_rows,]), 
                               t(cluster_hub_gene_expr[ordered_rows,]), 
                               method = 'pearson') ^ 2
# order dataframe
hub_and_eigengene_dists_ordered <- t(hub_and_eigengene_dists)

# make heatmap of hub gene/eigengene correlations
hub_and_eigengene_heatmap <- pheatmap::pheatmap(hub_and_eigengene_dists_ordered, cluster_rows = FALSE, 
                                                cluster_cols = FALSE, angle_col = 0, 
                                                main = 'Eigengene and hub gene correlations across modules',
                                                filename = 'figures/eigengene_hubgene_comparison.png', 
                                                width = 8, height = 5)
