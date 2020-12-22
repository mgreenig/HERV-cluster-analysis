suppressPackageStartupMessages(library(dynamicTreeCut))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(pheatmap))

if(!exists("cell_line")){
  cell_line <- commandArgs(trailingOnly = TRUE)
}

if(cell_line != 'Calu3' & cell_line != 'A549'){
  stop('Please input either Calu3 or A549 as the argument after the script name')
}

source(paste('DE_testing_', cell_line, '.R', sep = ''))

# set ggplot theme for increased text size
text_theme <- theme(plot.title = element_text(size = 18, hjust = 0.5), text = element_text(size = 14), plot.subtitle = element_text(size = 16, hjust = 0.5))

# combine significant HERVs with significant human genes
all_sig_genes_infection_expr <- rbind(sig_human_genes_infection_expr, sig_retro_genes_infection_expr) %>%
  as.data.frame

# get network adjacency matrix as 1 - correlation matrix squared
all_sig_genes_dist <- factoextra::get_dist(all_sig_genes_infection_expr, method = 'pearson')

print('Running clustering on human genes and retroelements...')

# cluster and use dynamic tree cut
all_sig_genes_dendro <- hclust(all_sig_genes_dist, method = 'average')
all_sig_genes_clusters <- dynamicTreeCut::cutreeDynamicTree(all_sig_genes_dendro, maxTreeHeight = 1.95, 
                                                            deepSplit = TRUE, minModuleSize = ifelse(cell_line == 'A549', 200, 300))

# add cluster to the dataframe
all_sig_genes_infection$cluster <- all_sig_genes_clusters + 1
all_sig_genes_infection_expr$cluster <- all_sig_genes_clusters + 1
all_sig_genes_infection_expr$gene_type <- ifelse(rownames(all_sig_genes_infection_expr) %in% rownames(sig_retro_genes_infection_expr),
                                                 'Retroelement', 'Human')

print(paste('Found', length(unique(all_sig_genes_clusters)), 'clusters'))

# get size of each cluster
cluster_sizes <- all_sig_genes_infection_expr %>% 
  dplyr::group_by(cluster, gene_type) %>% 
  dplyr::summarise(n = n()) 
cluster_sizes$cluster <- factor(cluster_sizes$cluster, levels = 1:length(unique(all_sig_genes_clusters)))

# plot size of each cluster
cluster_size_plot <- ggplot(cluster_sizes, aes(x = cluster, y = n)) + 
  geom_col(width = 0.75, fill = 'deepskyblue4', colour = 'black') + 
  scale_x_discrete(limits = rev(levels(cluster_sizes$cluster))) + coord_flip() + 
  labs(x = 'Cluster', y = 'Number of genes', title = 'Number of genes per cluster') +
  theme_minimal() + text_theme

if(!interactive()){
  ggsave(paste('figures/cluster_sizes_', cell_line, '.png', sep = ''), cluster_size_plot)
}

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
  cluster_gene_expr <- in_cluster[,colnames(in_cluster) != 'cluster'] %>% select_if(is.numeric)
  
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
  png(paste('figures/dendrogram_gene_clusters_', cell_line, '.png', sep = ''), width = 1080, height = 720, units = 'px')
  plot(eigengene_clustering, ylab = '', xlab = '', sub = '', main = 'Cluster eigengenes', cex.lab = 1.5, cex.main = 2, cex = 1.5)
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
  
  # convert from pearson distance to absolute correlation
  distance_matrix <- abs(1 - distance_matrix)
  
  # get genes in the cluster
  cluster_genes <- rownames(expr[expr$cluster == cluster_number,])
  
  # subset distance matrix
  cluster_dists <- distance_matrix[cluster_genes, cluster_genes]
  
  # iterate over rows and calculate the mean value for each gene
  avg_dists <- apply(cluster_dists, 1, mean)
  avg_dists_sorted <- sort(avg_dists, decreasing = TRUE)

  # get top n genes
  top_n_genes <- names(avg_dists_sorted)[1:n]
  
  return(top_n_genes)
  
}

# get hub genes for all clusters
cluster_hub_genes <- lapply(unique(all_sig_genes_infection_expr$cluster), get_hub_genes, 
                            expr = all_sig_genes_infection_expr, distance_matrix = all_sig_genes_dist, n = 1)
names(cluster_hub_genes) <- unique(all_sig_genes_infection_expr$cluster)

cluster_hub_genes_all <- lapply(unique(all_sig_genes_infection_expr$cluster), function(clus, expr, distance_matrix, n){
  get_hub_genes(cluster_number = clus, expr = expr, distance_matrix = distance_matrix, n = nrow(expr[expr$cluster == clus,]))
}, expr = all_sig_genes_infection_expr, distance_matrix = all_sig_genes_dist)
names(cluster_hub_genes_all) <- unique(all_sig_genes_infection_expr$cluster)

# get retroelement hub genes
retro_hub_genes <- lapply(cluster_hub_genes_all, function(genes, n) genes[genes %in% rownames(sig_retro_genes_infection)][1:n], n = 1)
names(retro_hub_genes) <- unique(all_sig_genes_infection_expr$cluster)

# plot module membership statistics for retro/human hub genes in 4 clusters
if(cell_line == 'Calu3'){
  clusters_to_plot <- c('5', '9', '11', '14')
} else {
  clusters_to_plot <- c('3', '5', '7', '8')
}

# get module membership statistics for each cluster to plot
module_membership_df <- lapply(clusters_to_plot, function(clus){
  eigengene <- cluster_eigengenes[clus,] %>% unlist
  human_hub_gene_expr <- all_sig_genes_infection_expr[cluster_hub_genes[[clus]],] %>% select(-cluster, -gene_type) %>% unlist
  retro_hub_gene_expr <- all_sig_genes_infection_expr[retro_hub_genes[[clus]],]  %>% select(-cluster, -gene_type) %>% unlist
  df <- data.frame(cluster = clus,
                   corr = c(abs(cor(eigengene, human_hub_gene_expr)),
                            abs(cor(eigengene, retro_hub_gene_expr))),
                   type = c('Human', 'Retro'),
                   gene = c(cluster_hub_genes[[clus]],
                            retro_hub_genes[[clus]]))
  return(df)
  }) %>%
  do.call(rbind, .)
module_membership_df$gene <- gsub('\\..*$', '', module_membership_df$gene)

# plot module memberships
module_membership_plot <- ggplot(module_membership_df, aes(x = cluster, y = corr, fill = type, label = gene)) + 
  geom_col(position = position_dodge(0.6), width = 0.5, colour = 'black') +
  labs(x = '\nCluster', y = 'Module membership statistic\n', fill = 'Hub gene',
       title = 'Module memberships of human and retroelement hub genes\n') +
  geom_label(position = position_dodge(0.75), vjust = -0.25, size = 4, show.legend = FALSE) + 
  theme_minimal() + theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5),
                          text = element_text(size = 18)) + ylim(c(0,1.05)) +
  scale_fill_manual(values = c('#00BFC4', '#F8766D'), labels = c('Human', 'Retroelement'))

if(!interactive()){
  ggsave(paste('figures/module_membership_plot_', cell_line, '.png', sep = ''),
         module_membership_plot, width = 12, height = 6, units = 'in')
}
