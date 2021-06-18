# Arguments for this script are:
# --counts : filepath to the counts object counts.Rdata outputted by differentialexpression.R
# --diffexp : filepath to the differential expression object DiffExp.Rdata outputted by differentialexpression.R
# --geneset : one of: ('LRT', 'All') defining the gene set to use in clustering. LRT only clusters genes/TEs with adj p < 0.05 from the LRT. Defaults to LRT.
# --minclustersize : integer specifying the minimum cluster size (default is 200, can be omitted)
# --outdir : directory to save data files (this argument can be omitted, in which case data files will be saved to the current directory)

options(warn = -1)

source('utils.R')

reqs <- c('dynamicTreeCut', 'factoextra', 'pheatmap', 'dplyr', 'stringr')
get_reqs(reqs)

suppressMessages(library(dynamicTreeCut, quietly = T))
suppressMessages(library(factoextra, quietly = T))
suppressMessages(library(pheatmap, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(stringr, quietly = T))

args <- commandArgs(trailingOnly = TRUE)
parsed_args <- parseArgs(args)
names(parsed_args) <- stringr::str_to_lower(names(parsed_args))
req_args <- c('counts', 'diffexp')

checkArgs(parsed_args, req_args, file_args = c('counts', 'diffexp'))

# check if min cluster size specified
if('minclustersize' %in% names(parsed_args)){
  min_cluster_size <- as.integer(parsed_args['minclustersize'])
} else {
  min_cluster_size <- 200
}

# check if genes to cluster is specified
if(!('geneset' %in% names(parsed_args))){
  parsed_args['geneset'] <- 'lrt'
}

# check if output directory specified, create if needed
if('outdir' %in% names(parsed_args)){
  if(!dir.exists(parsed_args['outdir'])){
    dir.create(parsed_args['outdir'])
  }
}

parsed_args['geneset'] <- str_to_lower(parsed_args['geneset'])
# get counts of significant genes identified by the LRT
if(parsed_args['geneset'] != 'lrt' & 
   parsed_args['geneset'] != 'all'){
  stop('Argument --geneset should be either "LRT" or "All"')
}

load(parsed_args['counts'])
load(parsed_args['diffexp'])

if(parsed_args['geneset'] == 'lrt'){
  sig_human_genes <- subset(DE_results[['human']][['LRT']], padj < 0.05)
  sig_retro_genes <- subset(DE_results[['TE']][['LRT']], padj < 0.05)
  human_genes_expr <- counts[['human_reg']][rownames(sig_human_genes),]
  retro_genes_expr <- counts[['TE_reg']][rownames(sig_retro_genes),]
  all_genes_DE <- rbind(sig_human_genes, sig_retro_genes)
} else {
  human_genes_expr <- counts[['human_reg']]
  retro_genes_expr <- counts[['TE_reg']]
  all_genes_DE <- rbind(DE_results[['human']][['LRT']], DE_results[['TE']][['LRT']])
}

# combine significant HERVs with significant human genes
all_genes_expr <- rbind(human_genes_expr, retro_genes_expr) %>%
  as.data.frame

writeLines(paste('\nRunning clustering on', nrow(human_genes_expr), 'human genes and', 
                 nrow(retro_genes_expr), 'transposable elements...\n'))

# get network adjacency matrix as 1 - correlation matrix squared
all_sig_genes_dist <- factoextra::get_dist(all_genes_expr, method = 'pearson')

# cluster and use dynamic tree cut
all_sig_genes_dendro <- hclust(all_sig_genes_dist, method = 'average')
all_sig_genes_clusters <- dynamicTreeCut::cutreeDynamicTree(all_sig_genes_dendro, maxTreeHeight = 1.95, 
                                                            deepSplit = TRUE, minModuleSize = min_cluster_size)

# add cluster to the dataframe
all_genes_expr$cluster <- all_sig_genes_clusters + 1
all_genes_expr$gene_type <- ifelse(rownames(all_genes_expr) %in% rownames(retro_genes_expr),
                                                 'Retroelement', 'Human')

writeLines(paste('Found', length(unique(all_sig_genes_clusters)), 'clusters\n'))

# write genes and their clusters to CSV
if('outdir' %in% names(parsed_args)){
  cluster_data_filepath <- paste0(parsed_args['outdir'], '/clustered_genes.csv')
} else {
  cluster_data_filepath <- 'clustered_genes.csv'
}

# save cluster data file
writeLines(paste('Saving cluster data to:', cluster_data_filepath, '\n'))
all_genes_DE$cluster <- all_genes_expr[rownames(all_genes_DE), 'cluster']
write.csv(all_genes_DE[,c('pvalue', 'padj', 'cluster')], cluster_data_filepath, row.names = T)

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
cluster_hub_genes <- lapply(unique(all_genes_expr$cluster), function(clus, expr, distance_matrix){
  get_hub_genes(cluster_number = clus, expr = expr, distance_matrix = distance_matrix, n = nrow(expr[expr$cluster == clus,]))
}, expr = all_genes_expr, distance_matrix = all_sig_genes_dist)
names(cluster_hub_genes) <- unique(all_genes_expr$cluster)

# get retroelement hub genes
retro_hub_genes <- lapply(cluster_hub_genes, function(genes, n) genes[genes %in% rownames(retro_genes_expr)][1:n], n = 1)
names(retro_hub_genes) <- unique(all_genes_expr$cluster)

# compile results into a single list
clusterdata <- list('expr' = all_genes_expr,
                    'hub_genes' = cluster_hub_genes, 
                    'TE_hub_genes' = retro_hub_genes)

if('outdir' %in% names(parsed_args)){
  clusterdata_filepath <- paste0(parsed_args['outdir'], '/ClusterData.Rdata')
} else {
  clusterdata_filepath <- '/ClusterData.Rdata'
}

writeLines(paste('Saving cluster data object to:', clusterdata_filepath, '\n'))

save(clusterdata, file = clusterdata_filepath)

