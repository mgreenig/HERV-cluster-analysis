# Arguments for this script are:
# --clusterdata : filepath to the cluster data object ClusterData.Rdata outputted by clustering.R
# --metadata : filepath to metadata file (comma-separated), with three columns: "sample", "condition", and "batch" (batch can be omitted)
# --outdir : directory to save data files (this argument can be omitted, in which case data files will be saved to the current directory)

options(warn = -1)

source('utils.R')

reqs <- c('reshape2', 'ggpubr', 'ggplot2', 'dplyr', 'stringr')
get_reqs(reqs)

suppressMessages(library(reshape2))
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = T)
parsed_args <- parseArgs(args)
names(parsed_args) <- stringr::str_to_lower(names(parsed_args))
req_args <- c('clusterdata', 'metadata')

checkArgs(parsed_args, req_args, file_args = c('clusterdata', 'metadata'))

load(parsed_args['clusterdata'])

# import metadata
metadata <- read.csv(parsed_args['metadata'], sep = ',', header = TRUE)
colnames(metadata) <- stringr::str_to_lower(colnames(metadata))

if(!('condition' %in% colnames(metadata))){
  stop('Condition column not found in metadata file.')
}

if(!('sample' %in% colnames(metadata))){
  stop('Sample column not found in metadata file.')
}
rownames(metadata) <- metadata$sample
metadata$sample <- NULL
metadata$condition <- as.factor(metadata$condition)

if('batch' %in% colnames(metadata)){
  metadata$batch <- as.factor(metadata$batch)
}

# if output directory does not exist, create
if('outdir' %in% names(parsed_args)){
  if(!dir.exists(parsed_args['outdir'])){
    dir.create(parsed_args['outdir'])
  }
}

# set ggplot theme for increased text size
text_theme <- theme(plot.title = element_text(size = 18, hjust = 0.5), text = element_text(size = 14), plot.subtitle = element_text(size = 16, hjust = 0.5))

# get size of each cluster
cluster_sizes <- clusterdata[['expr']] %>% 
  dplyr::group_by(cluster, gene_type) %>% 
  dplyr::summarise(n = n()) 
cluster_sizes$cluster <- factor(cluster_sizes$cluster, levels = 1:length(unique(clusterdata[['expr']][['cluster']])))

# plot size of each cluster
cluster_size_plot <- ggplot(cluster_sizes, aes(x = cluster, y = n)) + 
  geom_col(width = 0.75, fill = 'deepskyblue4', colour = 'black') + 
  scale_x_discrete(limits = rev(levels(cluster_sizes$cluster))) + coord_flip() + 
  labs(x = 'Cluster', y = 'Number of genes', title = 'Number of genes per cluster') +
  theme_minimal() + text_theme

if(!dir.exists('figures')){
  dir.create('figures')
}

writeLines('Saving cluster size plot to figures/cluster_sizes.png...')

# save cluster sizes figure
ggsave('figures/cluster_sizes.png', cluster_size_plot)

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
cluster_eigengenes <- lapply(unique(clusterdata[['expr']]$cluster), get_eigengene, expr = clusterdata[['expr']]) %>%
  do.call(rbind, .) %>%
  as.data.frame
rownames(cluster_eigengenes) <- unique(clusterdata[['expr']]$cluster) 

# get distance matrix for the eigengenes and cluster
eigengene_dists <- factoextra::get_dist(cluster_eigengenes, method = 'pearson')
eigengene_clustering <- hclust(eigengene_dists, method = 'average')

writeLines('Saving cluster dendrogram to figures/dendrogram_gene_clusters.png...')

png('figures/dendrogram_gene_clusters.png', width = 1080, height = 720, units = 'px')
plot(eigengene_clustering, ylab = '', xlab = '', sub = '', main = 'Cluster eigengenes', cex.lab = 1.5, cex.main = 2, cex = 1.5)
dev.off()

# get module membership statistics for each cluster to plot
module_membership_df <- lapply(unique(clusterdata[['expr']]$cluster), function(clus){
  eigengene <- cluster_eigengenes[clus,] %>% unlist
  human_hub_gene_expr <- clusterdata[['expr']][clusterdata[['hub_genes']][[clus]][1],] %>% select(-cluster, -gene_type) %>% unlist
  retro_hub_gene_expr <- clusterdata[['expr']][clusterdata[['TE_hub_genes']][[clus]],]  %>% select(-cluster, -gene_type) %>% unlist
  df <- data.frame(cluster = clus,
                   corr = c(abs(cor(eigengene, human_hub_gene_expr)),
                            abs(cor(eigengene, retro_hub_gene_expr))),
                   type = c('Human', 'Retro'),
                   gene = c(clusterdata[['hub_genes']][[clus]][1],
                            clusterdata[['TE_hub_genes']][[clus]]))
  return(df)
}) %>% do.call(rbind, .)
module_membership_df$gene <- gsub('\\..*$', '', module_membership_df$gene)
module_membership_df$cluster <- factor(module_membership_df$cluster, levels = unique(clusterdata[['expr']]$cluster))

# plot module memberships
module_membership_plot <- ggplot(module_membership_df, aes(x = cluster, y = corr, fill = type)) + 
  geom_col(position = position_dodge(0.6), width = 0.5, colour = 'black') +
  labs(x = '\nCluster', y = 'Module membership statistic\n', fill = 'Hub gene',
       title = 'Module memberships of human and retroelement hub genes\n') + 
  theme_minimal() + theme(plot.title = element_text(size = 20, face = 'bold', hjust = 0.5),
                          text = element_text(size = 18)) + ylim(c(0,1.05)) +
  scale_fill_manual(values = c('#00BFC4', '#F8766D'), labels = c('Human', 'Retroelement'))

ggsave('figures/module_membership_plot.png', width = 12, height = 6, units = 'in', module_membership_plot)


# function for getting the mean expression for each gene per condition for a gene expression data frame
get_mean_expr_per_group <- function(expr, coldata, condition_col = 'condition'){
  
  # Function takes in an expression data frame and a metadata dataframe 
  # and returns a dataframe containing the mean transcript count for each gene per condition
  # 
  # Arguments
  # ------------
  # 
  # expr : n x p dataframe with n rows (genes) and p columns (samples)
  # coldata : metadata frame with p rows corresponding to columns in expr, 
  #           and a condition column denoting the category of each row
  # condition_col : column of coldata that corresponds to the condition of interest, 
  #                 for which mean expression per group will be calculated
  # 
  # Returns
  # ------------
  # 
  # mean_expr_per_group : n x k dataframe with n rows (genes) and k columns 
  #                       (number of unique levels of condition_col)
  
  # get conditions from metadata
  
  conds <- unique(coldata[[condition_col]])
  
  # get mean absolute value of expression per condition group, put into data frame
  mean_expr_per_group <- lapply(conds, function(cond){
    
    cols <- rownames(coldata)[coldata[[condition_col]] == cond]
    in_group_expr <- expr[,cols] 
    
    # calculate mean expression per row (group)
    mean_expr_in_group <- apply(in_group_expr, 1, mean)

    return(mean_expr_in_group)
    
  }) %>%
    do.call(cbind, .) %>%
    as.data.frame
  
  # set column names to conditions
  colnames(mean_expr_per_group) <- conds
  
  return(mean_expr_per_group)
  
}

# function for plotting the expression of genes in a cluster
plot_trajectory <- function(mean_expr_df, cluster_number, retro_gene_list, 
                            cluster_hub_genes, retro_hub_genes,
                            conditon_var = 'condition'){
  
  # Function takes in an expression data frame, a cluster number, and a list of all HERV genes
  # and returns dot plot 
  # 
  # Arguments
  # ------------
  # 
  # mean_expr_df : n x k dataframe with n rows (genes) and k columns (condition levels), 
  #                output of get_mean_expr_per_group()
  # cluster_number : number of cluster being plotted, for plot title
  # retro_gene_list : list of all relevant retroelements, so retroelements in mean_expr_df can be identified
  # 
  # Returns
  # ------------
  # 
  # boxplot : box plot of differential expression for each condition level
  
  # add column for human/HERV genes
  mean_expr_df$gene_type <- ifelse(rownames(mean_expr_df) %in% retro_gene_list, 'Retroelement', 'Human') %>% 
    factor(levels = c('Human', 'Retroelement'))
  
  # compress data frame into two columns, one for infection group and the other for expression
  plot_df <- reshape2::melt(mean_expr_df, id.vars = 'gene_type', variable.name = conditon_var, value.name = 'expression')
  plot_df$HERV_expr[plot_df$gene_type == 'Retroelement'] <- plot_df$expression[plot_df$gene_type == 'Retroelement']
  
  # make plot
  boxplot <- ggplot(plot_df, aes(x = !!sym(conditon_var), y = expression)) + 
    geom_jitter(aes(colour = gene_type), width = 0.25, alpha = 0.75)  +
    geom_boxplot(aes(fill = gene_type), alpha = 0.75, outlier.shape = NA, coef = 0, show.legend = FALSE) +
    labs(y = 'Scaled expression\n', x = paste('\n', stringr::str_to_sentence(conditon_var), sep = ''), 
         title = paste('Cluster', cluster_number, '|', nrow(mean_expr_df), 'genes'),
         subtitle = paste('Hub gene:', cluster_hub_genes[[as.character(cluster_number)]], 
                          '\nHub retroelement:', retro_hub_genes[[as.character(cluster_number)]]),
         colour = 'Gene type') + 
    theme_minimal() + text_theme + theme(plot.title = element_text(face = 'bold')) + 
    scale_color_manual(values = c('#00BFC4', '#F8766D')) + 
    scale_fill_manual(values = c('#00BFC4', '#F8766D'))
  
  return(boxplot)
  
}

writeLines('\nPlotting cluster trajectories...')

# plot trajectory for every cluster
plots <- lapply(unique(clusterdata[['expr']]$cluster), function(cluster){
  
  # get expression of genes in the cluster
  in_cluster <- clusterdata[['expr']][clusterdata[['expr']]$cluster == cluster,]
  in_cluster_expr <- in_cluster[,colnames(in_cluster) != 'cluster']
  
  # get mean expression per condition level in the cluster
  mean_expr_per_group <- get_mean_expr_per_group(expr = in_cluster_expr, coldata = metadata)

  # plot the trajectory  
  plot_trajectory(cluster_number = cluster, mean_expr_df = mean_expr_per_group, 
                  retro_gene_list = rownames(subset(clusterdata[['expr']], gene_type == 'Retroelement')),
                  cluster_hub_genes = clusterdata[['hub_genes']],
                  retro_hub_genes = clusterdata[['TE_hub_genes']])
  
})
# set names to cluster numbers
names(plots) <- unique(clusterdata[['expr']]$cluster)

fig_dir <- 'clusters'
if(!dir.exists(paste0('figures/', fig_dir))){
  dir.create(paste0('figures/', fig_dir), showWarnings = FALSE)
}
for(plot in names(plots)){
  filepath <- paste0('figures/', fig_dir, '/', plot, '_trajectory.png')
  ggsave(filepath, plots[[plot]])
}