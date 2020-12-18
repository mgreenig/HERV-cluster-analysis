library(reshape2)
library(ggpubr)

source('full_clustering.R')

# function for getting the mean expression for each gene per condition for a gene expression data frame
get_mean_expr_per_group <- function(expr, coldata, condition_col = 'Infection'){
  
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
    
    cols <- rownames(coldata)[coldata$Infection == cond]
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
plot_trajectory <- function(mean_expr_df, cluster_number, retro_gene_list, conditon_var = 'Infection'){
  
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
    labs(y = 'Scaled expression\n', x = paste('\n', conditon_var, sep = ''), 
         title = paste('Cluster', cluster_number, '|', nrow(mean_expr_df), 'genes'),
         subtitle = paste('Hub gene (human):', cluster_hub_genes[[as.character(cluster_number)]], 
                          '\nHub gene (retro):', retro_hub_genes[[as.character(cluster_number)]]),
         colour = 'Gene type') + 
    theme_minimal() + text_theme + theme(plot.title = element_text(face = 'bold')) + 
    scale_color_manual(values = c('#00BFC4', '#F8766D')) + 
    scale_fill_manual(values = c('#00BFC4', '#F8766D'))
  
  return(boxplot)
   
}
 
# plot trajectory for every cluster
plots <- lapply(unique(all_sig_genes_infection_expr$cluster), function(cluster){
  
  # get expression of genes in the cluster
  in_cluster <- all_sig_genes_infection_expr[all_sig_genes_infection_expr$cluster == cluster,]
  in_cluster_expr <- in_cluster[,colnames(in_cluster) != 'cluster']

  # get mean expression per condition level in the cluster
  mean_expr_per_group <- get_mean_expr_per_group(expr = in_cluster_expr, coldata = metadata)
  
  # plot the trajectory  
  plot_trajectory(cluster_number = cluster, mean_expr_df = mean_expr_per_group, 
                  retro_gene_list = rownames(sig_retro_genes_infection_expr))
  
})
# set names to cluster numbers
names(plots) <- unique(all_sig_genes_infection_expr$cluster)

if(!interactive()){
  fig_dir <- paste(cell_line, '_clusters/', sep = '')
  dir.create(fig_dir, showWarnings = FALSE)
  for(plot in names(plots)){
    filepath <- paste('figures/', fig_dir, plot, '_trajectory.png', sep = '')
    ggsave(filepath, plots[[plot]])
  }
  
  legend_theme <- theme(legend.text = element_text(size = 14)) 
  legend_colors <- guides(colour = guide_legend(override.aes = list(size=6)))
  
  top_4_clusters <- ggarrange(plots[[clusters_to_plot[1]]] + legend_colors + legend_theme, NULL,
                              plots[[clusters_to_plot[2]]] + legend_colors + legend_theme, NULL,
                              plots[[clusters_to_plot[3]]] + legend_colors + legend_theme, NULL,
                              plots[[clusters_to_plot[4]]] + legend_colors + legend_theme, NULL,
                              common.legend = T, nrow = 1,
                              widths = c(1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1))
  
  ggsave(paste('figures/top_clusters_', cell_line, '.png', sep = ''), top_4_clusters, width = 18, height = 6, units = 'in')
}
