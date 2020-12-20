library(dplyr)

# import from script running differential expression analysis
cell_line <- 'Calu3'

source(paste('DE_testing_', cell_line, '.R', sep = ''))

# import HERV annotation table
retro_annotations <- read.table('data/Retro_annotation.tsv', col.names = c('Locus','Class','Family','Category','Chrom','Start','End','Strand'))

# function for going from producing a p-value table from annotations of significant genes
get_feature_enrichment <- function(annotation_df, sig_genes, feature, adjust.p = TRUE, adjust.p.method = 'BH'){
  
  # Function takes in an annotation dataframe, a list of DE genes, and a feature of interest
  # and returns a dataframe containing information on enrichment of the feature in the DE genes
  # 
  # Arguments
  # ------------
  # 
  # annotation_df : n x p dataframe with n rows (genes) and p columns (features/annotations)
  # sig_genes : list of differentially-expressed genes
  # feature : column of annotation_df that corresponds to the feature of interest, 
  #           for which enrichment within sig_genes will be rested
  # 
  # Returns
  # ------------
  # 
  # sig_genes_by_feature : k x 4 dataframe with k rows corresponding to levels of feature
  #                        and adjusted p-value for each level value
  
  # get annotations for statistically-significant loci
  sig_gene_annotations <- annotation_df[annotation_df$Locus %in% sig_genes,]
  
  # get number of significant elements per feature and the total number of elements in each family
  sig_genes_by_feature <- dplyr::group_by(sig_gene_annotations, !!sym(feature)) %>% 
    summarise(n_sig = n()) 
  
  # get total number of elements in each feature category
  sig_genes_by_feature$n_category <- sapply(sig_genes_by_feature[[feature]], function(category){
    n <- annotation_df$Locus[annotation_df[,feature] == category] %>%
      unique %>%
      length
    return(n)
  })

  # get an enrichment p-value for each family
  sig_genes_by_feature$p_value <- apply(sig_genes_by_feature, 1, function(row, all_sig, n_total){
    
    # get number of significant HERVs in each category and 
    overlap <- as.integer(row['n_sig'])
    family_total <- as.numeric(row['n_category'])
    
    # hypergeometric test for HERV enrichment
    p_value <- phyper(overlap - 1, family_total, n_total - family_total, all_sig, lower.tail = FALSE)
    return(p_value)
    
  }, all_sig = sum(sig_genes_by_feature$n_sig), n_total = length(unique(annotation_df$Locus)))
  
  # adjust p-values if necessary
  if(adjust.p == TRUE){
    sig_genes_by_feature$padj <- p.adjust(sig_genes_by_feature$p_value, method = adjust.p.method)
  }
  
  # sort by p-value
  sig_genes_by_feature <- dplyr::arrange(sig_genes_by_feature, p_value)
  
  return(sig_genes_by_feature)
  
}

# get enrichment for different families
retro_family_enrichment_Cov2 <- get_feature_enrichment(retro_annotations, sig_genes = rownames(sig_retro_genes_Cov2),
                                                       feature = 'Family')
retro_family_enrichment_infection <- get_feature_enrichment(retro_annotations, sig_genes = rownames(sig_retro_genes_infection),
                                                            feature = 'Family')

# get enrichment for different categories
retro_category_enrichment_Cov2 <- get_feature_enrichment(retro_annotations, sig_genes = rownames(sig_retro_genes_Cov2), 
                                                         feature = 'Category')
retro_category_enrichment_infection <-  get_feature_enrichment(retro_annotations, sig_genes = rownames(sig_retro_genes_infection),
                                                               feature = 'Category')