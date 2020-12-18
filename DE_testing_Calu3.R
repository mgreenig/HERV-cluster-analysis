library(VennDiagram)
library(EnhancedVolcano)
library(gridExtra)

source('preprocessing_Calu3.R')

# LRT to get genes for clustering
dds_human_LRT <- DESeq(dds_human, test = 'LRT', reduced = ~1)
dds_human_LRT_infection_results <- results(dds_human_LRT)
sig_human_genes_infection <- subset(dds_human_LRT_infection_results, padj < 0.05)

# get all sig genes for human LRT
all_sig_human_gene_names_LRT <- rownames(sig_human_genes_infection)

# LRT to get genes for clustering
dds_retro_LRT <- DESeq(dds_retro, test = 'LRT', reduced = ~1)
dds_retro_LRT_infection_results <- results(dds_retro_LRT)
sig_retro_genes_infection <- subset(dds_retro_LRT_infection_results, padj < 0.05)

# get all sig genes for retro LRT
all_sig_retro_gene_names_LRT <- rownames(sig_retro_genes_infection)

# combine DESeq tables for LRTs in retro/human genes
all_sig_genes_infection <- rbind(sig_human_genes_infection, sig_retro_genes_infection) %>%
  as.data.frame

# get counts of significant genes identified by the LRT
sig_human_genes_infection_expr <- non_zero_human_counts_scaled[rownames(sig_human_genes_infection),]
sig_retro_genes_infection_expr <- non_zero_retro_counts_scaled[rownames(sig_retro_genes_infection),]

# function to make p-value histogram from DESeq results
make_pvalue_histogram <- function(DESeq_results, title, save = FALSE){
  
  df <- as.data.frame(DESeq_results)
  
  plot <- ggplot(df, aes(x = padj)) + geom_histogram() + 
    theme_minimal() +
    labs(x = '\nAdjusted p-value', title = title) +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'), text = element_text(size = 14))
  
  if(save){
    filepath = paste('figures/', title, '.png', sep = '')
    ggsave(filepath, width = 33.8, height = 19.5, units = 'cm')
  }
  
  return(plot)
}

# p-value histograms
dds_human_LRT_pval_plot <- make_pvalue_histogram(dds_human_LRT_infection_results, 
                                                 title = 'p-value histogram for the likelihood ratio test for SARS-MERS-SARSCov2 infection (human genes)', 
                                                 save = TRUE)

dds_retro_LRT_pval_plot <- make_pvalue_histogram(dds_retro_LRT_infection_results, 
                                                 title = 'p-value histogram for the likelihood ratio test for SARS-MERS-SARSCov2 infection (retro genes)', 
                                                 save = TRUE)