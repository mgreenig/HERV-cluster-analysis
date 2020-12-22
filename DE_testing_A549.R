suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(gridExtra))

source('preprocessing_A549.R')

# Wald tests to identify significant human genes for SARS vs mock, IAV vs mock, and SARS vs IAV
dds_human_wald <- DESeq(dds_human, test = 'Wald', fitType = 'parametric') 
  
dds_human_Cov2_results <- results(dds_human_wald, name = 'Infection_SARS_Cov2_vs_MOCK')
sig_human_genes_Cov2 <- subset(dds_human_Cov2_results, padj < 0.05 & abs(log2FoldChange) > 1)

dds_human_IAV_results <- results(dds_human_wald, name = 'Infection_IAV_vs_MOCK')
sig_human_genes_IAV <- subset(dds_human_IAV_results, padj < 0.05 & abs(log2FoldChange) > 1)

dds_human_Cov2_IAV_results <- results(dds_human_wald, contrast = c('Infection', 'SARS_Cov2', 'IAV'))
sig_human_genes_Cov2_IAV <- subset(dds_human_Cov2_IAV_results, padj < 0.05 & abs(log2FoldChange) > 1)

# get all sig genes for human wald tests
all_sig_human_gene_names_wald <- c(rownames(sig_human_genes_Cov2), 
                                   rownames(sig_human_genes_IAV), 
                                   rownames(sig_human_genes_Cov2_IAV)) %>% unique

# LRT to get genes for clustering
dds_human_LRT <- DESeq(dds_human, test = 'LRT', reduced = ~Batch)
dds_human_LRT_infection_results <- results(dds_human_LRT)
sig_human_genes_infection <- subset(dds_human_LRT_infection_results, padj < 0.05)

# get all sig genes for human LRT
all_sig_human_gene_names_LRT <- rownames(sig_human_genes_infection)

# Wald tests to identify significant retro genes for SARS vs mock, IAV vs mock, and SARS vs IAV
dds_retro_wald <- DESeq(dds_retro, test = 'Wald', fitType = 'parametric')

dds_retro_Cov2_results <- results(dds_retro_wald, name = 'Infection_SARS_Cov2_vs_MOCK')
sig_retro_genes_Cov2 <- subset(dds_retro_Cov2_results, padj < 0.05 & abs(log2FoldChange) > 1)

dds_retro_IAV_results <- results(dds_retro_wald, name = 'Infection_IAV_vs_MOCK')
sig_retro_genes_IAV <- subset(dds_retro_IAV_results, padj < 0.05 & abs(log2FoldChange) > 1)

dds_retro_Cov2_IAV_results <- results(dds_retro_wald, contrast = c('Infection', 'SARS_Cov2', 'IAV'))
sig_retro_genes_Cov2_IAV <- subset(dds_retro_Cov2_IAV_results, padj < 0.05 & abs(log2FoldChange) > 1)

# get all sig genes for retro wald tests
all_sig_retro_gene_names_wald <- c(rownames(sig_retro_genes_Cov2), 
                                   rownames(sig_retro_genes_IAV), 
                                   rownames(sig_retro_genes_Cov2_IAV)) %>% unique

# LRT to get genes for clustering
dds_retro_LRT <- DESeq(dds_retro, test = 'LRT', reduced = ~Batch)
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
dds_human_Cov2_pval_plot <- make_pvalue_histogram(dds_human_Cov2_results, 
                                                  title = paste('pvalue_histogram_Cov2_mock_human_', cell_line, sep = ''), save = TRUE)
dds_human_IAV_pval_plot <- make_pvalue_histogram(dds_human_IAV_results, 
                                                 title = paste('pvalue_histogram_IAV_mock_human_', cell_line, sep = ''), save = TRUE)
dds_human_Cov2_IAV_pval_plot <- make_pvalue_histogram(dds_human_Cov2_IAV_results, 
                                                      title = paste('pvalue_histogram_Cov2_IAV_human_', cell_line, sep = ''), save = TRUE)

dds_human_LRT_pval_plot <- make_pvalue_histogram(dds_human_LRT_infection_results, 
                                                 title = paste('pvalue_histogram_LRT_human_', cell_line, sep = ''), save = TRUE)

dds_retro_Cov2_pval_plot <- make_pvalue_histogram(dds_retro_Cov2_results, 
                                                  title = paste('pvalue_histogram_Cov2_mock_retro_', cell_line, sep = ''), save = TRUE)
dds_retro_IAV_pval_plot <- make_pvalue_histogram(dds_retro_IAV_results, 
                                                 title = paste('pvalue_histogram_IAV_mock_retro_', cell_line, sep = ''), save = TRUE)
dds_retro_Cov2_IAV_pval_plot <- make_pvalue_histogram(dds_retro_Cov2_IAV_results, 
                                                      title = paste('pvalue_histogram_Cov2_IAV_retro_', cell_line, sep = ''), save = TRUE)

dds_retro_LRT_pval_plot <- make_pvalue_histogram(dds_retro_LRT_infection_results, 
                                                 title = paste('pvalue_histogram_LRT_retro_', cell_line, sep = ''), save = TRUE)

# function to make volcano plot from DESeq results
make_volcano_plot <- function(DESeq_results, title, transcript_type = c('Human genes', 'Retroelements')){
  
  volcano_plot <- EnhancedVolcano::EnhancedVolcano(DESeq_results, lab = rownames(DESeq_results), 
                                                   x = 'log2FoldChange', y = 'padj', caption = '', 
                                                   title = title, transcriptLabSize = 7, 
                                                   subtitle = transcript_type, titleLabSize = 24,
                                                   subtitleLabSize = 22) + xlim(-20, 20)
  
  return(volcano_plot)
  
}

if(!interactive()){
  
  # make volcano plots for all pairwise comparisons
  dds_human_Cov2_volcano_plot <- make_volcano_plot(dds_human_Cov2_results, 
                                                   title = 'SARS-Cov2 vs wild-type', 
                                                   transcript_type = 'Human genes')
  
  dds_human_IAV_volcano_plot <- make_volcano_plot(dds_human_IAV_results, 
                                                  title = 'Influenza A vs wild-type', 
                                                  transcript_type = 'Human genes')
  
  dds_human_Cov2_IAV_volcano_plot <- make_volcano_plot(dds_human_Cov2_IAV_results, 
                                                       title = 'SARS-Cov2 vs Influenza A', 
                                                       transcript_type = 'Human genes')
  
  dds_retro_Cov2_volcano_plot <- make_volcano_plot(dds_retro_Cov2_results, 
                                                   title = 'SARS-Cov2 vs wild-type', 
                                                   transcript_type = 'Retroelements')
  
  dds_retro_IAV_pval_plot <- make_volcano_plot(dds_retro_IAV_results, 
                                               title = 'Influenza A vs wild-type', 
                                               transcript_type = 'Retroelements')
  
  dds_retro_Cov2_IAV_pval_plot <- make_volcano_plot(dds_retro_Cov2_IAV_results, 
                                                    title = 'SARS-Cov2 vs Influenza A', 
                                                    transcript_type = 'Retroelements')
  
  all_volcano_plots <- grid.arrange(dds_human_Cov2_volcano_plot,
                                    dds_human_IAV_volcano_plot, 
                                    dds_human_Cov2_IAV_volcano_plot, 
                                    dds_retro_Cov2_volcano_plot,
                                    dds_retro_IAV_pval_plot,
                                    dds_retro_Cov2_IAV_pval_plot, 
                                    ncol = 3, nrow = 2)
  
  ggsave('figures/volcano_plots.png', plot = all_volcano_plots, width = 50.8, height = 28.6, units = 'cm')
}
