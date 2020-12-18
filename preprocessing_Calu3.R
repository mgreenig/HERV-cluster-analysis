library(org.Hs.eg.db)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(dendextend)

SARSCov2_Calu3_human <- read.table('data/Gene_Counts_Covid19_Calu3_4Jez.txt', sep = "\t", header = TRUE)
MERS_SARS_Calu3_human <- read.table('data/Gene_Counts_MERSSARS_Calu3_4Jez.txt', sep = "\t", header = TRUE)

SARSCov2_Calu3_retro <- read.table('data/Retro_Counts_Covid19_Calu3_4Jez.txt', sep = "\t", header = TRUE)
MERS_SARS_Calu3_retro <- read.table('data/Retro_Counts_MERSSARS_Calu3_4Jez.txt', sep = "\t", header = TRUE)

SARSCov2_Calu3_metadata <- read.table('data/Samples_Covid19_Calu3_4Jez.txt', sep = "\t", header = TRUE)
MERS_SARS_Calu3_metadata <- read.table('data/Samples_MERSSARS_Calu3_4Jez.txt', sep = "\t", header = TRUE)

# import human count data
human_counts <- cbind(SARSCov2_Calu3_human, MERS_SARS_Calu3_human)

# import retroelement count data
shared_retroelements <- rownames(SARSCov2_Calu3_retro)[rownames(SARSCov2_Calu3_retro) %in% rownames(MERS_SARS_Calu3_retro)]
retro_counts <- cbind(SARSCov2_Calu3_retro[shared_retroelements,], 
                      MERS_SARS_Calu3_retro[shared_retroelements,])

# remove retro genes with zero counts
all_zero_retro_mask <- apply(retro_counts, 1, function(row){all(row == 0)})
non_zero_retro_counts <- retro_counts[!all_zero_retro_mask,]

# get metadata
SARSCov2_Calu3_metadata$Batch <- 0
MERS_SARS_Calu3_metadata$Batch <- 1
metadata <- rbind(SARSCov2_Calu3_metadata[,c('Batch', 'Infection')], MERS_SARS_Calu3_metadata[,c('Batch', 'Infection')])
metadata$Infection <- as.character(metadata$Infection)
metadata$Infection[grep('MOCK', metadata$Infection)] <- 'MOCK'
metadata$Infection <- factor(metadata$Infection, levels = c('MOCK', 'SARS_Cov2', 'SARS', 'MERS'))
metadata$Batch <- factor(metadata$Batch)

# import gene ID database
human_IDs <- org.Hs.eg.db

# remove suffix from ensembl gene names
human_genes_cleaned <- rownames(human_counts) %>% 
  stringr::str_remove_all('\\..*$') %>% 
  unique

# map ensembl IDs to gene symbols
ids_mapped <- mapIds(human_IDs, human_genes_cleaned, column = 'SYMBOL', keytype='ENSEMBL')
human_gene_symbols <- ids_mapped[!is.na(ids_mapped)]
unique_human_genes <- make.unique(human_gene_symbols, sep = '_')
names(unique_human_genes) <- names(human_gene_symbols)

# get indices of mapped ensembl ids
unique_human_gene_idx <- match(names(unique_human_genes), human_genes_cleaned)

# isolate annotated genes and change rownames
annotated_human_counts <- human_counts[unique_human_gene_idx,]
rownames(annotated_human_counts) <- unname(unique_human_genes)

# remove human genes with zero counts
all_zero_human_mask <- apply(annotated_human_counts, 1, function(row){all(row == 0)})
non_zero_human_counts <- annotated_human_counts[!all_zero_human_mask,]

# build DESeq model for human genes
dds_human <- DESeq2::DESeqDataSetFromMatrix(non_zero_human_counts, metadata, design = ~Infection)

# get regularised human gene counts
non_zero_human_counts_reg <- DESeq2::varianceStabilizingTransformation(dds_human, blind = FALSE) %>% assay
# scale to unit variance per-gene using transpose
non_zero_human_counts_scaled <- t(scale(t(non_zero_human_counts_reg)))

# build DESeq model for retro genes
dds_retro <- DESeq2::DESeqDataSetFromMatrix(non_zero_retro_counts, metadata, design = ~Infection)

# filter for highly-expressed genes
dds_retro <- DESeq2::estimateSizeFactors(dds_retro)
idx <- rowSums(counts(dds_retro, normalized=TRUE) >= 5) >= 3
dds_retro <- dds_retro[idx,]

# get regularised retro gene counts
non_zero_retro_counts_reg <- DESeq2::varianceStabilizingTransformation(dds_retro, blind = FALSE) %>% assay
# scale to unit variance per-gene using transpose
non_zero_retro_counts_scaled <- t(scale(t(non_zero_retro_counts_reg)))

# function for making PCA plot from a transcript count table
plot_pca <- function(expr, title, coldata, condition_col = 'Infection', batch_col = 'Batch'){
  
  # run PCA and extract first two PCs
  PCA_results <- prcomp(expr, scale.=TRUE)
  PC1_PC2 <- PCA_results$rotation[,c('PC1', 'PC2')] %>% as.data.frame
  PC1_PC2$condition <- coldata[[condition_col]]
  PC1_PC2$batch <- coldata[[batch_col]]
  
  # get pct of total variance for first two PCs
  PC1_PC2_variance_pcts <- PCA_results$sdev[1:2]**2 / sum(PCA_results$sdev**2) * 100
  
  PC1_PC2_human_plot <- ggplot(PC1_PC2, aes(x = PC1, y = PC2, colour = condition, shape = batch)) + 
    geom_point(size = 3) + 
    labs(x = paste('PC1 (', PC1_PC2_variance_pcts[1] %>% round(2), '%)', sep = ''), 
         y =  paste('PC2 (', PC1_PC2_variance_pcts[2] %>% round(2), '%)', sep = ''),
         title = title)
  
  return(PC1_PC2_human_plot)
  
}


# function for making dendrogram from a transcript count table
plot_dendrogram <- function(expr, title){
  
  # get distance matrix
  samples_dist <- expr %>% 
    t %>%
    dist

  # cluster using average linkage
  samples_clustering <- hclust(samples_dist, method = 'average')
  dend <- as.dendrogram(samples_clustering)
  
  # return plot
  dend_labels <- labels(dend)
  labels(dend) <- ""
  plot(dend)
  text(x = 1:length(dend_labels), labels = dend_labels, srt = 45, adj = c(1,1), xpd = T)

}

# if run from the command line, make plots
if(!interactive()){
  
  PCA_plot_human <- plot_pca(non_zero_human_counts_reg, title = 'Human genes', coldata = metadata)
  ggsave(paste('figures/PCA_sample_plot_human_', cell_line, '.png', sep = ''), PCA_plot_human)
  
  PCA_plot_retro <- plot_pca(non_zero_retro_counts_reg, title = 'Retro genes', coldata = metadata) 
  ggsave(paste('figures/PCA_sample_plot_retro_', cell_line, '.png', sep = ''), PCA_plot_retro)
  
  # get cluster dendrograms
  png(paste('figures/sample_dendrogram_human_genes_', cell_line, '.png', sep = ''), width = 1280, height = 720, units = 'px')
  plot_dendrogram(non_zero_human_counts_reg, title = 'Cluster dendrogram - based on human genes')
  dev.off()
  
  png(paste('figures/sample_dendrogram_retro_genes_', cell_line, '.png', sep = ''), width = 1280, height = 720, units = 'px')
  plot_dendrogram(non_zero_retro_counts_reg, title = 'Cluster dendrogram - based on retro genes')
  dev.off()

}