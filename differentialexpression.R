# Arguments for this script are:
# --human_counts : filepath to human count matrix (tab-separated)
# --TE_counts : filepath to TE count matrix (tab-separated)
# --metadata : filepath to metadata file (comma-separated), with three columns: "sample", "condition", and "batch" (batch can be omitted)
# --compare : filepath to compare file (comma-separted), with two columns "A" and "B" for comparisons to be performed by DESeq2 (A = control, B = treatment)
# --outdir : directory to save data files (this argument can be omitted, in which case data files will be saved to the current directory)

options(warn = -1)

source('utils.R')

args <- commandArgs(trailingOnly = T)
parsed_args <- parseArgs(args)
names(parsed_args) <- stringr::str_to_lower(names(parsed_args))
req_args <- c('human_counts', 'te_counts', 'metadata', 'compare')

checkArgs(parsed_args, req_args, file_args = req_args)

# if output directory does not exist, create
if('outdir' %in% names(parsed_args)){
  if(!dir.exists(parsed_args['outdir'])){
    dir.create(parsed_args['outdir'])
  }
}

reqs <- c('ggplot2', 'dplyr', 'stringr')
get_reqs(reqs)
bioc_reqs <- c('DESeq2', 'org.Hs.eg.db')
get_bioc_reqs(bioc_reqs)

suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))

# import human count data
human_counts <- read.table(parsed_args['human_counts'], header = TRUE)

# import retroelement count data
retro_counts <- read.table(parsed_args['te_counts'], header = TRUE)

# remove retro genes with zero counts
all_zero_retro_mask <- apply(retro_counts, 1, function(row){all(row == 0)})
non_zero_retro_counts <- retro_counts[!all_zero_retro_mask,]

# import metadata
metadata <- read.csv(parsed_args['metadata'], sep = ',', header = TRUE)
colnames(metadata) <- str_to_lower(colnames(metadata))

if(!('condition' %in% colnames(metadata))){
  stop('Condition column not found in metadata file.')
}

if(!('sample' %in% colnames(metadata))){
  stop('Sample column not found in metadata file.')
}

rownames(metadata) <- metadata$sample
metadata$sample <- NULL
metadata$condition <- as.factor(metadata$condition)

if(!setequal(rownames(metadata), colnames(human_counts))){
  stop('Sample column of metadata file does not match column names of human count matrix. Try again.')
}

if(!setequal(rownames(metadata), colnames(retro_counts))){
  stop('Sample column of metadata file does not match column names of TE count matrix. Try again.')
}

if('batch' %in% colnames(metadata)){
  metadata$batch <- as.factor(metadata$batch)
}

# import comparison sheet
compare <- read.csv(parsed_args['compare'], sep = ',', header = TRUE)
colnames(compare) <- str_to_lower(colnames(compare))

# check all entries of comparison file are in the metadata file
if(!all(unlist(compare) %in% metadata$condition)){
  compare_conditions <- unlist(compare)
  conditions_not_found <- paste(compare_conditions[!(compare_conditions %in% metadata$condition)], collapse = ', ')
  stop(paste('The following entries in compare file could not be found in the conditon column of the metadata file: ', conditions_not_found))
}

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

if('batch' %in% colnames(metadata)){
  # build DESeq model for human genes
  dds_human <- DESeq2::DESeqDataSetFromMatrix(non_zero_human_counts, metadata, 
                                              design = ~condition+batch)
  # build DESeq model for retro genes
  dds_retro <- DESeq2::DESeqDataSetFromMatrix(non_zero_retro_counts, metadata, 
                                              design = ~condition+batch)
} else {
  dds_human <- DESeq2::DESeqDataSetFromMatrix(non_zero_human_counts, metadata, 
                                              design = ~condition)
  dds_retro <- DESeq2::DESeqDataSetFromMatrix(non_zero_retro_counts, metadata, 
                                              design = ~condition)
}

# get regularised human gene counts
non_zero_human_counts_reg <- DESeq2::varianceStabilizingTransformation(dds_human, blind = FALSE) %>% assay

# filter for highly-expressed genes
dds_retro <- DESeq2::estimateSizeFactors(dds_retro)
highly_expr_retro <- rowSums(counts(dds_retro, normalized=TRUE) >= 5) >= 3
dds_retro <- dds_retro[highly_expr_retro,]

# get regularised retro gene counts
non_zero_retro_counts_reg <- DESeq2::varianceStabilizingTransformation(dds_retro, blind = FALSE) %>% assay

# function for removing batch effects using OLS fits with a batch variable
remove_batch_effects <- function(expr, coldata, condition_col = 'condition', batch_col = 'batch'){
  
  # make data matrices for all genes
  Xs <- split(expr, 1:nrow(expr)) %>% 
    lapply(function(row) cbind(exp = row, coldata)) 
  
  # make linear model formulas with column names
  model_formula <- as.formula(paste('exp ~', condition_col, '+', batch_col))
  
  # generate linear models
  linear_models <- lapply(Xs, lm, formula = model_formula)
  
  # get the vector of batch memberships
  batch_levels <- sort(unique(coldata[,batch_col]))
  batch_vector <- coldata[,batch_col]
  
  # get the batch effects in a matrix
  batch_effects <- lapply(linear_models, function(model){
    coefs <- model$coefficients[paste(batch_col, batch_levels[2:length(batch_levels)], sep='')] 
    coefs <- c(0, coefs)
    names(coefs) <- batch_levels
    effect <- coefs[batch_vector]
    return(effect)
  }) %>%
    do.call(rbind, .)
  
  # subtract the batch effect from the raw expression values
  batch_adjusted_expr <- expr - batch_effects
  
  return(batch_adjusted_expr)
  
}

if('batch' %in% colnames(metadata)){
  # remove batch effects from human genes and retroelements
  non_zero_human_counts_reg <- remove_batch_effects(non_zero_human_counts_reg, metadata)
  non_zero_retro_counts_reg <- remove_batch_effects(non_zero_retro_counts_reg, metadata)
}

non_zero_human_counts_scaled <- t(scale(t(non_zero_human_counts_reg)))
non_zero_retro_counts_scaled <- t(scale(t(non_zero_retro_counts_reg)))

# function for making PCA plot from a transcript count table
plot_pca <- function(expr, title, coldata, condition_col = 'condition', batch_col = 'batch'){
  
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
  plot(dend)
}

if(!dir.exists('figures')){
  dir.create('figures')
}

PCA_plot_human <- plot_pca(non_zero_human_counts_reg, title = 'Human genes', coldata = metadata)
ggsave('figures/PCA_sample_plot_human.png', PCA_plot_human)

PCA_plot_retro <- plot_pca(non_zero_retro_counts_reg, title = 'Retro genes', coldata = metadata) 
ggsave('figures/PCA_sample_plot_retro.png', PCA_plot_retro)

# get cluster dendrograms
png('figures/sample_dendrogram_human_genes.png', width = 1280, height = 720, units = 'px')
plot_dendrogram(non_zero_human_counts_reg, title = 'Cluster dendrogram - based on human genes')
dev.off()

png('figures/sample_dendrogram_retro_genes.png', width = 1280, height = 720, units = 'px')
plot_dendrogram(non_zero_retro_counts_reg, title = 'Cluster dendrogram - based on retro genes')
dev.off()

writeLines('PCA and dendrogram plots produced, saved to figures/\n')

counts <- list('human' = non_zero_human_counts,
               'TE' = non_zero_retro_counts,
               'human_reg' = non_zero_human_counts_scaled,
               'TE_reg' = non_zero_retro_counts_scaled)

if('outdir' %in% names(parsed_args)){
  counts_filepath <- paste0(parsed_args['outdir'], '/counts.Rdata')
} else {
  counts_filepath <- 'counts.Rdata'
}

writeLines(paste('Saving counts object to:', counts_filepath, '\n'))

save(counts, file = counts_filepath)

writeLines('Preprocessing finished, starting differential expression analysis...\n')

# Wald tests to identify significant human genes for SARS vs mock, IAV vs mock, and SARS vs IAV
dds_human_wald <- DESeq(dds_human, test = 'Wald', fitType = 'parametric') 
dds_retro_wald <- DESeq(dds_retro, test = 'Wald', fitType = 'parametric')

human_DE_results <- list()
TE_DE_results <- list()

# conduct differential expression analysis for all comparisons specified in the compare sheet
for(i in 1:nrow(compare)){
  
  comparison <- paste(c(compare[i, 1], compare[i, 2]), collapse = '_vs_')
  
  # differential expression for human genes
  DE_human <- results(dds_human_wald, contrast = c('condition', compare[i, 1], compare[i, 2]))
  DE_human_p_values <- make_pvalue_histogram(DE_human, title = comparison, save = TRUE)
  volcano_human <- make_volcano_plot(DE_human, title = comparison, transcript_type = 'Human genes', save = TRUE)
  human_DE_results[[comparison]] <- DE_human
  
  # differential expression for transposable elements
  DE_retro <- results(dds_retro_wald, contrast = c('condition', compare[i, 1], compare[i, 2]))
  DE_retro_p_values <- make_pvalue_histogram(DE_retro, title = comparison, save = TRUE)
  volcano_TE <- make_volcano_plot(DE_retro, title = comparison, transcript_type = 'Transposable elements', save = TRUE)
  TE_DE_results[[comparison]] <- DE_retro
  
}

# perform LRT
if('batch' %in% colnames(metadata)){
  LRT_human <- DESeq(dds_human, test = 'LRT', reduced = ~batch)
  LRT_retro <- DESeq(dds_retro, test = 'LRT', reduced = ~batch)
} else {
  LRT_human <- DESeq(dds_human, test = 'LRT', reduced = ~1)
  LRT_retro <- DESeq(dds_retro, test = 'LRT', reduced = ~1)
}

LRT_human_results <- results(LRT_human)
human_DE_results[['LRT']] <- LRT_human_results

LRT_retro_results <- results(LRT_retro)
TE_DE_results[['LRT']] <- LRT_retro_results

writeLines('Finished differential expression analysis.\n')

if('outdir' %in% names(parsed_args)){
  DE_filepath <- paste0(parsed_args['outdir'], '/DiffExp.Rdata')
} else {
  DE_filepath <- 'DiffExp.Rdata'
}

writeLines(paste('Saving differential expression object to:', DE_filepath, '\n'))

DE_results <- list('human' = human_DE_results,
                   'TE' = TE_DE_results)

save(DE_results, file = DE_filepath)