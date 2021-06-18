# Arguments for this script are:
# --clusterdata : filepath to the cluster data object ClusterData.Rdata outputted by clustering.R
# --whichclusters : space-separated values indicating which clusters enriched terms to be plotted for (e.g. --whichclusters 1 2 3 4)
# --outdir : directory to save data files (this argument can be omitted, in which case data files will be saved to the current directory)

options(warn = -1)

source('utils.R')

reqs <- c('RColorBrewer', 'dplyr', 'stringr')
get_reqs(reqs)
bioc_reqs <- c('clusterProfiler', 'org.Hs.eg.db')
get_bioc_reqs(bioc_reqs)

suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(org.Hs.eg.db))

args <- commandArgs(trailingOnly = T)
parsed_args <- parseArgs(args)
names(parsed_args) <- stringr::str_to_lower(names(parsed_args))
req_args <- c('clusterdata', 'whichclusters')

checkArgs(parsed_args, req_args, file_args = c('clusterdata'))

load(parsed_args['clusterdata'])

# get clusters to plot from parsed arguments
clusters_to_plot <- stringr::str_trim(unlist(stringr::str_split(parsed_args['whichclusters'], ' ')))

# check that all specified clusters are present in the data
if(any(!(clusters_to_plot %in% clusterdata[['expr']]$cluster))){
  clusters_not_present <- clusters_to_plot[!(clusters_to_plot %in% clusterdata[['expr']]$cluster)]
  clusters_not_present <- paste(clusters_not_present, collapse = ', ')
  stop(paste('The following clusters were not found in the cluster data file provided:', clusters_not_present))
}

# import gene ID database
human_IDs <- org.Hs.eg.db

# function for getting either enriched GO or KEGG terms
get_enriched_terms <- function(gene_symbols, terms = c('GO', 'KEGG')){
  
  # Function that takes in a list of gene symbols and performs enrichment analysis
  # 
  # Arguments
  # ------------
  # 
  # gene_symbols : vector of gene names (symbols) to be analysed for enrichment
  # terms : type of gene annotations to be used for enrichment analysis, either GO or KEGG
  # 
  # 
  # Returns
  # ------------
  # 
  # enriched : S4 object containing enriched terms for gene_symbols
  
  # convert gene symbols to entrez IDs
  entrez_ids <- mapIds(human_IDs, gene_symbols, column = 'ENTREZID', keytype = 'SYMBOL')
  # remove NAs
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  if(terms == 'GO'){
    enriched <- clusterProfiler::enrichGO(entrez_ids, 'org.Hs.eg.db', ont = 'BP', qvalueCutoff = 0.1)
  }else{
    enriched <- clusterProfiler::enrichKEGG(entrez_ids, organism = 'hsa', qvalueCutoff = 0.1)
  }
  
  return(enriched)
  
}

# get enriched GO annotations for each cluster
enriched_GO_annotations <- lapply(unique(clusterdata[['expr']]$cluster), function(cluster){
  
  # isolate human genes in the cluster
  cluster_mask <- clusterdata[['expr']]$cluster == cluster
  human_gene_mask <- rownames(clusterdata[['expr']]) %in% rownames(subset(clusterdata[['expr']], gene_type == 'Human'))
  in_cluster <- clusterdata[['expr']][cluster_mask & human_gene_mask,]
  genes <- rownames(in_cluster)
  
  # get enriched gene ontologies for the cluster genes
  enriched_GO_terms <- get_enriched_terms(genes, terms = 'GO')
  enriched_GO_terms@result$cluster <- cluster
  return(enriched_GO_terms)
  
})

# set names of the list
names(enriched_GO_annotations) <- unique(clusterdata[['expr']]$cluster)

# mask for filtering results without any terms
no_sig_GO_found_mask <- sapply(enriched_GO_annotations, function(results){nrow(results@result) == 0})

# filter for clusters with enriched GO annotations found
enriched_GO_annotations <- enriched_GO_annotations[!no_sig_GO_found_mask]

# add new term column to each GO enrichment df
enriched_GO_annotations <- lapply(enriched_GO_annotations, function(ann) {
  # replace description for GO:0150115 - labelled NA
  ann@result$Description <- ifelse(ann@result$ID == 'GO:0150115', 'cell-substrate junction organisation', ann@result$Description)
  ann@result$term <- paste(ann@result$ID, ann@result$Description, sep = ': ')
  return(ann)
})

# combine enriched GO annotations into a single dataframe, write to file
enriched_GO_annotations_combined <- do.call(rbind, lapply(enriched_GO_annotations, function(ann) ann@result[1:10,]))
enriched_GO_annotations_to_export <- enriched_GO_annotations_combined[,c('term', 'GeneRatio', 'BgRatio', 'pvalue', 'p.adjust', 'qvalue', 'cluster')]

if('outdir' %in% names(parsed_args)){
  cluster_GO_filepath <- paste0(parsed_args['outdir'], '/enriched_GO_annotations.csv')
} else {
  cluster_GO_filepath <- 'enriched_GO_annotations.csv'
}

writeLines(paste('Saving cluster-enriched GO annotations to:', cluster_GO_filepath))

write.csv(enriched_GO_annotations_to_export, cluster_GO_filepath, row.names = F)

# get enriched KEGG annotations
enriched_KEGG_annotations <- lapply(unique(clusterdata[['expr']]$cluster), function(cluster){
  
  # isolate human genes in the cluster
  cluster_mask <- clusterdata[['expr']]$cluster == cluster
  human_gene_mask <- rownames(clusterdata[['expr']]) %in% rownames(subset(clusterdata[['expr']], gene_type == 'Human'))
  in_cluster <- clusterdata[['expr']][cluster_mask & human_gene_mask,]
  genes <- rownames(in_cluster)
  
  # get enriched KEGG pathways for the cluster genes
  enriched_KEGG_terms <- get_enriched_terms(genes, terms = 'KEGG')
  enriched_KEGG_terms@result$cluster <- cluster
  return(enriched_KEGG_terms)
  
})

# set names of the list
names(enriched_KEGG_annotations) <- unique(clusterdata[['expr']]$cluster)

# mask for filtering results without any terms
no_sig_KEGG_found_mask <- sapply(enriched_KEGG_annotations, function(results){nrow(results@result) == 0})

# filter for clusters with enriched GO annotations found
enriched_KEGG_annotations <- enriched_KEGG_annotations[!no_sig_KEGG_found_mask]

# add new term column to each KEGG enrichment df
enriched_KEGG_annotations <- lapply(enriched_KEGG_annotations, function(ann) {
  ann@result$term <- paste(ann@result$ID, ann@result$Description, sep = ': ')
  return(ann)
})

# combine enriched GO annotations into a single dataframe, write to file
enriched_KEGG_annotations_combined <- do.call(rbind, lapply(enriched_KEGG_annotations, function(ann) ann@result[1:10,]))
enriched_KEGG_annotations_to_export <- enriched_KEGG_annotations_combined[,c('term', 'GeneRatio', 'BgRatio', 'pvalue', 'p.adjust', 'qvalue', 'cluster')]

if('outdir' %in% names(parsed_args)){
  cluster_KEGG_filepath <- paste0(parsed_args['outdir'], '/enriched_GO_annotations.csv')
} else {
  cluster_KEGG_filepath <- 'enriched_GO_annotations.csv'
}

writeLines(paste('Saving cluster-enriched KEGG annotations to:', cluster_KEGG_filepath))

write.csv(enriched_KEGG_annotations_to_export, cluster_KEGG_filepath, row.names = F)

# create figures directory if it does not already exist
if(!dir.exists('figures')){
  dir.create('figures')
}

# get top n GO annotations for each cluster
n <- 5 
top_GO_annotations <- lapply(clusters_to_plot, function(cluster){
  cluster_annotations_GO <- enriched_GO_annotations[[cluster]]@result
  top_n_annotations <- cluster_annotations_GO[order(cluster_annotations_GO$p.adjust)[1:n],]
  top_n_annotations$cluster <- cluster
  top_n_annotations$logp <- -log10(top_n_annotations$pvalue)
  return(top_n_annotations)
}) %>% 
  do.call(rbind, .)

# get enriched GO terms
GO_enrichment_matrix <- matrix(0, ncol = length(clusters_to_plot), nrow = length(top_GO_annotations$term))
rownames(GO_enrichment_matrix) <- top_GO_annotations$term
colnames(GO_enrichment_matrix) <- clusters_to_plot
for(term in unique(rownames(GO_enrichment_matrix))){
  
  term_mask <- which(rownames(GO_enrichment_matrix) == term)
  logps <- sapply(clusters_to_plot, function(cluster){
    ann <- enriched_GO_annotations[[cluster]]@result
    if(term %in% ann$term){
      if(ann[ann$term == term,'pvalue'] < 0.05){
        logp <- -log10(ann[ann$term == term,'pvalue'])
        return(logp)
      } else {
        return(0)
      }
    } else {
      return(0)
    }
  })
  for(i in term_mask){
    GO_enrichment_matrix[i,] <- logps
  }
}

enriched_GO_heatmap <- pheatmap::pheatmap(GO_enrichment_matrix, cluster_rows = FALSE, 
                                          color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
                                          cluster_cols = FALSE, angle_col = 0, 
                                          filename = 'figures/enriched_GO_terms.png', 
                                          width = 9, height = 6, cellwidth = 25, cellheight = 15)

# get top n KEGG annotations for each cluster
n <- 5 
top_KEGG_annotations <- lapply(clusters_to_plot, function(cluster){
  cluster_annotations_KEGG <- enriched_KEGG_annotations[[cluster]]@result
  n <- min(n, nrow(subset(cluster_annotations_KEGG, p.adjust < 0.05)))
  top_n_annotations <- cluster_annotations_KEGG[order(cluster_annotations_KEGG$p.adjust)[1:n],]
  top_n_annotations$cluster <- cluster
  top_n_annotations$logp <- -log10(top_n_annotations$p.adjust)
  return(top_n_annotations)
}) %>% 
  do.call(rbind, .)

KEGG_enrichment_matrix <- matrix(0, ncol = length(clusters_to_plot), nrow = length(top_KEGG_annotations$term))
rownames(KEGG_enrichment_matrix) <- top_KEGG_annotations$term
colnames(KEGG_enrichment_matrix) <- clusters_to_plot
for(term in unique(rownames(KEGG_enrichment_matrix))){
  
  term_mask <- which(rownames(KEGG_enrichment_matrix) == term)
  logps <- sapply(clusters_to_plot, function(cluster){
    ann <- enriched_KEGG_annotations[[cluster]]@result
    if(term %in% ann$term){
      if(ann[ann$term == term,'pvalue'] < 0.05){
        logp <- -log10(ann[ann$term == term,'pvalue'])
        return(logp)
      } else {
        return(0)
      }
    } else {
      return(0)
    }
  })
  for(i in term_mask){
    KEGG_enrichment_matrix[i,] <- logps
  }
}

enriched_KEGG_heatmap <- pheatmap::pheatmap(KEGG_enrichment_matrix, cluster_rows = FALSE, 
                                            color = colorRampPalette(brewer.pal(n = 7, name = "Oranges"))(100),
                                            cluster_cols = FALSE, angle_col = 0, 
                                            filename = 'figures/enriched_KEGG_terms.png',  
                                            width = 7, height = 5, cellwidth = 25, cellheight = 15)