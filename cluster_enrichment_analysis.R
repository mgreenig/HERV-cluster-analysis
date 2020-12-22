library(clusterProfiler)
library(RColorBrewer)

cell_line <- commandArgs(trailingOnly = TRUE)

source('full_clustering.R')

# function for getting either enriched GO or KEGG terms
get_enriched_terms <- function(gene_symbols, terms = c('GO', 'KEGG')){
  
  # Function that takes in a list of gene symbols and performs enrichment analysis
  # 
  # Arguments
  # ------------
  # 
  # gene_symbols : list of gene names (symbols) to be analysed for enrichment
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
    enriched <- enrichGO(entrez_ids, 'org.Hs.eg.db', ont = 'BP', qvalueCutoff = 0.1)
  }else{
    enriched <- enrichKEGG(entrez_ids, organism = 'hsa', qvalueCutoff = 0.1)
  }
  
  return(enriched)
  
}

# get enriched GO annotations for each cluster
enriched_GO_annotations <- lapply(unique(all_sig_genes_infection$cluster), function(cluster){
  
  # isolate human genes in the cluster
  cluster_mask <- all_sig_genes_infection$cluster == cluster
  human_gene_mask <- rownames(all_sig_genes_infection) %in% all_sig_human_gene_names_LRT
  in_cluster <- all_sig_genes_infection[cluster_mask & human_gene_mask,]
  genes <- rownames(in_cluster)
  
  # get enriched gene ontologies for the cluster genes
  enriched_GO_terms <- get_enriched_terms(genes, terms = 'GO')
  return(enriched_GO_terms)
  
})

# set names of the list
names(enriched_GO_annotations) <- unique(all_sig_genes_infection$cluster)

# mask for filtering results without any terms
no_sig_GO_found_mask <- sapply(enriched_GO_annotations, function(results){nrow(results@result) == 0})

# filter for clusters with enriched GO annotations found
enriched_GO_annotations <- enriched_GO_annotations[!no_sig_GO_found_mask]

# add new term column to each GO enrichment df
enriched_GO_annotations <- lapply(enriched_GO_annotations, function(ann) {
  ann@result$term <- paste(ann@result$ID, 
                    ann@result$Description, sep = ': ')
  return(ann)
})

# get enriched KEGG annotations
enriched_KEGG_annotations <- lapply(unique(all_sig_genes_infection$cluster), function(cluster){
  
  # isolate human genes in the cluster
  cluster_mask <- all_sig_genes_infection$cluster == cluster
  human_gene_mask <- rownames(all_sig_genes_infection) %in% all_sig_human_gene_names_LRT
  in_cluster <- all_sig_genes_infection[cluster_mask & human_gene_mask,]
  genes <- rownames(in_cluster)
  
  # get enriched KEGG pathways for the cluster genes
  enriched_KEGG_terms <- get_enriched_terms(genes, terms = 'KEGG')
  return(enriched_KEGG_terms)
  
})

# set names of the list
names(enriched_KEGG_annotations) <- unique(all_sig_genes_infection$cluster)

# mask for filtering results without any terms
no_sig_KEGG_found_mask <- sapply(enriched_KEGG_annotations, function(results){nrow(results@result) == 0})

# filter for clusters with enriched GO annotations found
enriched_KEGG_annotations <- enriched_KEGG_annotations[!no_sig_KEGG_found_mask]

# add new term column to each KEGG enrichment df
enriched_KEGG_annotations <- lapply(enriched_KEGG_annotations, function(ann) {
  ann@result$term <- paste(ann@result$ID, 
                           ann@result$Description, sep = ': ')
  return(ann)
})

# function for parsing enrichment results, looking up genes in the DE results
parse_ann <- function(ann, DESeq_results, n_terms = 5, k_genes = 3){
  
  if(nrow(ann) == 0){
    return(NA)
  } else if(nrow(ann) < n_terms){
    n_terms <- nrow(ann)
  } 
  
  top_n <- ann[1:n_terms,]
  
  terms <- paste(rownames(top_n), top_n$Description, sep = ':')
  
  pvalues <- top_n$pvalue
  
  gene_sets <- strsplit(top_n$geneID, '/') %>% 
    lapply(mapIds, column = 'SYMBOL', keytype = 'ENTREZID')
  
  top_k_genes <- sapply(gene_sets, function(genes){
    res <- DESeq_results[genes,]
    top_k <- rownames(res)[order(res$padj),][1:k_genes]
    top_k_cat <- paste(top_k, collapse = ',')
    return(top_k_cat)
  })
  
  parsed_ann <- data.frame(term = terms, 
                           p = pvalues,
                           genes = top_k_genes)
  
  return(parsed_ann)
  
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

GO_enrichment_matrix <- matrix(0, ncol = 4, nrow = length(top_GO_annotations$term))
rownames(GO_enrichment_matrix) <- top_GO_annotations$term
colnames(GO_enrichment_matrix) <- clusters_to_plot
for(term in unique(rownames(GO_enrichment_matrix))){
  
  term_mask <- which(rownames(GO_enrichment_matrix) == term)
  logps <- sapply(clusters_to_plot, function(cluster){
    ann <- enriched_GO_annotations[[cluster]]@result
    if(term %in% ann$term){
      if(ann[ann$term == term,'p.adjust'] < 0.05){
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
                                          filename = paste('figures/enriched_GO_terms_', cell_line, '.png', sep = ''), 
                                          width = 8, height = 5, cellwidth = 25, cellheight = 15)

# get top n GO annotations for each cluster
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

KEGG_enrichment_matrix <- matrix(0, ncol = 4, nrow = length(top_KEGG_annotations$term))
rownames(KEGG_enrichment_matrix) <- top_KEGG_annotations$term
colnames(KEGG_enrichment_matrix) <- clusters_to_plot
for(term in unique(rownames(KEGG_enrichment_matrix))){
  
  term_mask <- which(rownames(KEGG_enrichment_matrix) == term)
  logps <- sapply(clusters_to_plot, function(cluster){
    ann <- enriched_KEGG_annotations[[cluster]]@result
    if(term %in% ann$term){
      if(ann[ann$term == term,'p.adjust'] < 0.05){
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
                                            filename = paste('figures/enriched_KEGG_terms_', cell_line, '.png', sep = ''), 
                                            width = 7, height = 5, cellwidth = 25, cellheight = 15)
