library(clusterProfiler)

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
    enriched <- enrichGO(entrez_ids, 'org.Hs.eg.db', ont = 'ALL')
  }else{
    enriched <- enrichKEGG(entrez_ids, organism = 'hsa')
  }
  
  return(enriched)
  
}

# get enriched GO annotations for each cluster
enriched_GO_annotations <- lapply(unique(all_sig_genes_infection$cluster), function(cluster){
  
  # isolate genes in the cluster
  cluster_mask <- all_sig_genes_infection$cluster == cluster
  in_cluster <- all_sig_genes_infection[cluster_mask,]
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

# get enriched KEGG annotations
enriched_KEGG_annotations <- lapply(unique(all_sig_genes_infection$cluster), function(cluster){
  
  # isolate genes in the cluster
  cluster_mask <- all_sig_genes_infection$cluster == cluster
  in_cluster <- all_sig_genes_infection[cluster_mask,]
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