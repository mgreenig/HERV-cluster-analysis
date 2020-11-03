library(ConsensusClusterPlus)

source('full_clustering.R')

# consensus clustering for HERVs
cc_results <- ConsensusClusterPlus::ConsensusClusterPlus(t(sig_retro_genes_infection_expr), maxK = 8, title = 'figures/HERV_consensus_cluster', 
                                                         reps = 1000, pItem = 0.8, seed = 42, plot = 'png')
# consensus clustering results for individual HERVs
icl <- ConsensusClusterPlus::calcICL(cc_results, title = 'figures/HERV_consensus_cluster', plot = 'png')