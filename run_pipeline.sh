#!/bin/bash

Rscript cluster_enrichment_analysis.R A549
Rscript cluster_trajectory_visualisation.R A549

Rscript cluster_enrichment_analysis.R Calu3
Rscript cluster_trajectory_visualisation.R Calu3
