# Network analysis of endogenous retroviruses in SARS-Cov2 infection

This repository contains a pipeline that performs a re-analysis of an RNA-sequencing data set published by 
[Blanco-Melo et al (2020)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507). We generate a retroelement transcript count data set by applying 
[Telescope](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453) to the raw read data generated 
by Blanco-Melo and colleagues. A human/HERV gene correlation network is constructed and clusters of HERVs/human genes are produced by applying hierarchical clustering with a 
[dynamic tree cut](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/BranchCutting/).

## Dependencies

This package uses R version 3.6 and requires the following package versions:

- dplyr == 0.8.3
- ggplot2 == 3.2.1
- stringr == 1.4.0
- org.Hs.eg.db == 3.8.2
- DESeq2 == 1.24.0
- gridExtra == 2.3
- EnhancedVolcano == 1.2.0
- factoextra == 1.0.5
- dynamicTreeCut == 1.63-1
- pheatmap == 1.0.12
- clusterProfiler == 3.12.0
- RColorBrewer == 1.1-2
- reshape2 == 1.4.3
- ggpubr == 0.4.0

All packages are installed automatically by running the scripts (if not already present).

## Reproducing the analysis

First, ensure that R version 3.6 is installed and added to your system's environment path variable. 
To confirm that the correct version of R is being used, run from the command line:

```
which Rscript
```

To reproduce the clustering, run from the command line as follows:

```
Rscript cluster_trajectory_visualisation.R A549
```

for the A549 cell data, or:

```
Rscript cluster_trajectory_visualisation.R Calu3
```

for the Calu3 cell data. To reproduce the enrichment analysis, run:

```
Rscript cluster_enrichment_analysis.R A549
```

Replace A549 with Calu3 as before for the other data set.
