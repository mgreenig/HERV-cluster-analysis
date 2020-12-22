# Network analysis of endogenous retroviruses in SARS-Cov2 infection

This repository contains a pipeline that performs a re-analysis of an RNA-sequencing data set published by 
[Blanco-Melo et al (2020)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507). We generate a retroelement transcript count data set by applying 
[Telescope](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453) to the raw read data generated 
by Blanco-Melo and colleagues. A human/HERV gene correlation network is constructed and clusters of HERVs/human genes are produced by applying hierarchical clustering with a 
[dynamic tree cut](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/BranchCutting/).





