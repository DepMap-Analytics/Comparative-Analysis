<h1>README</h1>
This directory contains scripts to reproduce Figures from Dempster et al, Agreement between two large pan-cancer CRISPR-Cas9 gene dependency datasets.


These scripts require the following R packages to be installed:

* CRISPRcleanR
* preprocessCore
* stringr
* sva
* tsne
* vcd
* reshape
* qvalue
* dendextend
* colorspace
* RColorBrewer
* matrixStats
* ggplot2
* ggrepel
* ggpubr
* beeswarm
* qusage

The scripts also require input data to be downloaded from Figshare: https://figshare.com/articles/Agreement_between_two_large_pan-cancer_CRISPR-Cas9_gene_dependency_datasets/7970993/1.
To run the scripts see the Initialise_script.R. This also requires the working directory to be set to the location of the Pipeline folder in this directory. The InputFolder directory needs to be set to the path where the downloaded data from Figshare is stored.
