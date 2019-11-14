setwd("/path/to/directory/containing/Pipelines/")
ExternalData<-"../ExternalData/"
ResultsFolder<-"../ResultsFolder/"
InputFolder<-"/path/to/data/downloaded/from/figshare"


library(CRISPRcleanR)
library(preprocessCore)
library(stringr)
library(sva)
library(tsne)
library(vcd)
library(reshape)
library(qvalue)
library(dendextend)
library(colorspace)
library(RColorBrewer)
library(matrixStats)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(beeswarm)
library(qusage)

data("BAGEL_essential")
data("BAGEL_nonEssential")

#source functions for use in analysis:
source("../Functions/CRISPRcollatio.R")

#script 1 to generate batch corrected data. Figures 2 d,e,f and Supplementary Figures 3a,b,c
source("BatchCorrection_and_Checks.R")

#script 2 to compare molecular biomarkers across the two datasets. Generates Figures 3 a,b,c
source("DepMarkers_agreement.R")

#script 3 for gene expression biomarkers. Generates Figures 3 d,e and Supplementary Figure 4
source("GeneExpressionBiomarkers.R")

#script 4 for Gene ontology analysis of Broad or Sanger exclusive dependencies. Figure 5c. 
source("DisagreementGOenrich.R")

#script 5 for GSEA analysis of Principal component loadings. Generates Supplementary Figure 5
source("GSEAPCloadings.R")