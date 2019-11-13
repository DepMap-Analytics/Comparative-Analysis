## read in the unscaled log FCs


broad_sgRNA<-read.csv(paste0(InputFolder,'LogfoldChange_Broad.csv'),header=TRUE,row.names=1,stringsAsFactors = FALSE)
sanger_sgRNA<-read.csv(paste0(InputFolder,'LogfoldChange_Sanger.csv'),header=TRUE,row.names=1,stringsAsFactors = FALSE)

#read in the guide maps:

broad_library<-read.csv(paste0(InputFolder,"GuideMap_Avana.csv"),header=T,stringsAsFactors = FALSE)
sanger_library<-read.csv(paste0(InputFolder,"GuideMap_KY10.csv"),header=T,stringsAsFactors = FALSE)

repmap_broad<-read.csv(paste0(InputFolder,"ReplicateMap_Broad.csv"),header=T,stringsAsFactors = F)
repmap_sanger<-read.csv(paste0(InputFolder,"ReplicateMap_Sanger.csv"),header=T,stringsAsFactors=F)

broad_sgRNA$gene<-broad_library[match(rownames(broad_sgRNA),as.character(broad_library$sgrna)),"gene"]
geneFC_broad<-aggregate(.~gene,broad_sgRNA,mean)

sanger_sgRNA$gene<-sanger_library[match(rownames(sanger_sgRNA),as.character(sanger_library$sgrna)),"gene"]
geneFC_sanger<-aggregate(.~gene,sanger_sgRNA,mean)

repmap_broad$mn<-make.names(repmap_broad$replicate_ID)
repmap_sanger$mn<-make.names(repmap_sanger$replicate_ID)

#get gene list for later:
genes_broad<-geneFC_broad$gene
genes_sanger<-geneFC_sanger$gene


#remove the gene list for second aggregation:
geneFC_sanger<-geneFC_sanger[,repmap_sanger$mn]
geneFC_broad<-geneFC_broad[,intersect(repmap_broad$mn,colnames(geneFC_broad))]

rownames(geneFC_sanger)<-genes_sanger
rownames(geneFC_broad)<-genes_broad

#now generate the cell line vector, rbind and aggregate over replicates as well
broad_cl<-repmap_broad[match(colnames(geneFC_broad),repmap_broad$mn),"cell_line"]
sanger_cl<-repmap_sanger[match(colnames(geneFC_sanger),repmap_sanger$mn),'cell_line']

dfFC_broad<-t(geneFC_broad)
dfFC_sanger<-t(geneFC_sanger)

dfFC_broad<-data.frame(dfFC_broad)
dfFC_sanger<-data.frame(dfFC_sanger)
dfFC_broad$cl<-broad_cl
dfFC_sanger$cl<-sanger_cl

geneFCcl_broad<-aggregate(.~cl,dfFC_broad,mean)
geneFCcl_sanger<-aggregate(.~cl,dfFC_sanger,mean)

rownames(geneFCcl_broad)<-geneFCcl_broad[,"cl"]
rownames(geneFCcl_sanger)<-geneFCcl_sanger[,"cl"]

write.csv(t(geneFCcl_broad[,2:ncol(geneFCcl_broad)]),file=paste0(ResultsFolder,"geneFCBroad.csv"),quote=F)

write.csv(t(geneFCcl_sanger[,2:ncol(geneFCcl_sanger)]),file=paste0(ResultsFolder,"geneFCSanger.csv"),quote=F)



broadExtraTP<-read.csv(paste0(InputFolder,'AlternateTimePointLFC_Broad.csv'),header=TRUE,row.names=1,stringsAsFactors = FALSE)
sangerExtraTP<-read.csv(paste0(InputFolder,'AlternateTimePointLFC_Sanger.csv'),header=TRUE,row.names=1,stringsAsFactors = FALSE)

broadExtraTP$gene<-broad_library[match(rownames(broadExtraTP),as.character(broad_library$sgrna)),"gene"]
geneFC_broadTP<-aggregate(.~gene,broadExtraTP,mean)

sangerExtraTP$gene<-sanger_library[match(rownames(sangerExtraTP),as.character(sanger_library$sgrna)),"gene"]
geneFC_sangerTP<-aggregate(.~gene,sangerExtraTP,mean)

HT29BroadD14<-rowMeans(geneFC_broadTP[,c("HT29.311Cas9_RepA_p4.AVANA","HT29.311Cas9_RepB_p4.AVANA")])
HT29BroadD21<-rowMeans(geneFC_broad[,c("HT29.311Cas9_RepA_p6.AVANA_batch3","HT29.311Cas9_RepB_p6.AVANA_batch3")])
names(HT29BroadD14)<-geneFC_broadTP$gene
names(HT29BroadD21)<-rownames(geneFC_broad)

HT29SangerD21<-rowMeans(geneFC_sangerTP[,c("HT29_c908R1_D21","HT29_c908R2_D21","HT29_c908R3_D21")])
HT29SangerD14<-rowMeans(geneFC_sangerTP[,c("HT29_c908R1_D14","HT29_c908R2_D14","HT29_c908R3_D14")])
names(HT29SangerD14)<-geneFC_sangerTP$gene
names(HT29SangerD21)<-geneFC_sangerTP$gene

HT29TP_Broad<-HT29BroadD21-HT29BroadD14[names(HT29BroadD21)]
HT29TP_Sanger<-HT29SangerD21-HT29SangerD14[names(HT29SangerD21)]

write.table(sort(HT29TP_Broad,decreasing=T),file=paste0(ResultsFolder,"BroadTP.rnk"),quote=F,sep="\t")
write.table(sort(HT29TP_Sanger,decreasing=T),file=paste0(ResultsFolder,"SangerTP.rnk"),quote=F,sep="\t")
#clear all variables after generating data
rm(list=ls())
