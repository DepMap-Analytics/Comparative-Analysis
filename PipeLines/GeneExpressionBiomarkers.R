
CombineRNA<-read.csv(file=paste0(ExternalData,"CombinedRNAseqData.csv"),stringsAsFactors = F,row.names = 1)
OverlapData<-read.table(paste0(ExternalData,"OverlapCLdata.tsv"),header=T,stringsAsFactors = FALSE)
cellmodelInfo<-read.csv(paste0(ExternalData,"model_list_cellmodelpassports.csv"),header=T,stringsAsFactors = FALSE,row.names=1)
#get COSMIC ids for cell lines in overlap:

OverlapData$COSMIC<-cellmodelInfo[rownames(OverlapData),"COSMIC_ID"]
clnames<-make.names(OverlapData[,"model_name"])
OverlapData$mn<-clnames

colnames(CombineRNA)<-sapply(colnames(CombineRNA),function(x) substr(x,2,nchar(x)))

dd<-dist(scale(t(CombineRNA)),method="euclidean")
hc<-hclust(dd,method="ward.D2")
dend<-as.dendrogram(hc)

#set labels to cell line names
sitelabels<-unlist(sapply(dend%>%labels,function(x) strsplit(x,"...",fixed=T)[[1]][2]))
idlabels<-unlist(sapply(dend%>%labels,function(x) strsplit(x,"...",fixed=T)[[1]][1]))
idlabels2<-OverlapData[match(idlabels,OverlapData$COSMIC),"model_name"]
idlabels2[seq(1,length(idlabels2),2)]<-""

#set colour according to tissue lineage
colgrp<-cellmodelInfo[match(idlabels,cellmodelInfo$COSMIC_ID),"tissue"]
ntissues<-length(unique(colgrp))

cols_4 <- rainbow_hcl(ntissues)
names(cols_4)<-unique(colgrp)
colvec <- cols_4[colgrp]
par(mar=c(3.1,0.1,2.1,0.1))
dendro<-dend%>%set("labels",idlabels2)%>%set("labels_cex",0.55)%>% set("leaves_pch", c(18, 18))%>%set("leaves_cex", 0.5)%>%set("leaves_col",c("red","blue"))
labels_colors(dendro)<-colvec
##Supplementary Figure 4 
pdf(paste0(ResultsFolder,"RNA_dendro.pdf"),width=5.5,height=10)
plot(dendro,axes=F,horiz=T,cex.lab=0.8,cex.axis=0.8)
dev.off()

#load CRISPR data
broad_data<-read.csv(paste0(ExternalData,"geneFCBroad.csv"),header=T,row.names=1,stringsAsFactors = F)
sanger_data<-read.csv(paste0(ExternalData,"geneFCSanger.csv"),header=T,row.names=1,stringsAsFactors = F)

dn<-dimnames(broad_data)
broad_data<-normalize.quantiles(as.matrix(broad_data))
dimnames(broad_data)<-dn

dn<-dimnames(sanger_data)
sanger_data<-normalize.quantiles(as.matrix(sanger_data))
dimnames(sanger_data)<-dn

#load strongly selective dependencies to test for gene expression biomarkers
SSDs_genes2<-read.csv(paste0(ExternalData,"StronglySelectiveDependencies.csv"),header=T,stringsAsFactors = F)


broad_gene<-broad_data[intersect(rownames(broad_data),SSDs_genes2[,1]),Code$mn]
sanger_gene<-sanger_data[intersect(rownames(sanger_data),SSDs_genes2[,1]),Code$mn]
colnames(broad_gene)<-Code$COSMIC
colnames(sanger_gene)<-Code$COSMIC

BiomarkerRNA<-RNAbiomarkers(CombineRNA,broad_gene,sanger_gene)
RNA_data<-BiomarkerRNA$PlotRNAdata
minT_B<-BiomarkerRNA$B_Sthresh[1]
minT_S<-BiomarkerRNA$B_Sthresh[2]

#try plotting a subset of data and with the ggplot label information and lines so can see what is what
RNA_BroadSignif<-RNA_data[abs(RNA_data$Broad)>=minT_B,]
RNA_SangerSignif<-RNA_data[abs(RNA_data$Sanger)>=minT_S,]

RNA_signifBoth<-RNA_data[abs(RNA_data$Broad)>=minT_B&abs(RNA_data$Sanger)>=minT_S,]
RNA_signifBoth<-RNA_signifBoth[order(abs(RNA_signifBoth$Broad),decreasing=T),]

selectedPoints<-read.csv(paste0(ExternalData,file="SelectedPoints_RNA_ssd.csv"),stringsAsFactors = F, skipNul = TRUE)
selectedPoints<-selectedPoints[!is.na(selectedPoints[,"X"]),]

##Figure 3 d  ###
pdf(paste0(ResultsFolder,"ExpressionBioMarker.pdf"),width=4,height=4)
par(mar=c(1,1,1.2,1)+0.1,mgp=c(0.1,0.1,0))
label_data<-RNA_signifBoth[as.character(selectedPoints[,"X"]),]
ggplot(RNA_data, aes(x=Broad,y=Sanger))+geom_point(alpha=0.5,col='purple')+
  xlim(-1,1)+
  ylim(-1,1)+
  geom_hline(yintercept=minT_B,linetype="dashed",col="grey")+
  geom_vline(xintercept = minT_S,linetype="dashed",col="grey")+
  geom_hline(yintercept= -1*minT_B,linetype="dashed",col="grey")+
  geom_vline(xintercept = -1*minT_S,linetype="dashed",col="grey")+
  theme_bw()+
  stat_density_2d(n=35,geom="polygon",aes(fill=..level..),show.legend = FALSE)+scale_fill_gradient(low = "purple", high = "white")+
  geom_text_repel(data=label_data,aes(x=Broad,y=Sanger,label=paste(label_data$gene,label_data$SSD,sep="\n")),size=3,inherit.aes = F,box.padding=0.25,point.padding = 0.25)+
  labs(x="Gene expression & fitness correlation, Broad",y="Gene expression & fitness correlation, Sanger",cex=0.7)

dev.off()

###Figure 3 e  ###

par(mar=c(2.5,2.5,2,1.5)+0.1,mgp=c(1.2,0.4,0.2))

Broad_RNA<-CombineRNA[,grep("Broad",colnames(CombineRNA))]
Sanger_RNA<-CombineRNA[,grep("Sanger",colnames(CombineRNA))]
colnames(Broad_RNA)<-sapply(colnames(Broad_RNA),function(x) strsplit(x,"...",fixed=TRUE)[[1]][1])
colnames(Sanger_RNA)<-sapply(colnames(Sanger_RNA),function(x) strsplit(x,"...",fixed=TRUE)[[1]][1])
d1<-Broad_RNA["ATP6V0E2",]
d1<-melt(d1)
colnames(d1)<-c("COSMIC","RNA")
d1$FC<-broad_gene["ATP6V0E1",colnames(Broad_RNA)]
pdf(paste0(ResultsFolder,"correlationATP_Broad.pdf"),width=3,height=3)
sp <- ggscatter(d1, x = "RNA", y = "FC",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0.5, label.y = 0.6)+labs(x="Gene expression",y="log FC",title="")

dev.off()
d1<-Sanger_RNA["ATP6V0E2",]
d1<-melt(d1)
colnames(d1)<-c("COSMIC","RNA")
d1$FC<-sanger_gene["ATP6V0E1",colnames(Sanger_RNA)]
pdf(paste0(ResultsFolder,"correlationATP_Sanger.pdf"),width=3,height=3)
sp <- ggscatter(d1, x = "RNA", y = "FC",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0.5, label.y = 0.5)+labs(x="Gene expression",y="log FC",title="")

dev.off()
d1<-Broad_RNA["ERBB2",]
d1<-melt(d1)
colnames(d1)<-c("COSMIC","RNA")
d1$FC<-broad_gene["ERBB2",colnames(Broad_RNA)]
pdf(paste0(ResultsFolder,"correlationERBB2_Broad.pdf"),width=3,height=3)
sp <- ggscatter(d1, x = "RNA", y = "FC",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 3, label.y = 0.5)+labs(x="",y="log FC",title="Broad")+scale_x_continuous(breaks=c(2,6,10))

dev.off()
d1<-Sanger_RNA["ERBB2",]
d1<-melt(d1)
colnames(d1)<-c("COSMIC","RNA")
d1$FC<-sanger_gene["ERBB2",colnames(Sanger_RNA)]
pdf(paste0(ResultsFolder,"correlationERBB2_Sanger.pdf"),width=3,height=3)
sp <- ggscatter(d1, x = "RNA", y = "FC",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 3, label.y = 0.5)+labs(x="",y="log FC",title="Sanger")

dev.off()

