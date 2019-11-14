
load(paste0(ExternalData,"OverlapInventory.RData"))
clinventory<-overlap_inventory

#unscaled logfc mean collapsed by gene and replicate

broad_data<-read.csv(paste0(ResultsFolder,"geneFCBroad.csv"),header=T,row.names=1,stringsAsFactors = F)
sanger_data<-read.csv(paste0(ResultsFolder,"geneFCSanger.csv"),header=T,row.names=1,stringsAsFactors = F)


## Attaching the study of origin identifier to each sample
colnames(broad_data)<-paste(colnames(broad_data),'---Broad',sep='')
colnames(sanger_data)<-paste(colnames(sanger_data),'---Sanger',sep='')


## individual sample-wise quantile normalisation of the two datasets
dn<-dimnames(broad_data)
d1<-normalize.quantiles(as.matrix(broad_data))
dimnames(d1)<-dn

dn<-dimnames(sanger_data)
d2<-normalize.quantiles(as.matrix(sanger_data))
dimnames(d2)<-dn

CombinedDataset<-cbind(d1,d2)
CombinedGDataset<-cbind(d1,d2)
## filtering out genes with a normLRT < 50 in both the datasets
broadtsk<-
  read.csv(paste0(InputFolder,'NormLRT_Broad.csv'),
           header = TRUE,stringsAsFactors = FALSE)
broadgenes<-intersect(broadtsk$gene[which(broadtsk$NormLRT>=100)],
                      rownames(CombinedDataset))
sangertsk<-
  read.csv(paste0(InputFolder,'NormLRT_Sanger.csv'),
           header = TRUE,stringsAsFactors = FALSE)
sangergenes<-intersect(sangertsk$gene[which(sangertsk$NormLRT>=100)],
                       rownames(CombinedDataset))
genes<-union(broadgenes,sangergenes)
CombinedDataset<-CombinedDataset[genes,]

site<-str_split(colnames(CombinedDataset),'---')
psite<-unlist(site)
site<-psite[seq(2,length(psite),2)]
cnames<-psite[seq(1,length(psite),2)]
## deriving COSMIC identifiers for each sample in the combined dataset
#
clinventory$MakeNames<-make.names(clinventory$Name)
cids<-clinventory$Cosmic.Identifier[match(cnames,clinventory$MakeNames)]

## loading Multiomic Cancer Functional Event binary matrix (MoBEM)
## and selecting columns corresponding to the cell lines in the inventory
load(paste0(ExternalData,'MoBEM.RData'))
includecids<-intersect(unique(cids),colnames(MoBEM))
MoBEM<-MoBEM[,includecids]
cname_inBEM<-cnames[!cnames=="MCAS"]
CombinedDataset<-CombinedDataset[,cnames%in%cname_inBEM]
CombinedGDataset<-CombinedGDataset[,colnames(CombinedDataset)]
## assembling study/cell-line of derivation for each sample in the combined dataset

site<-str_split(colnames(CombinedDataset),'---')
psite<-unlist(site)
site<-psite[seq(2,length(psite),2)]
cnames<-psite[seq(1,length(psite),2)]



## creating tissue dummy variables and attaching them to the MoBEM
tissues<-clinventory$Tissue[match(colnames(MoBEM),clinventory$Cosmic.Identifier)]
allt<-unique(tissues)
TissueBEM<-t(do.call(rbind,lapply(tissues,function(x) (allt==x)+0)))
rownames(TissueBEM)<-allt
MoBEM<-rbind(MoBEM,TissueBEM)

## creating MSS_status dummy variables and attaching them to the MoBEM
MSS<-clinventory$MSS_status[match(colnames(MoBEM),clinventory$Cosmic.Identifier)]
MSS[MSS=='MSI-H']<-1
MSS[MSS!='1']<-0
MSS<-as.numeric(MSS)
MSS[is.na(MSS)]<-0
MoBEM<-rbind(MoBEM,MSS)
rownames(MoBEM)[nrow(MoBEM)]<-'MSI'

## decode CNA identifiers in the MoBEM

load(paste0(ExternalData,"CNAdecode.RData"))
decodeCNA_cp<-function(MoBEM){

  rn <- rownames(MoBEM)
  ii <- grep("cna", rownames(MoBEM))
  cnaId <- rownames(MoBEM)[ii]
  containedGenes <- unlist(lapply(str_split(cnaId, " "), function(x) {
    x[2]
  }))
  containedGenes[is.na(containedGenes)] <- ""
  segments <- unlist(lapply(str_split(unlist(lapply(str_split(cnaId, 
                                                              ":"), function(x) {
                                                                x[2]
                                                              })), " "), function(x) {
                                                                x[1]
                                                              }))
  loci <- as.character(CNAdecode$locus[match(segments, CNAdecode$Identifier)])
  altType <- as.character(CNAdecode$Recurrent.Amplification..Deletion[match(segments, 
                                                                            CNAdecode$Identifier)])
  altType[altType == "Amplification"] <- "G:"
  altType[altType == "Deletion"] <- "L:"
  rownames(MoBEM)[ii] <- paste(altType, loci, " ", containedGenes, 
                               sep = "")
  return(MoBEM)
}

MoBEM<-decodeCNA_cp(MoBEM)



## Performing systematic interaction test

## Filtering out cancer functional events present in less than 2 and more than 144 cell lines
MoBEM<-MoBEM[rowSums(MoBEM)>2,]
MoBEM<-MoBEM[rowSums(MoBEM)<144,]


RES_S<-matrix(NA,(nrow(CombinedDataset))*nrow(MoBEM),6,
              dimnames = list(NULL,c('CFE','GENE','delta','effect_size','p','fdr')))
RES_B<-matrix(NA,(nrow(CombinedDataset))*nrow(MoBEM),6,
              dimnames = list(NULL,c('CFE','GENE','delta','effect_size','p','fdr')))
RES_S<-as.data.frame(RES_S)
RES_B<-as.data.frame(RES_B)
nfet<-nrow(MoBEM)
count<-1
## deriving COSMIC identifiers for each sample in the combined dataset
#
clinventory$MakeNames<-make.names(clinventory$Name)
cids<-clinventory$Cosmic.Identifier[match(cnames,clinventory$MakeNames)]
for (i in 1:nfet){

  print(paste('Cancer functional event tested n.', i,sep=''))
  CFE<-rownames(MoBEM)[i]

  Sa<-diffDep(MoBEM = MoBEM,site = site,dataset = as.matrix(CombinedDataset),
              CFE = CFE,s = 'Sanger',display = FALSE)
  Br<-diffDep(MoBEM = MoBEM,site = site,dataset = as.matrix(CombinedDataset),
              CFE = CFE,s = 'Broad',display = FALSE)

  RES_S[count:(count+nrow(Sa)-1),]<-as.data.frame(Sa)
  RES_B[count:(count+nrow(Br)-1),]<-as.data.frame(Br)

  count<-count+nrow(Sa)
}


## Correct pvalues overall tests
fdr_oa<-p.adjust(RES_S$p,'fdr')
RES_S<-cbind(RES_S,fdr_oa)

fdr_oa<-p.adjust(RES_B$p,'fdr')
RES_B<-cbind(RES_B,fdr_oa)

## saving results
save(RES_S,file=paste0(ResultsFolder,'RES_S_red.RData'))
save(RES_B,file=paste0(ResultsFolder,'RES_B_red.RData'))

## interaction comparison plot, to interactively identify points, set IDENTIFY
## parameter to TRUE


iSanger_s <- which(RES_S$fdr_oa < 0.05 & RES_S$delta < -1)
iBroad_s <- which(RES_B$fdr_oa < 0.05 & RES_B$delta < -1)
iboth_s <- intersect(iSanger_s, iBroad_s)
iEither <- union(iSanger_s, iBroad_s)

selectedPoints<-read.csv(paste0(ExternalData,file="SelectedPoints_BioMarker.csv"),stringsAsFactors = F)
selectedPoints<-selectedPoints[!is.na(selectedPoints[,"X"]),]

plot_data<-RES_S
plot_data$Broad<-RES_B$delta
plot_data$Sanger<-RES_S$delta
both_data<-plot_data[iboth_s,]
label_data<-plot_data[as.character(selectedPoints[,"X"]),]
label_data$Label<-selectedPoints[match(rownames(label_data),selectedPoints[,"X"]),"Label"]
###Figure 3 a #### 
pdf(paste0(ResultsFolder,"BioMarker.pdf"),width=3.5,height=3.5)
par(mar=c(1,1,1.2,1)+0.1,mgp=c(0.1,0.1,0))
ggplot(plot_data, aes(x=Broad,y=Sanger))+geom_point(alpha=0.5,col='grey')+
  xlim(-2.5,1)+
  ylim(-2.5,1)+
  geom_hline(yintercept= -1,linetype="dashed",col="gray")+
  geom_vline(xintercept = -1,linetype="dashed",col="gray")+
  geom_hline(yintercept= 0,linetype="dashed",col="gray")+
  geom_vline(xintercept = 0,linetype="dashed",col="gray")+
  theme_bw()+
  geom_text_repel(data=label_data,aes(x=Broad,y=Sanger,label=paste(label_data$Label,label_data$GENE,sep="\n")),size=2,inherit.aes = F,box.padding=0.5,point.padding = 0.5)+
  stat_density_2d(n=50,geom="polygon",aes(fill=..level..),show.legend = FALSE)+scale_fill_gradient(low = "grey", high = "white")+
  geom_point(data=both_data,aes(x=Broad,y=Sanger),col='purple')+
  labs(x="Differential essentiality [Broad]",y="Differential essentiality [Sanger]",cex=0.8)

dev.off()


SAvsBR<-PrRc_and_ROC_Curves(loc_RES_S = RES_S,loc_RES_B = RES_B,es_th = 0,
                            fdr_th = 0.05,study1name = 'Sanger',study2name = 'Broad')
BRvsSA<-PrRc_and_ROC_Curves(loc_RES_S = RES_B,loc_RES_B = RES_S,es_th = 0,
                            fdr_th = 0.05,study1name = 'Broad',study2name = 'Sanger')


###Figure 3 b### 
##SA is Broad association vs Sanger ranked p-value
par(mar=c(0.1,0.1,0.1,0.1)+0.1,mgp=c(0.1,0.1,0))
pdf(paste0(ResultsFolder,"SA_error1.pdf"),height=4.5,width=4.5)
SA_pr<-PR(loc_RES_S = RES_S,loc_RES_B = RES_B,es_th = 0,
                            fdr_th = 0.05,study1name = 'Sanger',study2name = 'Broad')
dev.off()
pdf(paste0(ResultsFolder,"SA_error2.pdf"),height=4.5,width=4.5)
SA_roc<-ROC(loc_RES_S = RES_S,loc_RES_B = RES_B,es_th = 0,
          fdr_th = 0.05,study1name = 'Sanger',study2name = 'Broad')
dev.off()
##Br is Sanger association vs Broad ranked p-value
pdf(paste0(ResultsFolder,"Br_error1.pdf"),height=4.5,width=4.5)
BR_pr<-PR(loc_RES_S = RES_B,loc_RES_B = RES_S,es_th = 0,
          fdr_th = 0.05,study1name = 'Broad',study2name = 'Sanger')
dev.off()
pdf(paste0(ResultsFolder,"Br_error2.pdf"),height=4.5,width=4.5)
BR_roc<-ROC(loc_RES_S = RES_B,loc_RES_B = RES_S,es_th = 0,
            fdr_th = 0.05,study1name = 'Broad',study2name = 'Sanger')
dev.off()


#get associations passing significance threshold for each data set:
RES_B_signif<-RES_B[RES_B$fdr_oa<0.05&RES_B$delta< (-1),]
RES_S_signif<-RES_S[RES_S$fdr_oa<0.05&RES_S$delta< (-1),]

inboth<-RES_B_signif[intersect(rownames(RES_S_signif),rownames(RES_B_signif)),]


## Plot examples of consistent associations


###Figure 3c ###
pdf(paste0(ResultsFolder,"ExampleBiomarkers.pdf"),width=6.5,height=2.25)
par(mfrow=c(1,8))
par(mar=c(0.3,2,0.5,0.1)+0.1,mgp=c(1.25,0.5,0))
plotAssociationExamplesCP(gene = 'MYCN',feat = 'Peripheral Nervous System',
                        loc_CombinedDataset = CombinedGDataset,loc_site = site,loc_cids = cids)
plotAssociationExamplesCP('ERBB2','G:17q12 (CDK12,ERBB2,MED24)',
                        CombinedGDataset,site,cids)
plotAssociationExamplesCP('CTNNB1','APC_mut',
                        CombinedGDataset,site,cids)

plotAssociationExamplesCP('CTNNB1','chr1:120835961-120839391(FAM72B)_HypMET',
                        CombinedGDataset,site,cids)
dev.off()



