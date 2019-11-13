

## Generating gene logFCs 
source("Generate_unscaledLogFC.R")
broad_data<-read.csv(paste0(ResultsFolder,"geneFCBroad.csv"),header=T,row.names=1,stringsAsFactors = F)
sanger_data<-read.csv(paste0(ResultsFolder,"geneFCSanger.csv"),header=T,row.names=1,stringsAsFactors = F)
broadqn<-normalize.quantiles(as.matrix(broad_data))
sangerqn<-normalize.quantiles(as.matrix(sanger_data))
dimnames(broadqn)<-dimnames(broad_data)
dimnames(sangerqn)<-dimnames(sanger_data)
broad_scaled<-t(scale(t(broad_data)))
sanger_scaled<-t(scale(t(sanger_data)))
broad_scaledq<-t(scale(t(broadqn)))
sanger_scaledq<-t(scale(t(sangerqn)))
## Attaching the study of origin identifier to each sample
colnames(broad_scaled)<-paste(colnames(broad_scaled),'---Broad',sep='')
colnames(sanger_scaled)<-paste(colnames(sanger_scaled),'---Sanger',sep='')
colnames(broad_scaledq)<-paste(colnames(broad_scaledq),'---Broad',sep='')
colnames(sanger_scaledq)<-paste(colnames(sanger_scaledq),'---Sanger',sep='')
corrected_scaled<-cbind(broad_scaled,sanger_scaled[rownames(broad_scaled),])
REScorrected_scaled<-classPerf(corrected_scaled)
dimnames(corrected_scaled)<-dimnames(corrected_scaled)
corrected_scaledq<-cbind(broad_scaledq,sanger_scaledq[rownames(broad_scaledq),])
corrected_scaledq<-normalize.quantiles(corrected_scaledq)
dimnames(corrected_scaledq)<-dimnames(corrected_scaled)
REScorrected_scaledq<-classPerf(corrected_scaledq)
## Attaching the study of origin identifier to each sample
colnames(broad_data)<-paste(colnames(broad_data),'---Broad',sep='')
colnames(sanger_data)<-paste(colnames(sanger_data),'---Sanger',sep='')

## Collating the two dataset on the domain of the genes screened in both studies
mergedDataset<-
  cbind(broad_data,sanger_data[rownames(broad_data),])

## Quantile-normalising the combined dataset
dn<-dimnames(mergedDataset)
mergedDataset<-normalize.quantiles(as.matrix(mergedDataset))
dimnames(mergedDataset)<-dn

## Storing the site of origin of each sample in a vector
site<-str_split(colnames(mergedDataset),'---')
site<-unlist(site)
site<-site[seq(2,length(site),2)]
names(site)<-colnames(mergedDataset)

## Batch correcting the combined dataset using the site of origin of the
## samples as batch-covariate
corrected<-ComBat(as.matrix(mergedDataset),batch = site[colnames(mergedDataset)])

## Quantile-normalising the corrected dataset
dn<-dimnames(corrected)
corrected<-normalize.quantiles(as.matrix(corrected))
dimnames(corrected)<-dn

REScorrected<-classPerf(corrected)

SSDs_genes<-read.csv(paste0(InputFolder,"StronglySelectiveDependencies.csv"),header=T,stringsAsFactors = F)

## Storing identifiers of screened cell lines in a vector
CLs<-unlist(lapply(str_split(colnames(mergedDataset),'---'),function(x) x[1]))
CLs<-CLs[1:147]
GeneSetOnesignif<-GetSigProfiles(mergedDataset,CLs,BAGEL_essential,BAGEL_nonEssential)
OneSignif<-union(GeneSetOnesignif[[1]],GeneSetOnesignif[[2]])
write.csv(OneSignif,file=paste0(ExternalData,"GenesVariable.csv"))


dataset<-mergedDataset[OneSignif,]
clnames<-str_split(colnames(dataset),'---')
clnames<-unlist(lapply(clnames,function(x){x[1]}))
uc<-unique(clnames)
ncl<-length(uc)
site <- str_split(colnames(dataset), "---")
site <- unlist(lapply(site, function(x) {
  x[[2]]
}))
colnames(dataset) <- clnames
cdist <- as.dist(1 - cor(dataset))
set.seed(679661)
tfitM <- tsne(cdist, max_iter = 1000, perplexity = 100)
dataset<-corrected[OneSignif,]
cdist <- as.dist(1 - cor(dataset))
set.seed(679661)
tfitC <- tsne(cdist, max_iter = 1000, perplexity = 100)

color = grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), 
                                 invert = T)]
uc <- unique(clnames)
COLS <- sample(color, length(uc))
names(COLS) <- uc
tCOLS <- makeTransparent(COLS, alpha = 140)
SYMBOLS <- rep(21, length(site))
SYMBOLS[site == "Broad"] <- 23
names(SYMBOLS) <- site
CEX <- rep(1.2, length(site))
CEX[site == "Broad"] <- 1.4
#Figures 2d and e ### MAIN FIGURE NEEDS UPDATING TO THIS ONE ####
pdf(paste0(ResultsFolder,"tSNE_OneSignif.pdf"),width=6,height=3.25)
par(mfrow=c(1,2))
par(mar=c(1,1,1.2,1)+0.1)
plot(tfitM[, 1], tfitM[, 2], bg = tCOLS[clnames], col = "black", 
     cex = CEX, pch = SYMBOLS[site], frame.plot = FALSE, xaxt = "n", 
     yaxt = "n", xlab = "tSNE AU", ylab = "tSNE AU", main = "Uncorrected fitness profiles",cex.main=0.9)
plot(tfitC[, 1], tfitC[, 2], bg = tCOLS[clnames], col = "black", 
     cex = CEX, pch = SYMBOLS[site], frame.plot = FALSE, xaxt = "n", 
     yaxt = "n", xlab = "tSNE AU", ylab = "tSNE AU", main = "Batch corrected fitness profiles",cex.main=0.9)
dev.off()

## knn classifier performances' plots
RESmerged<-classPerf(mergedDataset[OneSignif,])
REScorrected<-classPerf(corrected)
REScorrected_SSDs<-classPerf(corrected[SSDs_genes[,1],])
REScorrected_SSDs2<-classPerf(corrected[OneSignif,])


###Figure 2f #### 
pdf(paste0(ResultsFolder,"batchcorrected_ROCcurves.pdf"),width=3,height=3)

par(mfrow=c(1,1))
par(mar=c(2.5,2.5,0.1,0.2)+0.1,mgp=c(1,0.25,0),xpd=TRUE)
plot(1:293,100*RESmerged$CURV,frame.plot = FALSE,col=makeTransparent('brown',150),lwd=4,type='l',
     xlab='k-neighbourhood',ylab='% cell lines matching with counterpart',cex.axis=0.7,cex.lab=0.7,tcl=0.5,tck=-0.01)
lines(1:293,100*1:293*1/293,col=makeTransparent('black',150))
lines(1:293,100*REScorrected$CURV,col=makeTransparent('red',150),lwd=4,type='l')
lines(1:293,100*REScorrected_SSDs2$CURV,col=makeTransparent('purple',150),lwd=4,type='l')
lines(1:293,100*REScorrected_SSDs$CURV,col=makeTransparent('orange',150),lwd=4,type='l')
nAUCmerged=trapz(1:293,100*RESmerged$CURV)/(100*293)
nAUCcorrected=trapz(1:293,100*REScorrected$CURV)/(100*293)
nAUCcorrected_SSD2=trapz(1:293,100*REScorrected_SSDs2$CURV)/(100*293)
nAUCcorrected_SSD=trapz(1:293,100*REScorrected_SSDs$CURV)/(100*293)

legend('bottomright',legend=paste(c('Uncorrected','Corrected all','Corrected variable',
                                    'Corrected SSD','Random'),
                                  ' (',
                                  format(c(nAUCmerged,
                                           nAUCcorrected,
                                           nAUCcorrected_SSD2,
                                           nAUCcorrected_SSD,0.5),digits=2),
                                  ')',sep=''),
       col=c(makeTransparent('brown',150),
             makeTransparent('red',150),
             makeTransparent('purple',150),
             makeTransparent('orange',150),makeTransparent('black',150)),lwd=2,cex=0.6,bty="n")



dev.off()


###Supplementary Figure 2a #### 
pdf(paste0(ResultsFolder,"Unprocessed_batchcorrected_tsneAllgenes.pdf"),width=7,height=3.5)

par(mfrow=c(1,2))
par(mar=c(1,1,1.2,1)+0.1)
tSNEplot(mergedDataset,title = 'Uncorrected dataset',perplexity = 100)
tSNEplot(corrected,title='Batch corrected combined dataset',perplexity = 100)
dev.off()

#### Supplementary Figure 2b #### 
pdf(paste0(ResultsFolder,"Unprocessed_batchcorrected_tsneSSDgenes.pdf"),width=7,height=3.5)
par(mfrow=c(1,2))
par(mar=c(1,1,1.2,1)+0.1)
tSNEplot(mergedDataset[SSDs_genes[,1],],title='Uncorrected SSD genes',perplexity = 100)

tSNEplot(corrected[SSDs_genes[,1],],title='Batch corrected SSD genes',perplexity = 100)
dev.off()


## Comparing n.significantly essential genes before/after batch correction
## across the two datasets
Ndepletions<-vector()
Recalls<-vector()

Ndepletions_corrected<-vector()
Recalls_corrected<-vector()

for (i in 1:ncol(mergedDataset)){
  FCs<-mergedDataset[,i]
  names(FCs)<-rownames(mergedDataset)
  RES<-nDepletions(FCs)
  Ndepletions[i]<-RES$Ndep
  Recalls[i]<-RES$Recall
  FCs<-corrected[,i]
  names(FCs)<-rownames(corrected)
  RES<-nDepletions(FCs)
  Ndepletions_corrected[i]<-RES$Ndep
  Recalls_corrected[i]<-RES$Recall
}

names(Ndepletions)<-colnames(mergedDataset)
names(Recalls)<-colnames(Recalls)
names(Ndepletions_corrected)<-colnames(corrected)
names(Recalls_corrected)<-colnames(Recalls)

SangerId<-grep('Sanger',names(Ndepletions))
BroadId<-grep('Broad',names(Ndepletions))
##Supplementary Figure 3a #### 
pdf(paste0(ResultsFolder,"Unprocessed_batchcorrected_NumberDepletions.pdf"),width=4.5,height=2.5)
par(mfrow=c(1,2))
par(mar=c(2.5,2.5,0.5,0.2)+0.1,mgp=c(1.2,0.25,0))
plot(Ndepletions[SangerId],Ndepletions[BroadId],
     xlab='Sanger n. dependencies',
     ylab='Broad n. dependencies',
     xlim=c(200,3000),ylim=c(200,3000),bg=makeTransparent('darkgreen',100),col='darkgreen',pch=21,
     main='Merged dataset',frame.plot=F,cex.main=0.8,cex.axis=0.7,cex.lab=0.8)
abline(0,1,col='gray')
abline(v=median(Ndepletions[SangerId]),col='blue',lwd=2,lty=2)
abline(h=median(Ndepletions[BroadId]),col='purple',lwd=2,lty=2)
plot(Ndepletions_corrected[SangerId],Ndepletions_corrected[BroadId],
     xlab='Sanger n. dependencies',ylab='',
     xlim=c(200,3000),ylim=c(200,3000),bg=makeTransparent('darkgreen',100),col='darkgreen',pch=21,
     main='Batch corrected',frame.plot=F,cex.main=0.8,cex.axis=0.7,cex.lab=0.8)
abline(0,1,col='gray')
abline(v=median(Ndepletions_corrected[SangerId]),col='blue',lwd=2,lty=2)
abline(h=median(Ndepletions_corrected[BroadId]),col='purple',lwd=2,lty=2)
dev.off()


paste('Median number depletions pre correction Broad:',median(Ndepletions[BroadId]))
paste('Median number depletions pre correction Sanger:',median(Ndepletions[SangerId]))
paste('Median number depletions post correction Broad:',median(Ndepletions_corrected[BroadId]))
paste('Median number depletions post correction Sanger:',median(Ndepletions_corrected[SangerId]))



## Supplementary Figure 3b #### 
pdf(paste0(ResultsFolder,"Unprocessed_batchcorrected_DistancePlots.pdf"),width=11,height=3)

par(mfrow=c(2,3))
par(mar=c(2.5,2.5,1.5,0.2)+0.1,mgp=c(1,0.25,0))
distPlot(mergedDataset,title='All genes',XLIM = c(0.2,1),YLIMS = c(0,8.5))
distPlot(mergedDataset[OneSignif,],title='Variable Genes',XLIM = c(-0.3,1),YLIMS = c(0,9))
distPlot(mergedDataset[SSDs_genes[,1],],title='SSD Genes',XLIM = c(-0.3,1),YLIMS = c(0,5.5))

distPlot(corrected,title='All genes (corrected)',XLIM = c(0.2,1),YLIMS = c(0,8.5))
distPlot(corrected[OneSignif,],title='Variable Genes (corrected)',XLIM = c(-0.3,1),YLIMS = c(0,9))
distPlot(corrected[SSDs_genes[,1],],title='SSD Genes (corrected)',XLIM = c(-0.3,1),YLIMS = c(0,6))
dev.off()



#Dependency agreement:

## Loading reference sets of essential/non-essential genes
essential_genes<-read.csv(paste0(ExternalData,"EssentialGenes.csv"),header=T,stringsAsFactors = FALSE)
Nonessential_genes<-read.csv(paste0(ExternalData,"NonessentialGenes.csv"),header=T,stringsAsFactors = FALSE)



## identifying and plotting significant dependencies (at 5% FDR) before/after batch corretion, across cell lines
clnames<-str_split(colnames(mergedDataset),'---')
clnames<-unlist(lapply(clnames,function(x){x[1]}))
uc<-unique(clnames)
ncl<-length(uc)

cols<-rep(makeTransparent('gray',150),nrow(mergedDataset))
names(cols)<-rownames(mergedDataset)

BroadOnly<-list()
SangerOnly<-list()
common<-list()

BroadOnly_corrected<-list()
SangerOnly_corrected<-list()
common_corrected<-list()

for (i in 1:ncl){
  #print(c(clnames[i],i,ncl))
  
  id<-which(clnames==uc[i])
  
  FCs<-mergedDataset[,id[1]]
  names(FCs)<-rownames(mergedDataset)
  RES<-PrRc_Curve(FCs,essential_genes[,1],Nonessential_genes[,1],display = FALSE,FDRth = 0.05)
  depBroad<-names(which(FCs<=RES$sigthreshold))
  
  FCs<-mergedDataset[,id[2]]
  names(FCs)<-rownames(mergedDataset)
  RES<-PrRc_Curve(FCs,essential_genes[,1],Nonessential_genes[,1],display = FALSE,FDRth = 0.05)
  depSanger<-names(which(FCs<=RES$sigthreshold))
  
  BroadOnly[[i]]<-setdiff(depBroad,depSanger)
  SangerOnly[[i]]<-setdiff(depSanger,depBroad)
  common[[i]]<-intersect(depSanger,depBroad)
  
  FCs<-corrected[,id[1]]
  names(FCs)<-rownames(corrected)
  RES<-PrRc_Curve(FCs,essential_genes[,1],Nonessential_genes[,1],display = FALSE,FDRth = 0.05)
  depBroad<-names(which(FCs<=RES$sigthreshold))
  
  FCs<-corrected[,id[2]]
  names(FCs)<-rownames(corrected)
  RES<-PrRc_Curve(FCs,essential_genes[,1],Nonessential_genes[,1],display = FALSE,FDRth = 0.05)
  depSanger<-names(which(FCs<=RES$sigthreshold))
  
  BroadOnly_corrected[[i]]<-setdiff(depBroad,depSanger)
  SangerOnly_corrected[[i]]<-setdiff(depSanger,depBroad)
  common_corrected[[i]]<-intersect(depSanger,depBroad)
  
}

preCorrection<-cbind(unlist(lapply(BroadOnly,'length')),unlist(lapply(common,'length')),unlist(lapply(SangerOnly,'length')))
colnames(preCorrection)<-c('broad','both','sanger')
rownames(preCorrection)<-uc

print(paste("Average number common essentials pre correction:",mean(preCorrection[,"both"])))


postCorrection<-cbind(unlist(lapply(BroadOnly_corrected,'length')),unlist(lapply(common_corrected,'length')),unlist(lapply(SangerOnly_corrected,'length')))
colnames(postCorrection)<-c('broad','both','sanger')
rownames(postCorrection)<-uc
print(paste("Average number common essentials post correction:",mean(postCorrection[,"both"])))
O<-order(preCorrection[,2])

#Cohen's kappa for inter-rater agreement

PreCorrectionMat<-matrix(c(sum(preCorrection[,"both"]),sum(preCorrection[,"sanger"]),sum(preCorrection[,"broad"]),nrow(mergedDataset)*147-sum(unlist(preCorrection))),nrow=2)
PostCorrectionMat<-matrix(c(sum(postCorrection[,"both"]),sum(postCorrection[,"sanger"]),sum(postCorrection[,"broad"]),nrow(mergedDataset)*147-sum(unlist(postCorrection))),nrow=2)

K_pre<-Kappa(PreCorrectionMat)
K_post<-Kappa(PostCorrectionMat)

## Computing dependency agreement across cell-lines before/after batch correction
preAgreement<-preCorrection[,2]/rowSums(preCorrection)
postAgreement<-postCorrection[,2]/rowSums(postCorrection)
#CHECKED OK:
print(paste('Average study agreement pre correction:', 100*(round(mean(preAgreement),4)),'%',sep=''))

print(paste('Average study agreement batch corrected :', 100*(round(mean(postAgreement),4)),'%',sep=''))


## Estimating data quality before/after batch correction
preDataQ<-dataQuality(mergedDataset)
postDataQ<-dataQuality(corrected)

write.csv(preDataQ,file="./ResultsFolder/PreCorrectionDataQuality.csv")
write.csv(postDataQ,file="./ResultsFolder/PostCorrectionDataQuality.csv")


## Plotting screen agreement in relation to data quality
###Supplementary Figure 3 c ### 
pdf(paste0(ResultsFolder,"Unprocessed_batchcorrected_ScreenAgreement.pdf"),width=3,height=3)

par(mfrow=c(1,1))
par(mar=c(2.5,2.5,0.1,0.2)+0.1,mgp=c(1,0.25,0))
plot(colMeans(rbind(preDataQ[1:147],preDataQ[148:294])),
     100*preAgreement,xlab='Average data quality score (QS)',
     ylab='Screen Agreement %',col='darkgreen',
     ylim=c(0,100),xlim=c(0,1),pch=21,bg=makeTransparent('darkgreen'),cex.lab=0.9,cex.axis=0.8)
abline(lm(100*preAgreement ~ colMeans(rbind(preDataQ[1:147],preDataQ[148:294]))),col='darkgreen',lwd=2)
points(colMeans(rbind(preDataQ[1:147],preDataQ[148:294])),
       100*postAgreement,col='cyan',
       ylim=c(0,100),xlim=c(0,1),pch=21,bg=makeTransparent('cyan'))
abline(lm(100*postAgreement ~ colMeans(rbind(preDataQ[1:147],preDataQ[148:294]))),col='cyan',lwd=2)
abline(0,100,col='gray')
legend('topleft',col=c(makeTransparent('darkgreen'),
                       makeTransparent('cyan'),'gray'),legend=c('Merged dataset','Batch corrected','x=y'),
       lty=1,lwd=c(3,3,1),bty = 'n')
dev.off()

lmPre<-lm(100*preAgreement ~ colMeans(rbind(preDataQ[1:147],preDataQ[148:294])))
coefPvalPre<-summary(lmPre)$coefficients[2,4]
lmPost<-lm(100*postAgreement ~ colMeans(rbind(preDataQ[1:147],preDataQ[148:294])))
coefPvalPots<-summary(lmPost)$coefficients[2,4]

write.csv(corrected,file=paste0(ResultsFolder,"Unprocessed_CombatMerged.csv"),quote=F)

