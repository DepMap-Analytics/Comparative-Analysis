
broad_data<-read.csv(paste0(ResultsFolder,"geneFCBroad.csv"),header=T,row.names=1,stringsAsFactors = F)
sanger_data<-read.csv(paste0(ResultsFolder,"geneFCSanger.csv"),header=T,row.names=1,stringsAsFactors = F)


## Selecting genes screened in both studies
cgenes<-intersect(rownames(broad_data),rownames(sanger_data))

## Attaching the study of origin identifier to each sample
colnames(broad_data)<-paste(colnames(broad_data),'---Broad',sep='')
colnames(sanger_data)<-paste(colnames(sanger_data),'---Sanger',sep='')

## Collating the two dataset on the domain of the genes screened in both studies
mergedDataset<-
  cbind(broad_data[cgenes,],sanger_data[cgenes,])



## Loading GO:terms collection
data(GOterms)

## Storing identifiers of screened cell lines in a vector
CLs<-unlist(lapply(str_split(colnames(mergedDataset),'---'),function(x) x[1]))
CLs<-CLs[1:147]

## Performing GO:enrichment analysis of study-exclusive dependencies, saving plots and tables

AllDisagreementCarCP<-function (dataset, CLs, GOterms, ref_Essential, ref_nonEssential, 
          ncbroad = 39, ncsanger = 20, filePrefix1, filePrefix2) 
{
  ncls <- length(CLs)
  FoundGenes <- lapply(1:length(GOterms$NAMES), function(x) "")
  EnrichmentsBroad <- list()
  EnrichmentsSanger <- list()
  for (i in 1:ncls) {
    print(paste("Characterising the study-exclusive dependencies found in cell line: ", 
                CLs[i], " (", i, " out of ", ncls, ")", sep = ""))
    RES <- DisCar(dataset = dataset, CL = CLs[i], ref_Essential, 
                  ref_nonEssential)
    tmpBroad <- RES$RES_Broad$Adj.p[RES$RES_Broad$Adj.p < 
                                      0.05]
    names(tmpBroad) <- RES$RES_Broad$GOterm[RES$RES_Broad$Adj.p < 
                                              0.05]
    EnrichmentsBroad[[i]] <- tmpBroad
    tmpSanger <- RES$RES_Sanger$Adj.p[RES$RES_Sanger$Adj.p < 
                                        0.05]
    names(tmpSanger) <- RES$RES_Sanger$GOterm[RES$RES_Sanger$Adj.p < 
                                                0.05]
    EnrichmentsSanger[[i]] <- tmpSanger
  }
  toSave <- qpcR:::cbind.na(names(EnrichmentsBroad[[1]]), EnrichmentsBroad[[1]])
  for (i in 2:length(CLs)) toSave <- qpcR:::cbind.na(toSave, 
                                                     qpcR:::cbind.na(names(EnrichmentsBroad[[i]]), EnrichmentsBroad[[i]]))
  rownames(toSave) <- NULL
  colnames(toSave) <- c(rbind(CLs, rep("", length(CLs))))
  write.table(toSave, file = paste(filePrefix1, ".txt"), sep = "\t", 
              quote = FALSE, row.names = FALSE)
  toSave <- qpcR:::cbind.na(names(EnrichmentsSanger[[1]]), 
                            EnrichmentsSanger[[1]])
  for (i in 2:length(CLs)) toSave <- qpcR:::cbind.na(toSave, 
                                                     qpcR:::cbind.na(names(EnrichmentsSanger[[i]]), EnrichmentsSanger[[i]]))
  rownames(toSave) <- NULL
  colnames(toSave) <- c(rbind(CLs, rep("", length(CLs))))
  write.table(toSave, file = paste(filePrefix2, ".txt"), sep = "\t", 
              quote = FALSE, row.names = FALSE)
  enrichmentMapBroad <- do.call(cbind, lapply(EnrichmentsBroad, 
                                              function(x) {
                                                sigPattenrn <- rep(1, length(GOterms$NAMES))
                                                ii <- which(is.element(GOterms$NAMES, names(x)))
                                                sigPattenrn[ii] <- x
                                                return(sigPattenrn)
                                              }))
  rownames(enrichmentMapBroad) <- GOterms$NAMES
  colnames(enrichmentMapBroad) <- CLs
  enrichmentMapSanger <- do.call(cbind, lapply(EnrichmentsSanger, 
                                               function(x) {
                                                 sigPattenrn <- rep(1, length(GOterms$NAMES))
                                                 ii <- which(is.element(GOterms$NAMES, names(x)))
                                                 sigPattenrn[ii] <- x
                                                 return(sigPattenrn)
                                               }))
  rownames(enrichmentMapSanger) <- GOterms$NAMES
  colnames(enrichmentMapSanger) <- CLs
  enrichmentMapBroad <- -log10(enrichmentMapBroad)
  enrichmentMapSanger <- -log10(enrichmentMapSanger)
  enrichmentMapBroad <- enrichmentMapBroad[which(rowSums(enrichmentMapBroad > 
                                                           0) > 0), ]
  enrichmentMapSanger <- enrichmentMapSanger[which(rowSums(enrichmentMapSanger > 
                                                             0) > 0), ]
  enrichmentMapBroad <- enrichmentMapBroad[order(apply(enrichmentMapBroad, 
                                                       MARGIN = 1, "median"), decreasing = TRUE), ]
  enrichmentMapSanger <- enrichmentMapSanger[order(apply(enrichmentMapSanger, 
                                                         MARGIN = 1, "median"), decreasing = TRUE), ]
  pdf(paste(filePrefix1, ".pdf", sep = ""), width=8, height=4)
  par(mfrow = c(1, 1))
  #par(mar = c(5, 25, 4, 1))
  par(mar=c(2.5,15,2,2)+0.1,mgp=c(1.2,0.4,0.2))
  toPlot <- enrichmentMapBroad[which(rowSums(enrichmentMapBroad > 
                                               0) > ncbroad), ]
  toBarplot <- rowSums(toPlot > 0)
  broadToBarplot <- toBarplot
  O <- order(toBarplot)
  names(toBarplot) <- str_replace_all(names(toBarplot), "_", 
                                      " ")
  names(toBarplot) <- str_replace(names(toBarplot), " ", ":")
  names(toBarplot) <- str_replace(names(toBarplot), "GO:", " ")
  #names(toBarplot)<-tolower(names(toBarplot))
  print(names(toBarplot))
  toplotBroad<-toBarplot[0]/ncls
  barplot((100 * toBarplot[O]/ncls), horiz = TRUE, las = 2, cex.names = 0.6, 
          xlim = c(0, 100), col = "#80479C", border = FALSE, main =  "GO:Terms enriched in Broad-exclusive \n depleted gene sets", 
          xlab = "% cell lines",cex.axis = 0.8,cex.main=0.8,cex.lab=0.6)
  #par(mar = c(5, 0, 4, 4))
  #boxplot((t(toPlot[O, ])), las = 2, horizontal = TRUE, cex.axis = 1, 
  #        names = rep("", length(O)), frame.plot = FALSE, xlab = "-log10 (enrichment p)")
  dev.off()
  pdf(paste(filePrefix2, ".pdf", sep = ""), width=4, height=4)
  par(mfrow = c(1, 1))
  #par(mar = c(5, 25, 4, 1))
  par(mar=c(2.5,9.5,2,2)+0.1,mgp=c(1.2,0.4,0.2))
  toPlot <- enrichmentMapSanger[which(100 * rowSums(enrichmentMapSanger > 
                                                      0)/ncls > ncsanger), ]
  nEnrichments <- length(which(rowSums(enrichmentMapSanger > 
                                         0) > 0))
  toBarplot <- rowSums(toPlot > 0)
  O <- order(toBarplot)
  names(toBarplot) <- str_replace_all(names(toBarplot), "_", 
                                      " ")
  names(toBarplot) <- str_replace(names(toBarplot), " ", ":")
  names(toBarplot) <- str_replace(names(toBarplot), "GO:", " ")
  barplot((100 * toBarplot[O]/ncls), horiz = TRUE, las = 2, cex.names = 0.6, 
          xlim = c(0, 100), col = "#3A53A4", border = FALSE, main = "GO:Terms enriched in Sanger-exclusive \n depleted gene sets", 
          xlab = "% cell lines",cex.axis = 0.8,cex.main=0.8,cex.lab=0.6)
  #par(mar = c(5, 0, 4, 4))
  #boxplot((t(toPlot[O, ])), las = 2, horizontal = TRUE, cex.axis = 1, 
          #names = rep("", length(O)), frame.plot = FALSE, xlab = "-log10 (enrichment p)")
  dev.off()
  return(list(EMB = enrichmentMapBroad, EMS = enrichmentMapSanger,toplotBroad=toplotBroad))
}


#Filter genes to remove those with significantly different sgRNA - assume these are library effects
## Storing cell line name and study of origin for each sample into two vectors
clnames<-str_split(colnames(mergedDataset),'---')
clnames<-unlist(lapply(clnames,function(x){x[1]}))

## identifying significant dependencies (at 5% FDR) before batch corretion, across cell lines
uc<-unique(clnames)
ncl<-length(uc)

BroadOnly<-list()
SangerOnly<-list()
common<-list()

for (i in 1:ncl){
  #print(c(clnames[i],i,ncl))
  
  id<-which(clnames==uc[i])
  
  FCs<-mergedDataset[,id[1]]
  names(FCs)<-rownames(mergedDataset)
  RES<-PrRc_Curve(FCs,BAGEL_essential,BAGEL_nonEssential,display = FALSE,FDRth = 0.05)
  depBroad<-names(which(FCs<=RES$sigthreshold))
  
  FCs<-mergedDataset[,id[2]]
  names(FCs)<-rownames(mergedDataset)
  RES<-PrRc_Curve(FCs,BAGEL_essential,BAGEL_nonEssential,display = FALSE,FDRth = 0.05)
  depSanger<-names(which(FCs<=RES$sigthreshold))
  
  BroadOnly[[i]]<-setdiff(depBroad,depSanger)
  SangerOnly[[i]]<-setdiff(depSanger,depBroad)
  common[[i]]<-intersect(depSanger,depBroad)
  
}


Broad_eff<-read.csv(paste0(InputFolder,"GuideMap_Avana.csv"))
Broad_eff$efficacy<-Broad_eff$efficacy-mean(Broad_eff$efficacy,na.rm=TRUE)
Broad_g_eff<-aggregate(Broad_eff$efficacy,by=list(Broad_eff$gene),FUN='sd',na.rm=TRUE)
nn<-Broad_g_eff$Group.1
Broad_g_eff<-Broad_g_eff$x
names(Broad_g_eff)<-nn

Sanger_eff<-read.csv(paste0(InputFolder,"GuideMap_KY10.csv"))
Sanger_eff$efficacy<-Sanger_eff$efficacy-mean(Sanger_eff$efficacy,na.rm=TRUE)
Sanger_g_eff<-aggregate(Sanger_eff$efficacy,by=list(Sanger_eff$gene),FUN='sd',na.rm=TRUE)
nn<-Sanger_g_eff$Group.1
Sanger_g_eff<-Sanger_g_eff$x
names(Sanger_g_eff)<-nn

allBroadExclusive<-unique(unlist(BroadOnly))
bo_map<-do.call(cbind,lapply(BroadOnly,function(x){is.element(allBroadExclusive,x)}))+0
rownames(bo_map)<-allBroadExclusive

allSangerExclusive<-unique(unlist(SangerOnly))
so_map<-do.call(cbind,lapply(SangerOnly,function(x){is.element(allSangerExclusive,x)}))+0
rownames(so_map)<-allSangerExclusive

soGenes<-names(which((rowSums(so_map)/ncl)>=0.50))
boGenes<-names(which((rowSums(bo_map)/ncl)>=0.50))

par(mfrow=c(2,2))
boxplot(Broad_g_eff[boGenes],Broad_g_eff[setdiff(names(Broad_g_eff),boGenes)],
        Sanger_g_eff[boGenes],Sanger_g_eff[setdiff(names(Sanger_g_eff),boGenes)],
        outline=FALSE,
        main='Broad exclusive dependencies',col=makeTransparent(c('purple','pink','darkblue','blue')),
        ylab='Mean norm. sgRNA efficiency')
t.test(Broad_g_eff[boGenes],Broad_g_eff[setdiff(names(Broad_g_eff),boGenes)])$p.value
cohens_d(Broad_g_eff[boGenes],Broad_g_eff[setdiff(names(Broad_g_eff),boGenes)])
t.test(Sanger_g_eff[boGenes],Sanger_g_eff[setdiff(names(Sanger_g_eff),boGenes)])$p.value
cohens_d(Sanger_g_eff[boGenes],Sanger_g_eff[setdiff(names(Sanger_g_eff),boGenes)])

boxplot(Broad_g_eff[soGenes],Broad_g_eff[setdiff(names(Broad_g_eff),soGenes)],
        Sanger_g_eff[soGenes],Sanger_g_eff[setdiff(names(Sanger_g_eff),soGenes)],
        outline=FALSE,
        main='Sanger exclusive dependencies',col=makeTransparent(c('purple','pink','darkblue','blue')),
        ylab='Mean norm. sgRNA efficiency')
t.test(Broad_g_eff[soGenes],Broad_g_eff[setdiff(names(Broad_g_eff),soGenes)])$p.value
cohens_d(Broad_g_eff[soGenes],Broad_g_eff[setdiff(names(Broad_g_eff),soGenes)])
t.test(Sanger_g_eff[boGenes],Sanger_g_eff[setdiff(names(Sanger_g_eff),soGenes)])$p.value
cohens_d(Sanger_g_eff[boGenes],Sanger_g_eff[setdiff(names(Sanger_g_eff),soGenes)])

#Remove all sanger only or broad only genes that have significantly different efficiency.

#Null distributions is there is no difference in sgRNA efficiency. Test all the broad only differences against this to find signif difference in efficiency
#Null distn differences are genes across institutes that are not in either institute exclusive sets
soGenes2<-names(which((rowSums(so_map)/ncl)>=0.25))
boGenes2<-names(which((rowSums(bo_map)/ncl)>=0.25))
exclusivegenes<-c(boGenes2,soGenes2)
nullgenes<-setdiff(names(Broad_g_eff),exclusivegenes)
nullgene_diff<-Broad_g_eff[nullgenes]-Sanger_g_eff[nullgenes]
mean_null<-mean(nullgene_diff,na.rm=T)
sd_null<-sd(nullgene_diff,na.rm=T)
Bonly_diff<-Broad_g_eff[boGenes2]-Sanger_g_eff[boGenes2]
Bonly_diff<-na.omit(Bonly_diff)
min_cutoff<-mean_null-2*sd_null
max_cutoff<-mean_null+2*sd_null
remove_genes<-names(Bonly_diff)[which(Bonly_diff>max_cutoff|Bonly_diff<min_cutoff)]


Sonly_diff<-Broad_g_eff[soGenes2]-Sanger_g_eff[soGenes2]
Sonly_diff<-na.omit(Sonly_diff)

remove_genes2<-names(Sonly_diff)[which(Sonly_diff>max_cutoff|Sonly_diff<min_cutoff)]

###Figure 5 c ####

ALLenrichTP_allHalf<-AllDisagreementCarCP(dataset = mergedDataset[setdiff(rownames(mergedDataset),intersect(remove_genes,remove_genes2)),],
                                      CLs = CLs,
                                      GOterms = GOterms,
                                      ref_Essential = BAGEL_essential,
                                      ref_nonEssential = BAGEL_nonEssential,ncbroad=74,ncsanger=44,
                                      filePrefix1 = paste0(ResultsFolder,'GOenrichBroad_TP_allHalf'),
                                      filePrefix2 = paste0(ResultsFolder,'GOenrichSanger_TP_allHalf'))

