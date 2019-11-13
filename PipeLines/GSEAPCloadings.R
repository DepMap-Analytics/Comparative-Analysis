
kegg<-read.gmt(file=paste0(ExternalData,"/kegg.gmt"))


pcloadings<-read.csv(file=paste0(ExternalData,"/gene_PC_loadings.csv"),header=T,stringsAsFactors = F,row.names=1)
pc1<-pcloadings[,1]
names(pc1)<-rownames(pcloadings)
rankPC1a<-sort(abs(pc1),decreasing=TRUE)
rankPC1<-sort(pc1,decreasing=TRUE)
pc2<-pcloadings[,2]
names(pc2)<-rownames(pcloadings)
rankPC2<-sort(pc2,decreasing=TRUE)


PC1kegg<-lapply(kegg,function(x) GSEAfunction(names(rankPC1),x))



PC2kegg<-lapply(kegg,function(x) GSEAfunction(names(rankPC2),x))
set.seed(31243)

#random permutation test can be loaded from ResultsFolder
RandPC1<-list()
RandPC2<-list()
for(i in 1:length(kegg)){
  RandPC1[[i]]<-list()
  RandPC2[[i]]<-list()
  for(j in 1:1000){
    RandPC1[[i]][[j]]<-GSEAfunction(sample(names(rankPC1)),kegg[[i]])[[1]]
    RandPC2[[i]][[j]]<-GSEAfunction(sample(names(rankPC2)),kegg[[i]])[[1]]

  }
}
#save(RandPC1,file=paste0(ResultsFolder,"PC1_Random.Rdata"))
#save(RandPC2,file=paste0(ResultsFolder,"PC2_Random.Rdata"))
load(file=paste0(ResultsFolder,"PC1_Random.Rdata"))
load(file=paste0(ResultsFolder,"PC2_Random.Rdata"))
pvalPC1<-list()
pvalPC2<-list()
for(i in 1:length(kegg)){

    pvalPC1[[i]]<-sum(abs(unlist(RandPC1[[i]]))>=abs(PC1kegg[[i]][[1]]))/1000

    pvalPC2[[i]]<-sum(abs(unlist(RandPC2[[i]]))>=abs(PC2kegg[[i]][[1]]))/1000

  
}
pvalPC1<-unlist(pvalPC1)
names(pvalPC1)<-names(kegg)
pvalPC2<-unlist(pvalPC2)
names(pvalPC2)<-names(kegg)

sigPC1<-names(pvalPC1[pvalPC1==0])
sigPC2<-names(pvalPC2[pvalPC2==0])

ciKeggPC1<-lapply(RandPC1,function(x) sd(unlist(x)))
ciKeggPC2<-lapply(RandPC2,function(x) sd(unlist(x)))
names(ciKeggPC1)<-names(kegg)
names(ciKeggPC2)<-names(kegg)

colPC1<-rainbow_hcl(length(sigPC1))
#supplementary Figures 
pdf(paste0(ResultsFolder,"/GSEA_PC1.pdf"))
for(i in 1:length(sigPC1)){
  if(i==1){
    plot(PC1kegg[[sigPC1[1]]]$RunningSum,col=colPC1[1],ylim=c(-0.5,0.5),type="l",ylab="Enrichment Score Running Sum",lwd=2,xlab="Gene Rank")
  }else{
    lines(PC1kegg[[sigPC1[i]]]$RunningSum,col=colPC1[i],ylim=c(-0.5,0.5),lwd=2)
    
  }
}
legend("topright",legend=sigPC1,col=colPC1,cex=0.5,lty=1)
dev.off()

colPC2<-rainbow_hcl(length(sigPC2))
pdf(paste0(ResultsFolder,"/GSEA_PC2.pdf"))
for(i in 1:length(sigPC2)){
  if(i==1){
    plot(PC2kegg[[sigPC2[1]]]$RunningSum,col=colPC2[1],ylim=c(-0.75,0.75),type="l",ylab="Enrichment Score Running Sum",lwd=2,xlab="Gene Rank")
  }else{
    lines(PC2kegg[[sigPC2[i]]]$RunningSum,col=colPC2[i],ylim=c(-0.75,0.75),lwd=2)
    
  }
}
legend("topright",legend=sigPC2,col=colPC2,cex=0.5,lty=1)
dev.off()