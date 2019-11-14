my.hypTest<-function(x,k,n,N){


  PVALS<-phyper(x-1,n,N-n,k,lower.tail=FALSE)

  return(PVALS)
}
decodeCNA<-function(MoBEM){
  data(CNAdecode)

  rn<-rownames(MoBEM)

  ii<-grep('cna',rownames(MoBEM))

  cnaId<-rownames(MoBEM)[ii]

  containedGenes<-unlist(lapply(str_split(cnaId,' '),function(x) {x[2]}))

  containedGenes[is.na(containedGenes)]<-''

  segments<-unlist(lapply(str_split(unlist(lapply(str_split(cnaId,':'),function(x) {x[2]})),' '),function(x){x[1]}))

  loci<-as.character(CNAdecode$locus[match(segments,CNAdecode$Identifier)])
  altType<-as.character(CNAdecode$Recurrent.Amplification..Deletion[match(segments,CNAdecode$Identifier)])

  altType[altType=='Amplification']<-'G:'
  altType[altType=='Deletion']<-'L:'

  rownames(MoBEM)[ii]<-paste(altType,loci,' ',containedGenes,sep='')
  return(MoBEM)
}
diffDep<-function(MoBEM,site,dataset,CFE,s,display=TRUE,labels=TRUE){

  current_cids<-cids[which(site==s)]
  cfePattern<-MoBEM[CFE,as.character(current_cids)]

  subdataset<-dataset[,which(site==s)]

  subSetWT<-subdataset[,cfePattern==0]
  subSetALT<-subdataset[,cfePattern==1]

  effect_size<-
    unlist(lapply(seq_len(nrow(subSetALT)), function(i) {
      cohens_d(subSetALT[i,],subSetWT[i,])
    }))

  delta<-
    unlist(lapply(seq_len(nrow(subSetALT)), function(i) {
      mean(subSetALT[i,])-mean(subSetWT[i,])
    }))

  pval<-
    unlist(lapply(seq_len(nrow(subSetALT)), function(i) {
      RES<-t.test(subSetALT[i,],subSetWT[i,],var.equal = TRUE)
      RES$p.value}))

  FDR<-p.adjust(pval,'fdr')

  if(display){

    pval_at_0.10FDR<-max(pval[which(FDR<=max(FDR[which(FDR<0.10)]))])
    pval_at_0.05FDR<-max(pval[which(FDR<=max(FDR[which(FDR<0.05)]))])
    pval_at_0.01FDR<-max(pval[which(FDR<=max(FDR[which(FDR<0.01)]))])

    col<-rep('darkgray',length(pval))
    col[which(effect_size > 1 & FDR<0.10 & delta < -1)]<-'purple'
    col[which(effect_size > 1 & FDR<0.10 & delta > 1)]<-'orange'

    plot(effect_size*sign(delta),-log10(pval),
         bg=makeTransparent(col),pch=21,col=col,frame.plot=FALSE,
         xlab='signed effect size')
    abline(v=0,col='darkgray')
    abline(h=-log10(pval_at_0.10FDR),col='darkgray',lty=3)
    abline(h=-log10(pval_at_0.05FDR),col='darkgray',lty=2)
    abline(h=-log10(pval_at_0.01FDR),col='darkgray',lty=1)
    abline(v=c(-1,1),col='gray')

    if(labels){
      id<-which(FDR<0.10 & effect_size > 1 & abs(delta)>1)
      if(length(id)>0){
        text(effect_size[id]*sign(delta[id]),-log10(pval[id]),rownames(subdataset)[id],cex=0.5,pos = 4)
      }
    }

  }


  res<-data.frame(CFE=CFE,GENE=rownames(dataset),
                  delta=delta,
                  effect_size=effect_size,
                  p=pval,
                  fdr=FDR,
                  stringsAsFactors = FALSE)
  return(res)
}
cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1

  md  <- abs(mean(x,na.rm = TRUE) - mean(y,na.rm=TRUE))        ## mean difference (numerator)
  csd <- lx * var(x,na.rm = TRUE) + ly * var(y,na.rm=TRUE)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation

  cd  <- md/csd                        ## cohen's d

  return(cd)
}
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
compareInteractionPlot<-function(loc_RES_S,loc_RES_B,IDENTIFY=FALSE,
                                 FDRth=0.01,
                                 diffDepTh= -1,
                                 ESth=1,
                                 XLIM=NULL,YLIM=NULL,
                                 study1name='study1',study2name='study2'){
  options(warn=-1)
  col<-rep('gray',nrow(loc_RES_S))

  tth<-FDRth
  #iSanger_s<-which(loc_RES_S$fdr_oa<quantile(loc_RES_S$fdr_oa,na.rm = TRUE,tth))
  #iBroad_s<-which(loc_RES_B$fdr_oa<quantile(loc_RES_B$fdr_oa,na.rm = TRUE,tth))

  iSanger_s<-which(loc_RES_S$fdr_oa<tth & loc_RES_S$effect_size > ESth & loc_RES_S$delta < diffDepTh)
  iBroad_s<-which(loc_RES_B$fdr_oa<tth & loc_RES_B$effect_size > ESth & loc_RES_B$delta < diffDepTh)

  iboth_s<-intersect(iSanger_s,iBroad_s)
  iEither<-union(iSanger_s,iBroad_s)

  col[iSanger_s]<-'darkgray'
  col[iBroad_s]<-'darkgray'
  col[iboth_s]<-'purple'

  if(length(XLIM)==0){
    XLIM<-range(loc_RES_S$delta,na.rm = TRUE)
  }

  if(length(YLIM)==0){
    YLIM<-range(loc_RES_B$delta,na.rm = TRUE)
  }


  #par(mfrow=c(1,2))
  smoothScatter(loc_RES_S$delta,
                loc_RES_B$delta,
                nrpoints = 0,
                colramp = colorRampPalette(colors=c('white','black')),nbin = 100,
                xlab=paste('delta log FC [',study1name,']',sep=''),
                ylab=paste('delta log FC [',study2name,']',sep=''),
                xlim=XLIM,ylim=YLIM)

  #plot(hexbin(loc_RES_S$delta,loc_RES_B$delta))

  points(loc_RES_S$delta[c(iboth_s,setdiff(iSanger_s,iBroad_s),setdiff(iBroad_s,iSanger_s))],
         loc_RES_B$delta[c(iboth_s,setdiff(iSanger_s,iBroad_s),setdiff(iBroad_s,iSanger_s))],
         cex=1.5,pch=21,bg=makeTransparent(col[c(iboth_s,setdiff(iSanger_s,iBroad_s),setdiff(iBroad_s,iSanger_s))]),
         col=col[c(iboth_s,setdiff(iSanger_s,iBroad_s),setdiff(iBroad_s,iSanger_s))])

  abline(h= -1,col='gray',lty=2)
  abline(v= -1,col='gray',lty=2)

  abline(h= 0,col='gray')
  abline(v= 0,col='gray')


  if(IDENTIFY){
    identify(loc_RES_S$delta[c(iboth_s,setdiff(iSanger_s,iBroad_s),setdiff(iBroad_s,iSanger_s))],
             loc_RES_B$delta[c(iboth_s,setdiff(iSanger_s,iBroad_s),setdiff(iBroad_s,iSanger_s))],
             paste(loc_RES_B$CFE[c(iboth_s,setdiff(iSanger_s,iBroad_s),setdiff(iBroad_s,iSanger_s))],
                   loc_RES_B$GENE[c(iboth_s,setdiff(iSanger_s,iBroad_s),setdiff(iBroad_s,iSanger_s))],
                   sep='\n'),
             cex=0.5)
  }
  options(warn=0)

  N<-nrow(RES_S)
  n<-length(iSanger_s)
  k<-length(iBroad_s)
  x<-length(iboth_s)


  print(paste('total number of performed tests:',N))
  print(paste('n. of significant associations in the first study at ',100*tth,'% FDR',': ',n,sep=''))
  print(paste('n. of significant associations in the second study at ',100*tth,'% FDR',': ',k,sep=''))
  print(paste('n. of significant associations shared by the two studies at ',100*tth,'% FDR',': ',x,sep=''))
  print(paste(format(100*x/n,digits=4),'% of those in the first study'))
  print(paste(format(100*x/k,digits=4),'% of those in the second study'))
  print(paste('Fisher exact test p =',my.hypTest(x,k,n,N)))

  print(paste('Correlation betwenn differential essentialities = ',
              format(cor(loc_RES_S$delta,loc_RES_B$delta,use = 'complete'),
                     digits=3)))


# smoothScatter(RES_S$effect_size,RES_B$effect_size,nrpoints = 0,
#                colramp = colorRampPalette(colors=c('white','black')),
#                nbin = 100,
#                xlab=paste('Signed effect size [',study1name,']',sep=''),
#                ylab=paste('Signed effect size [',study2name,']',sep=''))

#  concoPoints<-which(RES_B$fdr_oa<tth | RES_S$fdr_oa<tth)
#  points(RES_S$signed_es[concoPoints],RES_B$signed_es[concoPoints],
#         cex=1.5,pch=21,bg=makeTransparent('darkgray',120),
#         col='darkgray')

#  concoPoints<-which(RES_B$fdr_oa<tth & RES_S$fdr_oa<tth)
#  points(RES_S$signed_es[concoPoints],RES_B$signed_es[concoPoints],
#         cex=1.5,pch=21,bg=makeTransparent('darkgreen'),
#         col='darkgreen')
#  abline(v=0,col='black')
#  abline(h=0,col='black')

  XLAB<-RES_S$delta
  YLAB<-RES_B$delta

  xloc<-max(range(XLAB[XLAB>0],na.rm = TRUE))/2+
    min(XLAB[XLAB>0],na.rm = TRUE)
  yloc<-max(range(YLAB[YLAB>0],na.rm = TRUE))/2+
    min(YLAB[YLAB>0],na.rm = TRUE)

  xlocneg<-min(range(XLAB[XLAB<0],na.rm = TRUE))/2
  ylocneg<-min(range(YLAB[YLAB<0],na.rm = TRUE))/2

  negneg<-length(which(RES_B$fdr_oa<FDRth & RES_S$fdr_oa<FDRth &
                         RES_B$delta<0 & RES_S$delta<0))

  negnegall<-length(which((RES_B$fdr_oa<FDRth | RES_S$fdr_oa<FDRth) &
                            RES_B$delta<0 & RES_S$delta<0))

  negpos<-length(which(RES_B$fdr_oa<FDRth & RES_S$fdr_oa<FDRth &
                         RES_B$delta<0 & RES_S$delta>0))

  negposall<-length(which((RES_B$fdr_oa<FDRth | RES_S$fdr_oa<FDRth) &
                            RES_B$delta<0 & RES_S$delta>0))

  posneg<-length(which(RES_B$fdr_oa<FDRth & RES_S$fdr_oa<FDRth &
                         RES_B$delta>0 & RES_S$delta<0))

  posnegall<-length(which((RES_B$fdr_oa<tth | RES_S$fdr_oa<tth) &
                            (RES_B$signed_es>0 & RES_S$signed_es<0)))

  pospos<-length(which(RES_B$fdr_oa<FDRth & RES_S$fdr_oa<FDRth &
                         RES_B$delta>0 & RES_S$delta>0))
  posposall<-length(which((RES_B$fdr_oa<tth | RES_S$fdr_oa<tth) &
                            (RES_B$signed_es>0 & RES_S$signed_es>0)))

  # smoothScatter(loc_RES_S$delta,
  #               loc_RES_B$delta,
  #               nrpoints = 0,
  #               colramp = colorRampPalette(colors=c('white','black')),nbin = 100,
  #               xlab=paste('delta log FC [',study1name,']',sep=''),
  #               ylab=paste('delta log FC [',study2name,']',sep=''),
  #               xlim=XLIM,ylim=YLIM)
  # text(xlocneg,ylocneg,paste(negneg,negnegall),col='black',cex=1.5)
  # text(xloc,yloc,paste(pospos,posposall),col='black',cex=1.5)
  # text(xlocneg,yloc,paste(negpos,negposall),col='black',cex=1.5)
  # text(xloc,ylocneg,paste(posneg,posnegall),col='black',cex=1.5)
  # abline(h=0)
  # abline(v=0)

  print(paste('Sign concordance of the associations shared by the studies: ',
        format(100*(pospos+negneg)/(pospos+negneg+negpos+posneg),digits = 3),'%',sep=''))
  print(paste('Fisher exact test p =',
              fisher.test(matrix(c(negpos,negneg,pospos,posneg),2,2))$p.value))

  text(xloc,ylocneg,paste(posneg,posnegall),col='black',cex=1.5)

  print(paste('Sign concordance of the associations detected in at least one study: ',
              format(100*(posposall+negnegall)/(posposall+negnegall+negposall+posnegall),
                     digits=3),'%',sep=''))

  print(paste('Fisher exact test p =',
              fisher.test(matrix(c(negposall,negnegall,posposall,posnegall),2,2))$p.value))

  print(paste('Correlation betwenn deltas = ',
              format(cor(loc_RES_S$delta,loc_RES_B$delta,use='complete'),digits=3)))
}


PrRc_and_ROC_Curves<-function(loc_RES_S,loc_RES_B,es_th=1,fdr_th=0.01,
                        study1name = 'study1',
                        study2name = 'study2'){

  WW<-loc_RES_S
  WA<-loc_RES_B

  ii<-which(WW$fdr_oa<fdr_th)
  O<-ii[order(WW$p[ii])]

  prctl<-floor(quantile(1:length(O),c(0.20,0.40,0.60,0.80,1)))

  cols<-makeTransparent(colorRampPalette(c('purple','orange'))(5),120)

  Npos<-vector()
  AUPREC<-vector()
  rndAUPREC<-vector()
  AUROC<-vector()
  rndAUROC<-rep(0.5,5)

  for (i in 1:length(prctl)){

    positives<-O[1:prctl[i]]
    Npos[i]<-length(positives)

    observed<-rep(0,nrow(WA))
    names(observed)<-paste('a',1:nrow(WA))
    observed[positives]<-1
    predictions<- -log10(WA$p)
    names(predictions)<-paste('a',1:nrow(WA))
    RES<-ccr.PrRc_Curve(FCsprofile = -predictions,
                        positives=names(which(observed>0)),
                        negatives=names(which(observed==0)),display = FALSE)
    AUPREC[i]<-RES$AUC
    rndAUPREC[i]<-RES$RND

    if(i == 1){
      plot(RES$curve[,"recall"],
           RES$curve[,"precision"],type='l',lwd=4,
           xlab='Recall',
           ylab='Precision',col=cols[i],
           main=paste(study1name,'associations\nvs',study2name,'ranked p-value'))
    }else{
      lines(RES$curve[,"recall"],
            RES$curve[,"precision"],
            col=cols[i],
            lwd=4)
    }

    abline(h=min(RES$curve[,'precision']),col=cols[i],lty=1)

  }

  legend('left','random',lty=1,bty = 'n')

  for (i in 1:length(prctl)){
    positives<-O[1:prctl[i]]
    observed<-rep(0,nrow(WA))
    names(observed)<-paste('a',1:nrow(WA))
    observed[positives]<-1
    predictions<- -log10(WA$p)
    names(predictions)<-paste('a',1:nrow(WA))
    RES<-ROC_curve(-predictions,
                   positives=names(which(observed>0)),
                   negatives=names(which(observed==0)))
    AUROC[i]<-RES$auc
    if(i == 1){
      plot(RES$specificities,
           RES$sensitivities,type='l',lwd=4,
           xlab='Specificity',
           ylab='Recall',col=cols[i],
           xlim=c(1,0))
    }else{
      lines(RES$specificities,
            RES$sensitivities,
            col=cols[i],
            lwd=4)
    }
  }
  lines(c(0,1),c(1,0),col='gray')

  legend('bottomright',c('20%','40%','60%','80%','100%'),title = 'Significance quantile',
         col=cols,lwd=4)

  names(Npos)<-paste(100*c(0.20,0.40,0.60,0.80,1),'%',sep='')
  names(AUPREC)<-paste(100*c(0.20,0.40,0.60,0.80,1),'%',sep='')
  names(rndAUPREC)<-paste(100*c(0.20,0.40,0.60,0.80,1),'%',sep='')
  names(AUROC)<-paste(100*c(0.20,0.40,0.60,0.80,1),'%',sep='')
  names(rndAUROC)<-paste(100*c(0.20,0.40,0.60,0.80,1),'%',sep='')

  RES<-cbind(Npos,AUPREC,rndAUPREC,AUROC,rndAUROC)

  return(RES)

}
ROC_curve<-function(FCsprofile,positives,negatives,display=TRUE,FDRth=NULL){

  FCsprofile<-FCsprofile[intersect(c(positives,negatives),names(FCsprofile))]

  predictions<-FCsprofile
  observations<-is.element(names(FCsprofile),positives)+0
  names(observations)<-names(predictions)

  RES<-roc(observations,predictions,direction = '>')

  return(RES)
}
PrRc_Curve<-function(FCsprofile,positives,negatives,display=TRUE,FDRth=NULL){

  FCsprofile<-FCsprofile[intersect(c(positives,negatives),names(FCsprofile))]

  predictions<- -FCsprofile
  observations<-is.element(names(FCsprofile),positives)+0
  names(observations)<-names(predictions)


  prc<-pr.curve(scores.class0 = predictions,weights.class0 = observations,
                curve = TRUE,sorted = TRUE)

  PRECISION<-prc$curve[,2]
  RECALL<-prc$curve[,1]

  if(display){
    plot(RECALL,PRECISION,col='blue',lwd=3,xlab='Recall',ylab='Precision',type='l',xlim=c(0,1),ylim=c(0,1))
  }

  SENS<-NULL
  threshold<-NULL
  if(length(FDRth)>0){

    FDR5percTh<- -prc$curve[min(which(prc$curve[,2]>= 1-FDRth)),3]
    SENS<- prc$curve[min(which(prc$curve[,2]>= 1-FDRth)),1]
    threshold<-FDR5percTh
    if(display){
      abline(h=1-FDRth,lty=2)

      abline(v=SENS,lty=1)
    }
  }

  if(display){
    if(length(SENS)==0){
      legend('bottomleft',paste('AUC = ',format(prc$auc.integral,digits=3)),bty = 'n')
    }else{
      legend('bottomleft',c(paste('Recall ',100*FDRth,'%FDR = ',format(SENS,digits=3),sep=''),
                            paste('AUC = ',format(prc$auc.integral,digits=3))),bty = 'n')
    }

    abline(h=sum(observations)/length(observations))
  }
  #
  curve<-prc$curve
  colnames(curve)<-c('recall','precision','threshold')
  RES<-list(AUC=prc$auc.integral,Recall=SENS,sigthreshold=threshold,curve=curve)
  # ### threshold, and recall at fixed FDR to be returned
  return(RES)
}
plotAssociationExamples<-function(gene,feat,
                                  loc_CombinedDataset,
                                  loc_site,
                                  loc_cids){
  essS<-loc_CombinedDataset[gene,loc_site=='Sanger']
  essB<-loc_CombinedDataset[gene,loc_site=='Broad']
  pattern<-MoBEM[feat,as.character(loc_cids[loc_site=='Sanger'])]

  cols<-rep('darkgray',length(essB))
  cols[which(pattern==1)]<-'darkgreen'

#  par(mfrow=c(1,2))

  beeswarm(essB~pattern,corral='wrap',bg=c(makeTransparent('gray'),makeTransparent('darkgreen')),
           pch=21,col=c('gray','darkgreen'),cex=1.5,ylim=range(essB),las=2,labels=c('absent','present'),
           xlab=feat,ylab='Broad fitness score')
  par(new=TRUE)
  boxplot(essB~pattern,col=NA,ylim=range(essB),frame.plot=FALSE,xaxt='n',yaxt='n',outline=FALSE)

  beeswarm(essS~pattern,corral='wrap',bg=c(makeTransparent('gray'),makeTransparent('darkgreen')),
           pch=21,col=c('gray','darkgreen'),cex=1.5,ylim=range(essS),las=2,labels=c('absent','present'),
           xlab=feat,ylab='Sanger fitness score')
  par(new=TRUE)
  boxplot(essS~pattern,col=NA,ylim=range(essS),frame.plot=FALSE,xaxt='n',yaxt='n',outline=FALSE)


}
decodeCNA<-function(MoBEM){
  load('../../___OpenTargets_p015/_BROADcompData/CNAdecode.RData')

  rn<-rownames(MoBEM)

  ii<-grep('cna',rownames(MoBEM))

  cnaId<-rownames(MoBEM)[ii]

  containedGenes<-unlist(lapply(str_split(cnaId,' '),function(x) {x[2]}))

  containedGenes[is.na(containedGenes)]<-''

  segments<-unlist(lapply(str_split(unlist(lapply(str_split(cnaId,':'),function(x) {x[2]})),' '),function(x){x[1]}))

  loci<-as.character(CNAdecode$locus[match(segments,CNAdecode$Identifier)])
  altType<-as.character(CNAdecode$Recurrent.Amplification..Deletion[match(segments,CNAdecode$Identifier)])

  altType[altType=='Amplification']<-'G:'
  altType[altType=='Deletion']<-'L:'

  rownames(MoBEM)[ii]<-paste(altType,loci,' ',containedGenes,sep='')
  return(MoBEM)
}
decodeCFEs<-function(CFEs){
  load('../../___OpenTargets_p015/_BROADcompData/CNAdecode.RData')

  rn<-CFEs

  ii<-grep('cna',rn)

  cnaId<-rn[ii]

  containedGenes<-unlist(lapply(str_split(cnaId,' '),function(x) {x[2]}))

  containedGenes[is.na(containedGenes)]<-''

  segments<-unlist(lapply(str_split(unlist(lapply(str_split(cnaId,':'),function(x) {x[2]})),' '),function(x){x[1]}))

  loci<-as.character(CNAdecode$locus[match(segments,CNAdecode$Identifier)])
  altType<-as.character(CNAdecode$Recurrent.Amplification..Deletion[match(segments,CNAdecode$Identifier)])

  altType[altType=='Amplification']<-'G:'
  altType[altType=='Deletion']<-'L:'

  rn[ii]<-paste(altType,loci,' ',containedGenes,sep='')
  return(rn)
}
nDepletions<-function(FCs){
  
  data(BAGEL_essential)
  data(BAGEL_nonEssential)
  
  RES<-ccr.PrRc_Curve(FCs,BAGEL_essential,BAGEL_nonEssential,display = FALSE,FDRth = 0.05)
  
  ndep<-length(which(FCs< RES$sigthreshold))
  
  return(list(Ndep=ndep,Recall=RES$Recall,sigDep=FCs<RES$sigthreshold+0))
}
tSNEplot<-function(dataset,nit=1000,title='',perplexity=30){
  clnames<-str_split(colnames(dataset),'---')
  clnames<-unlist(lapply(clnames,function(x){x[1]}))
  
  site<-str_split(colnames(dataset),'---')
  site<-unlist(lapply(site,function(x){x[[2]]}))
  
  colnames(dataset)<-clnames
  cdist<-as.dist(1-cor(dataset))
  
  set.seed(0xA5EED)
  
  tfit<-tsne(cdist,max_iter = nit,perplexity=perplexity)
  
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  uc<-unique(clnames)
  COLS<-sample(color,length(uc))
  names(COLS)<-uc
  
  tCOLS<-makeTransparent(COLS,alpha = 140)
  
  SYMBOLS<-rep(21,length(site))
  SYMBOLS[site=='Broad']<-23
  names(SYMBOLS)<-site
  CEX<-rep(1.3,length(site))
  CEX[site=='Broad']<-1.8
  
  plot(tfit[,1],tfit[,2],bg=tCOLS[clnames],col='black',cex=CEX,pch=SYMBOLS[site],frame.plot = FALSE,xaxt='n',yaxt='n',
       xlab='tSNE AU',ylab='tSNE AU',main=title)
  
}

multDensPlot<-function(TOPLOT,COLS,XLIMS,YLIMS,TITLE,LEGentries,XLAB=''){

  YM<-vector()
  for (i in 1:length(TOPLOT)){
    YM[i]<-max(TOPLOT[[i]]$y,na.rm = TRUE)
  }

  Ymax<-max(YM,na.rm=TRUE)

  plot(0,0,col=NA,ylab='density',xlab=XLAB,
       xlim=XLIMS,ylim=YLIMS,type='l',main=TITLE)

  for (i in 1:length(TOPLOT)){
    cord.x <- c(TOPLOT[[i]]$x)
    cord.y <- c(TOPLOT[[i]]$y)
    rgbc<-col2rgb(COLS[i])
    currCol<-rgb(rgbc[1],rgbc[2],rgbc[3],alpha = 100,maxColorValue = 255)
    polygon(cord.x,cord.y,col=currCol,border = NA)
    lines(TOPLOT[[i]],col=COLS[i],lwd=3)
    if (i == 1){
      legend('topleft',legend = LEGentries,col=COLS,lwd=3,bty = 'n')
    }
  }
}
distPlot<-function(dataset,title,XLIM,YLIMS){
  clnames<-str_split(colnames(dataset),'---')
  clnames<-unlist(lapply(clnames,function(x){x[1]}))
  
  site<-str_split(colnames(dataset),'---')
  site<-unlist(lapply(site,function(x){x[[2]]}))
  
  
  site1dataset<-dataset[,which(site=='Broad')]
  site2dataset<-dataset[,which(site=='Sanger')]
  
  cdist1<-c(as.dist(cor(site1dataset)))
  cdist2<-c(as.dist(cor(site2dataset)))
  
  cdist<-cor(dataset)
  
  ucl<-unique(clnames)
  
  cdistSAME<-NULL
  
  for (i in 1:length(ucl)){
    id<-which(clnames==ucl[i])
    
    cdistSAME<-c(cdistSAME,cdist[id[1],id[2]])
    cdist[id[1],id[2]]<-NA
    cdist[id[2],id[1]]<-NA
    
    
  }
  
  cdistALL<-as.dist(cdist)
  
  multDensPlot(list(density(cdist1),
                    density(cdist2),
                    density(cdistALL,na.rm = TRUE),
                    density(cdistSAME,na.rm =TRUE)),XLIMS=XLIM,TITLE = title,YLIMS=YLIMS,
               COLS = c('purple','blue','gray','darkgreen'),LEGentries = c('Broad',
                                                                           'Sanger',
                                                                           'Overall',
                                                                           'Same Cell Line'))
  
  print(paste('Expected R within screen (Broad)',median(cdist1)))
  print(paste('Expected R within screen (Sanger)',median(cdist2)))
  print(paste('Expected overall R',median(cdistALL,na.rm = TRUE)))
  print(paste('R withn cell line across screens',median(cdistSAME)))
  
  print(t.test(cdistSAME,cdist1)$p.value)
  print(t.test(cdistSAME,cdist2)$p.value)
  
}

dataQuality<-function(dataset){
  
  data(BAGEL_essential)
  data(BAGEL_nonEssential)
  
  nsamples<-ncol(dataset)
  
  qc<-vector()
  for (i in 1:nsamples){
    FCs<-dataset[,i]
    names(FCs)<-rownames(dataset)
    RES<-PrRc_Curve(FCs,positives = BAGEL_essential,negatives = BAGEL_nonEssential,display = FALSE,FDRth = 0.05)
    qc[i]<-RES$Recall
  }
  names(qc)<-colnames(dataset)
  return(qc)
}
DisCar<-function(dataset,ref_Essential,ref_nonEssential,CL,th=0.05,SubSet=NULL){
  id<-grep(CL,colnames(dataset))

  FC_broad<-dataset[,id[1]]
  names(FC_broad)<-rownames(dataset)
  RES<-ccr.PrRc_Curve(FC_broad,ref_Essential,ref_nonEssential,FDRth = th,display = FALSE)
  DEPgenes_Broad<-names(which(FC_broad<RES$sigthreshold))

  if(length(SubSet)>0){
    DEPgenes_Broad<-intersect(DEPgenes_Broad,SubSet)
  }

  FC_sanger<-dataset[,id[2]]
  names(FC_sanger)<-rownames(dataset)
  RES<-ccr.PrRc_Curve(FC_sanger,ref_Essential,ref_nonEssential,FDRth = th,display = FALSE)
  DEPgenes_Sanger<-names(which(FC_sanger<RES$sigthreshold))

  if(length(SubSet)>0){
    DEPgenes_Sanger<-intersect(DEPgenes_Sanger,SubSet)
  }

  BroadOnly<-setdiff(DEPgenes_Broad,DEPgenes_Sanger)
  SangerOnly<-setdiff(DEPgenes_Sanger,DEPgenes_Broad)

  singleTermEnrBroad<-function(i){
    inC<-intersect(BroadOnly,GOterms$GS[[i]])
    N<-nrow(dataset)
    n<-length(GOterms$GS[[i]])
    k<-length(BroadOnly)
    x<-length(inC)

    p<-my.hypTest(x,k,n,N)
    return(c(i,p,paste(inC,collapse=', ')))
  }
  singleTermEnrSanger<-function(i){
    inC<-intersect(SangerOnly,GOterms$GS[[i]])
    N<-nrow(dataset)
    n<-length(GOterms$GS[[i]])
    k<-length(BroadOnly)
    x<-length(inC)

    p<-my.hypTest(x,k,n,N)
    return(c(i,p,paste(inC,collapse=', ')))
  }


  RES_Broad<-do.call(rbind,lapply(1:length(GOterms$NAMES),function(x) singleTermEnrBroad(x)))
  RES_Sanger<-do.call(rbind,lapply(1:length(GOterms$NAMES),function(x) singleTermEnrSanger(x)))

  RES_Broad<-data.frame(Id=as.numeric(RES_Broad[,1]),
                        GOterm=GOterms$NAMES[as.numeric(RES_Broad[,1])],
                        p = as.numeric(RES_Broad[,2]),
                        Genes = RES_Broad[,3],
                        stringsAsFactors = FALSE)

  sl<-str_length(RES_Broad$Genes)
  RES_Broad<-RES_Broad[sl>0,]

  RES_Broad<-cbind(RES_Broad,Adj.p=p.adjust(RES_Broad$p))
  RES_Broad<-RES_Broad[order(RES_Broad$Adj.p),]


  RES_Sanger<-data.frame(Id=as.numeric(RES_Sanger[,1]),
                         GOterm=GOterms$NAMES[as.numeric(RES_Sanger[,1])],
                         p = as.numeric(RES_Sanger[,2]),
                         Genes = RES_Sanger[,3],
                         stringsAsFactors = FALSE)

  sl<-str_length(RES_Sanger$Genes)
  RES_Sanger<-RES_Sanger[sl>0,]

  RES_Sanger<-cbind(RES_Sanger,Adj.p=p.adjust(RES_Sanger$p))
  RES_Sanger<-RES_Sanger[order(RES_Sanger$Adj.p),]


  return(list(RES_Sanger=RES_Sanger,RES_Broad=RES_Broad))
}
AllDisagreementCar<-function(dataset,CLs,GOterms,ref_Essential,ref_nonEssential,
                             ncbroad=39,ncsanger=20,
                             filePrefix1,filePrefix2){

  ncls<-length(CLs)

  FoundGenes<-lapply(1:length(GOterms$NAMES),function(x) '')

  EnrichmentsBroad<-list()
  EnrichmentsSanger<-list()

  for (i in 1:ncls){
    print(paste('Characterising the study-exclusive dependencies found in cell line: ',
          CLs[i],' (',i,' out of ',ncls,')',sep=''))
    RES<-DisCar(dataset = dataset,CL = CLs[i],ref_Essential,ref_nonEssential)

    tmpBroad<-RES$RES_Broad$Adj.p[RES$RES_Broad$Adj.p<0.05]
    names(tmpBroad)<-RES$RES_Broad$GOterm[RES$RES_Broad$Adj.p<0.05]
    EnrichmentsBroad[[i]]<-tmpBroad

    tmpSanger<-RES$RES_Sanger$Adj.p[RES$RES_Sanger$Adj.p<0.05]
    names(tmpSanger)<-RES$RES_Sanger$GOterm[RES$RES_Sanger$Adj.p<0.05]

    EnrichmentsSanger[[i]]<-tmpSanger
  }

  toSave<-qpcR:::cbind.na(names(EnrichmentsBroad[[1]]),EnrichmentsBroad[[1]])
  for(i in 2:length(CLs)) toSave <- qpcR:::cbind.na(toSave,
                                                    qpcR:::cbind.na(names(EnrichmentsBroad[[i]]),EnrichmentsBroad[[i]]))
  rownames(toSave)<-NULL
  colnames(toSave)<-c(rbind(CLs,rep('',length(CLs))))
  write.table(toSave,file=paste(filePrefix1,'.txt'),sep='\t',quote=FALSE,row.names = FALSE)


  toSave<-qpcR:::cbind.na(names(EnrichmentsSanger[[1]]),EnrichmentsSanger[[1]])
  for(i in 2:length(CLs)) toSave <- qpcR:::cbind.na(toSave,
                                                    qpcR:::cbind.na(names(EnrichmentsSanger[[i]]),EnrichmentsSanger[[i]]))
  rownames(toSave)<-NULL
  colnames(toSave)<-c(rbind(CLs,rep('',length(CLs))))
  write.table(toSave,file=paste(filePrefix2,'.txt'),sep='\t',quote=FALSE,row.names = FALSE)


  enrichmentMapBroad<-do.call(cbind,lapply(EnrichmentsBroad,function(x) {
    sigPattenrn<-rep(1,length(GOterms$NAMES))
    ii<-which(is.element(GOterms$NAMES,names(x)))
    sigPattenrn[ii]<-x
    return(sigPattenrn)
  }))
  rownames(enrichmentMapBroad)<-GOterms$NAMES
  colnames(enrichmentMapBroad)<-CLs

  enrichmentMapSanger<-do.call(cbind,lapply(EnrichmentsSanger,function(x) {
    sigPattenrn<-rep(1,length(GOterms$NAMES))
    ii<-which(is.element(GOterms$NAMES,names(x)))
    sigPattenrn[ii]<-x
    return(sigPattenrn)
  }))
  rownames(enrichmentMapSanger)<-GOterms$NAMES
  colnames(enrichmentMapSanger)<-CLs

  enrichmentMapBroad<- -log10(enrichmentMapBroad)
  enrichmentMapSanger<- -log10(enrichmentMapSanger)

  enrichmentMapBroad<-enrichmentMapBroad[which(rowSums(enrichmentMapBroad>0)>0),]
  enrichmentMapSanger<-enrichmentMapSanger[which(rowSums(enrichmentMapSanger>0)>0),]

  enrichmentMapBroad<-enrichmentMapBroad[order(apply(enrichmentMapBroad,MARGIN = 1,'median'),decreasing=TRUE),]
  enrichmentMapSanger<-enrichmentMapSanger[order(apply(enrichmentMapSanger,MARGIN = 1,'median'),decreasing=TRUE),]

  pdf(paste(filePrefix1,'.pdf',sep=''),20,10)
  par(mfrow=c(1,2))
  par(mar=c(5,25,4,1))
  toPlot<-enrichmentMapBroad[which(rowSums(enrichmentMapBroad>0)>ncbroad),]
  toBarplot<-rowSums(toPlot>0)
  broadToBarplot<-toBarplot
  O<-order(toBarplot)
  names(toBarplot)<-str_replace_all(names(toBarplot),'_',' ')
  names(toBarplot)<-str_replace(names(toBarplot),' ',':')
  barplot((100*toBarplot[O]/79),horiz = TRUE,las=2,cex.names = 0.8,xlim=c(0,100),col='#80479C',border = FALSE,
          main=paste(length(O), 'GO:Terms enriched in more than ',ncbroad,'\nBroad-exclusive essential gene sets (Adj.p < 0.05)'),
          xlab='% cell lines')
  par(mar=c(5,0,4,4))
  boxplot((t(toPlot[O,])),las=2,horizontal=TRUE,cex.axis=1,names=rep('',length(O)),frame.plot=FALSE,
          xlab='-log10 (enrichment p)')
  dev.off()

  pdf(paste(filePrefix2,'.pdf',sep=''),20,2.623)
  par(mfrow=c(1,2))
  par(mar=c(5,25,4,1))
  toPlot<-enrichmentMapSanger[which(100*rowSums(enrichmentMapSanger>0)/79>ncsanger),]
  nEnrichments<-length(which(rowSums(enrichmentMapSanger>0)>0))
  toBarplot<-rowSums(toPlot>0)
  O<-order(toBarplot)
  names(toBarplot)<-str_replace_all(names(toBarplot),'_',' ')
  names(toBarplot)<-str_replace(names(toBarplot),' ',':')
  barplot((100*toBarplot[O]/79),horiz = TRUE,las=2,cex.names = 0.8,xlim=c(0,100),col='#3A53A4',border = FALSE,
          main=paste(length(O), 'GO:Terms enriched in more than ',ncsanger,'\nSanger-exclusive essential gene sets (Adj.p < 0.05)'),
          xlab='% cell lines')

  par(mar=c(5,0,4,4))
  boxplot((t(toPlot[O,])),las=2,horizontal=TRUE,cex.axis=1,names=rep('',length(O)),frame.plot=FALSE,
          xlab='-log10 (enrichment p)')
  dev.off()



  return(list(EMB=enrichmentMapBroad,EMS=enrichmentMapSanger))

}

PR<-function (loc_RES_S, loc_RES_B, es_th = 1, fdr_th = 0.01, study1name = "study1", 
              study2name = "study2") 
{
  WW <- loc_RES_S
  WA <- loc_RES_B
  ii <- which(WW$fdr_oa < fdr_th)
  O <- ii[order(WW$p[ii])]
  prctl <- floor(quantile(1:length(O), c(0.2, 0.4, 0.6, 0.8, 
                                         1)))
  cols <- makeTransparent(colorRampPalette(c("purple", "orange"))(5), 
                          120)
  Npos <- vector()
  AUPREC <- vector()
  rndAUPREC <- vector()
  AUROC <- vector()
  rndAUROC <- rep(0.5, 5)
  for (i in 1:length(prctl)) {
    positives <- O[1:prctl[i]]
    Npos[i] <- length(positives)
    observed <- rep(0, nrow(WA))
    names(observed) <- paste("a", 1:nrow(WA))
    observed[positives] <- 1
    predictions <- -log10(WA$p)
    names(predictions) <- paste("a", 1:nrow(WA))
    RES <- ccr.PrRc_Curve(FCsprofile = -predictions, positives = names(which(observed > 
                                                                               0)), negatives = names(which(observed == 0)), display = FALSE)
    AUPREC[i] <- RES$AUC
    rndAUPREC[i] <- RES$RND
    if (i == 1) {
      plot(RES$curve[, "recall"], RES$curve[, "precision"], 
           type = "l", lwd = 4, xlab = "Recall", ylab = "Precision", 
           col = cols[i], main = "",cex.axis=1,cex.lab=1,tcl=0.5,tck=-0.01)
    }
    else {
      lines(RES$curve[, "recall"], RES$curve[, "precision"], 
            col = cols[i], lwd = 4)
    }
    #abline(h = min(RES$curve[, "precision"]), col = cols[i], 
    #       lty = 1)
  }
  #legend("left", "random", lty = 1, bty = "n")
  
  
  
}

ROC<-function (loc_RES_S, loc_RES_B, es_th = 1, fdr_th = 0.01, study1name = "study1", 
               study2name = "study2") {
  WW <- loc_RES_S
  WA <- loc_RES_B
  ii <- which(WW$fdr_oa < fdr_th)
  O <- ii[order(WW$p[ii])]
  prctl <- floor(quantile(1:length(O), c(0.2, 0.4, 0.6, 0.8, 
                                         1)))
  cols <- makeTransparent(colorRampPalette(c("purple", "orange"))(5), 
                          120)
  Npos <- vector()
  AUPREC <- vector()
  rndAUPREC <- vector()
  AUROC <- vector()
  rndAUROC <- rep(0.5, 5)
  for (i in 1:length(prctl)) {
    positives <- O[1:prctl[i]]
    observed <- rep(0, nrow(WA))
    names(observed) <- paste("a", 1:nrow(WA))
    observed[positives] <- 1
    predictions <- -log10(WA$p)
    names(predictions) <- paste("a", 1:nrow(WA))
    RES <- ROC_curve(-predictions, positives = names(which(observed > 
                                                             0)), negatives = names(which(observed == 0)))
    AUROC[i] <- RES$auc
    if (i == 1) {
      plot(RES$specificities, RES$sensitivities, type = "l", 
           lwd = 4, xlab = "Specificity", ylab = "Recall", 
           col = cols[i], xlim = c(1, 0),cex.axis=1,cex.lab=1,tcl=0.5,tck=-0.01,cex.main=0.9)
    }
    else {
      lines(RES$specificities, RES$sensitivities, col = cols[i], 
            lwd = 4)
    }
  }
  lines(c(0, 1), c(1, 0), col = "gray")
  legend("bottomright", c("20%", "40%", "60%", "80%", "100%"), 
         title = "Significance quantile", col = cols, lwd = 4)
  
}

plotAssociationExamplesCP<-function (gene, feat, loc_CombinedDataset, loc_site, loc_cids) 
{
  essS <- loc_CombinedDataset[gene, loc_site == "Sanger"]
  essB <- loc_CombinedDataset[gene, loc_site == "Broad"]
  pattern <- MoBEM[feat, as.character(loc_cids[loc_site == 
                                                 "Sanger"])]
  cols <- rep("darkgray", length(essB))
  cols[which(pattern == 1)] <- "darkgreen"
  beeswarm(essB ~ pattern, corral = "wrap", bg = c(makeTransparent("gray"), 
                                                   makeTransparent("darkgreen")), pch = 21, col = c("gray", 
                                                                                                    "darkgreen"), cex = 1.5, ylim = range(essB), las = 2, 
           labels = c("absent", "present"), xlab = feat, ylab = "",axes=F)
  par(new = TRUE)
  boxplot(essB ~ pattern, col = NA, ylim = range(essB), frame.plot = FALSE, 
          xaxt = "n",  outline = FALSE,tcl=0.5,tck=-0.01)
  beeswarm(essS ~ pattern, corral = "wrap", bg = c(makeTransparent("gray"), 
                                                   makeTransparent("darkgreen")), pch = 21, col = c("gray", 
                                                                                                    "darkgreen"), cex = 1.5, ylim = range(essS), las = 2, 
           labels = c("absent", "present"), xlab = feat, ylab = "",axes=F)
  par(new = TRUE)
  boxplot(essS ~ pattern, col = NA, ylim = range(essS), frame.plot = FALSE, 
          xaxt = "n",yaxt="n", outline = FALSE,tcl=0.5,tck=-0.01)
}


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

###Read in all data ###
GSEAfunction<-function(rankedlist,geneset){
  guse<-intersect(rankedlist,geneset)
  Inset<-(rankedlist%in%guse)*1
  sizeset<-length(guse)
  N_r<-sizeset
  PosScores<-Inset/N_r
  NegScores<-(!rankedlist%in%guse)*1
  N_Nh<-length(rankedlist)-sizeset
  NegScores<-NegScores/N_Nh
  Allscores<-PosScores-NegScores
  RunningSum<-cumsum(Allscores)
  maxdev<-RunningSum[which.max(abs(RunningSum))]
  return(list(ESscore=maxdev,RunningSum=RunningSum))
}
GetSigProfiles<-function(dataset,CLs,ref_Essential,ref_nonEssential,thresh=0.05){
  ncls <- length(CLs)
  
  DepletedBroad <- list()
  DepletedSanger <- list()
  th<-thresh
  for (i in 1:ncls) {
    CL<-CLs[i]
    id <- grep(CL, colnames(dataset))
    FC_broad <- dataset[, id[1]]
    names(FC_broad) <- rownames(dataset)
    RES <- ccr.PrRc_Curve(FC_broad, ref_Essential, ref_nonEssential,
                          FDRth = th, display = FALSE)
    DepletedBroad[[i]] <- names(which(FC_broad < RES$sigthreshold))
    FC_sanger <- dataset[, id[2]]
    names(FC_sanger) <- rownames(dataset)
    RES <- ccr.PrRc_Curve(FC_sanger, ref_Essential, ref_nonEssential,
                          FDRth = th, display = FALSE)
    DepletedSanger[[i]] <- names(which(FC_sanger < RES$sigthreshold))
    
  }
  return(list(DepletedBroad<-unique(unlist(DepletedBroad)),DepletedSanger<-unique(unlist(DepletedSanger))))
}

classPerf<-function(dataset,qualityTH=Inf,QC=NULL){
  
  
  if(length(QC)>0){
    clnames<-str_split(colnames(dataset),'---')
    clnames<-unlist(lapply(clnames,function(x){x[1]}))
    
    uc<-unique(clnames)
    for (i in 1:length(uc)){
      QC[which(clnames==uc[i])]<-rep(min(QC[which(clnames==uc[i])]),2)
    }
    
    dataset<-dataset[,QC>=qualityTH]
  }
  
  clnames<-str_split(colnames(dataset),'---')
  clnames<-unlist(lapply(clnames,function(x){x[1]}))
  
  site<-str_split(colnames(dataset),'---')
  site<-unlist(lapply(site,function(x){x[[2]]}))
  
  cdist<-as.matrix(1-cor(dataset))
  
  #cdist<-as.matrix(dist(t(dataset),method = 'manhattan'))
  
  ncls<-ncol(cdist)
  
  MATCH<-NULL
  for (i in 1:ncls){
    neighbours<-clnames[setdiff(order(cdist[i,]),i)]
    
    MATCH<-rbind(MATCH,neighbours==clnames[i])
  }
  
  MATCH<-MATCH+0
  rownames(MATCH)<-rownames(cdist)
  curv<-cumsum(colSums(MATCH))/ncol(cdist)
  
  return(list(MATCHmat=MATCH,CURV=curv))
}

RNAbiomarkers<-function(CombineRNA,BroadDep,SangerDep){
  #number of genes in Sanger RNA expression file:
  N<-21669
  
  
  CLid<-sapply(colnames(CombineRNA),function(x) strsplit(x,"...",fixed=TRUE)[[1]][2])
  Broad_Var<-CombineRNA[,names(CLid)[CLid=="Broad"]]
  Sanger_Var<-CombineRNA[,names(CLid)[CLid=="Sanger"]]
  colnames(Broad_Var)<-sapply(colnames(Broad_Var),function(x) strsplit(x,"...",fixed=TRUE)[[1]][1])
  colnames(Sanger_Var)<-sapply(colnames(Sanger_Var),function(x) strsplit(x,"...",fixed=TRUE)[[1]][1])
  #can do this the long way round:
  fullmatrixBroad<-rbind(Broad_Var,BroadDep[,colnames(Broad_Var)])
  
  fullmatrixSanger<-rbind(Sanger_Var,SangerDep[,colnames(Sanger_Var)])
  
  #now do the correlation matrices, although only want a subset of all these correlations!
  
  allCorBroad<-cor(t(fullmatrixBroad))
  allCorSanger<-cor(t(fullmatrixSanger))
  
  ExtractCorBroad<-allCorBroad[,(nrow(Broad_Var)+1):ncol(allCorBroad)]
  ExtractCorBroad<-ExtractCorBroad[1:nrow(Broad_Var),]
  
  BroadCorInput<-as.data.frame(ExtractCorBroad)
  BroadCorInput$geneRNA<-rownames(BroadCorInput)
  BroadCorInput<-melt(BroadCorInput)
  
  ExtractCorSanger<-allCorSanger[,(nrow(Sanger_Var)+1):ncol(allCorSanger)]
  ExtractCorSanger<-ExtractCorSanger[1:nrow(Sanger_Var),]
  
  OverallCor<-cor(as.vector(ExtractCorBroad),as.vector(ExtractCorSanger))
  
  print(paste("Overall correlation between correlations!",OverallCor))
  
  SangerCorInput<-as.data.frame(ExtractCorSanger)
  SangerCorInput$geneRNA<-rownames(SangerCorInput)
  SangerCorInput<-melt(SangerCorInput)
  
  
  colnames(Broad_Var)<-paste0(colnames(Broad_Var),"---Broad")
  colnames(Sanger_Var)<-paste0(colnames(Sanger_Var),"---Sanger")
  BothRNA<-cbind(Broad_Var,Sanger_Var)
  
  ## want to do the correlation tests between the gene expression vs dependency scores
  ##systematic correlation tests on ExtractCorBroad and ExtractCorSanger
  r2Broad<-ExtractCorBroad*ExtractCorBroad
  r2Sanger<-ExtractCorSanger*ExtractCorSanger
  Broad_denom<-sqrt((1-r2Broad)/(ncol(Broad_Var)-2))
  Sanger_denom<-sqrt((1-r2Sanger)/(ncol(Sanger_Var)-2))
  
  tStat_Broad<-ExtractCorBroad/Broad_denom
  tStat_Sanger<-ExtractCorSanger/Sanger_denom
  
  pval_Broad<-2*pt(abs(tStat_Broad),ncol(Broad_Var)-2,lower.tail=FALSE)
  pval_Sanger<-2*pt(abs(tStat_Sanger),ncol(Sanger_Var)-2,lower.tail=FALSE)
  
  inputtest<-as.vector(pval_Broad)
  fdr_Broad<-qvalue(as.vector(pval_Broad))$qvalues
  fdrthresh_Broad<-qvalue(as.vector(pval_Broad),fdr.level = 0.05)$significant
  fdr_Sanger<-qvalue(as.vector(pval_Sanger))$qvalues
  fdrthresh_Sanger<-qvalue(as.vector(pval_Sanger),fdr.level=0.05)$significant
  
  fdrMat_Broad<-matrix(fdr_Broad,nrow=nrow(pval_Broad),ncol=ncol(pval_Broad),byrow=FALSE)
  fdrMat_Sanger<-matrix(fdr_Sanger,nrow=nrow(pval_Sanger),ncol=ncol(pval_Sanger),byrow=FALSE)
  dimnames(fdrMat_Broad)<-dimnames(pval_Broad)
  dimnames(fdrMat_Sanger)<-dimnames(pval_Sanger)
  
  fdrMatT_Broad<-matrix(fdrthresh_Broad,nrow=nrow(pval_Broad),ncol=ncol(pval_Broad),byrow=FALSE)
  fdrMatT_Sanger<-matrix(fdrthresh_Sanger,nrow=nrow(pval_Sanger),ncol=ncol(pval_Sanger),byrow=FALSE)
  dimnames(fdrMatT_Broad)<-dimnames(pval_Broad)
  dimnames(fdrMatT_Sanger)<-dimnames(pval_Sanger)
  
  #from threshold can work out the minimum correlation value required to get 5% fdr.
  minT_B<-min(abs(ExtractCorBroad[fdrMatT_Broad]))
  minT_S<-min(abs(ExtractCorSanger[fdrMatT_Sanger]))
  
  RNA_data<-BroadCorInput
  colnames(RNA_data)<-c("gene","SSD","Broad")
  RNA_data$Sanger<-SangerCorInput[,"value"]
  
  RES_cP<-classPerf(BothRNA)
  
  ##Fishers exact test for 2x2 contingency, signif S, B, not signif S, B
  SignifBoth<-sum(fdrMatT_Broad*fdrMatT_Sanger)
  SignifSanger<-sum(fdrMatT_Sanger)-SignifBoth
  SignifBroad<-sum(fdrMatT_Broad)-SignifBoth
  NotSignif<-sum((fdrMatT_Broad+fdrMatT_Sanger)==0)
  
  FTMatrix<-matrix(c(NotSignif,SignifSanger,SignifBroad,SignifBoth),nrow=2)
  FTtest<-fisher.test(FTMatrix)
  
  FTtest<-my.hypTest(SignifBoth,SignifSanger+SignifBoth,SignifBroad+SignifBoth,N)
  
  return(list(BroadData=Broad_Var,SangerData=Sanger_Var,PerfAnalysis=RES_cP,CombineRNA=BothRNA,PlotRNAdata=RNA_data,B_Sthresh=c(minT_B,minT_S),FTTest=FTtest))
}


