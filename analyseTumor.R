options <- commandArgs(trailingOnly = T)
tumor = options[1]
control = options[2]
outF = options[3]
outF2=outF

#require(MASS,quietly=T)
require(DNAcopy,quietly=T)

createCor<-function(c){
  corMat<-matrix(-1.2,ncol=ncol(c)-4,nrow=ncol(c)-4)
  
  for(i in 1:(ncol(c)-4)){
    for(j in 1:(ncol(c)-4)){
      if(i!=j){
        corMat[i,j]<-cor(c[,i+4],c[,j+4])
      }
    }
  }
  return(corMat)
}

selectControl<-function(t,c,thresh=0.90){
  cMat<-createCor(c)
  corT<-rep(-1.2,ncol(cMat))
  for(i in 5:ncol(c)){
    corT[i-4]<-cor(t,c[,i])
  }
  id<-which.max(corT)
  ids<-which(cMat[id,]>thresh)
  if(length(ids)<2){
    ids<-which(cMat[id,]>median(cMat[id,-id]))
  }
  ids<-append(ids,id)
  return(ids)
}

rowMads<-function(x){
  m<-rep(0,nrow(x))
  for (i in c(1:nrow(x))){
    m[i]<-mad(as.numeric(x[i,]))
  }
  return (m)
}

rowMedians<-function(x){
  m<-rep(0,nrow(x))
  mads<-rep(0,nrow(x))
  for (i in c(1:nrow(x))){
    m[i]<-median(as.numeric(x[i,]))
    mads[i]<-mad(as.numeric(x[i,]))
  }
  dt<-NULL
  dt$med<-m
  dt$mad<-mads
  dt<-data.frame(dt)
  return (dt)
}

madWins <- function(x,tau,k){
  xhat <- medianFilter(x,k)
  d <- x-xhat
  SD <- mad(d)
  z <- tau*SD
  xwin <- xhat + psi(d, z)
  outliers <- rep(0, length(x))
  outliers[x > xwin] <- 1
  outliers[x < xwin] <- -1
  return(list(ywin=xwin,sdev=SD,outliers=outliers))
}

psi <- function(x,z){
  xwin <- x
  xwin[x < -z] <- -z
  xwin[x > z] <- z
  return(xwin)
}


#########################################################################
# Function to calculate running median for a given a window size
#########################################################################

##Input:
### x: vector of numeric values
### k: window size to be used for the sliding window (actually half-window size)

## Output:
### runMedian : the running median corresponding to each observation

##Required by:
### getMad
### medianFilter


##Requires:
### none

medianFilter <- function(x,k){
  n <- length(x)
  filtWidth <- 2*k + 1
  
  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }
  
  runMedian <- runmed(x,k=filtWidth,endrule="median")
  
  return(runMedian)
  
}

createControl<-function(c,ids,dev=2){
  ids<-ids+4
  #print(ids)
  cc<-c[,ids]
  rM<-rowMedians(cc)
  #rMads<-rowMads(cc)
  #rMedians<-rowMedians(cc)
  rMads<-rM$mad
  rMedians<-rM$med
  upT<-rMedians+dev*rMads
  lowT<-rMedians-dev*rMads
  cont<-NULL
  DOC<-rep(NA,nrow(c))
  for(i in 1:nrow(c)){
    f<-!(cc[i,]>upT[i] & cc[i,]<lowT[i])
    if(sum(cc[i,])==0){
      DOC[i]=0
    }
    else{
      DOC[i]<-rowMeans(cc[i,f])
    }
  }
  cont$DOC<-DOC
  ratio<-log2(cc)-log2(DOC)
  cont$MAD<-rowMads(ratio)
  #resi<-cc
  #for (i in 1:ncol(cc)){
   # y<-rlm(cc[,i]~cont$DOC)
   # resi[,i]<-y$residuals
  #}
  #cont$resiMads<-rowMads(resi)
  #cont$resiMean<-rowMeans(resi)
  #cont$resiMedian<-rowMedians(resi)
  cont<-data.frame(cont)
  return(cont)
}

smoothData<-function(x){
  if(is.vector(x)){
    x<-madWins(x,2.5,5)$ywin
  } else{
    for(i in 1:ncol(x)){
      x[,i]<-madWins(x[,i],2.5,5)$ywin
    }
  }
  return(x)
}

segProfile<-function(ratios,sampleName,chr,loc){
  CNA.object<-CNA(genomdat=as.vector(ratios),chrom=chr,maploc=loc,
                  data.type="logratio",sampleid=sampleName)
  smoothed.CNA.object<-smooth.CNA(CNA.object)
  segments<-segment(smoothed.CNA.object,verbose=1,alpha=0.01)
  #plot(segment.smoothed.CNA.object,plot.type="chrombysample",xmaploc=TRUE)
  
  return(segments)
}

control<-read.delim(control,check.names=F)
tumor<-read.delim(tumor,check.names=F)

print("Data smoothing...")
control[,5:ncol(control)]<-smoothData(control[,5:ncol(control)])
tumor[,5:ncol(tumor)]<-smoothData(tumor[,5:ncol(tumor)])

thresh<-0.15 # threshold on control coverage

for(i in 5:ncol(tumor)){
#for(i in 8:ncol(tumor)){
  print(paste("Predicting CNV in sample",names(tumor)[i],sep=" "))
  ids<-selectControl(tumor[,i],control)
  cont<-createControl(control,ids)
  #y<-rlm(tumor[,i]~cont$DOC)
  #resi<-y$residuals
  resi<-log2(tumor[,i])-log2(cont$DOC)
  f<-cont$DOC>thresh&cont$DOC<quantile(cont$DOC,prob=0.8)
  sdLOESS<-loess(cont$MAD[f]~log2(cont$DOC[f]),
                 control=loess.control(surface="direct"))
  print("SDs calculated")
  
  ###Finding starting points of regions
  changes<-c(1,which(control$Start[-1]!=control$End[-nrow(control)])+1)
  ends<-c(which(control$Start[-1]!=control$End[-nrow(control)]),nrow(control))
  nChanges<-NULL
  nEnds<-NULL
  pval<-NULL
  tumorDOC<-NULL
  controlDOC<-NULL
  ratios<-NULL
  pvalMean<-NULL
  
  mad.thresh<-quantile(cont$MAD,na.rm=T,prob=0.99)
  
  ###Perform t-test
  for (j in 1:length(changes)){
    f1<-cont$DOC[changes[j]:ends[j]]>thresh
    f2<-!is.infinite(resi[changes[j]:ends[j]])
    f3<-(cont$MAD[changes[j]:ends[j]]<mad.thresh)&!is.na(cont$MAD[changes[j]:ends[j]])
    filt<-f1&f2&f3
    
    #####include the doc dependent MADs instead of residual MADs
    ff<-resi[changes[j]:ends[j]][filt]
    #print(ff)
    if(length(ff)>=1){
      coverage<-(log2(cont$DOC[changes[j]:ends[j]][filt])+log2(tumor[changes[j]:ends[j],i][filt]))/2
      mads<-predict(sdLOESS,coverage)
      outlier<-which(mads>quantile(cont$MAD,prob=0.95,na.rm=T))
      mads[outlier]<-mean(cont$MAD,na.rm=T)
      t_doc<-mean(tumor[changes[j]:ends[j],i][filt])
      c_doc<-mean(cont$DOC[changes[j]:ends[j]][filt])
      tumorDOC<-append(tumorDOC,t_doc)
      controlDOC<-append(controlDOC,c_doc)
      ratios<-append(ratios,log2(t_doc)-log2(c_doc))
      nChanges<-append(nChanges,changes[j])
      nEnds<-append(nEnds,ends[j])
      if(length(ff)==1){
        x<-pnorm(ff,sd=mads)
        p<-2*ifelse(x>0.5,1-x,x)
        p2<-p
      }
      else{
        testValues<-ff/mads
        ##checking if values are constant
        if(all(abs(max(testValues)-min(testValues))<1e-10)){
          x<-pnorm(mean(testValues),sd=mean(mads))
          p<-2*ifelse(x>0.5,1-x,x)
          p2<-p
        }
        else{
          p<-t.test(testValues,mu=0)$p.value
          sd<-predict(sdLOESS,(log2(t_doc)+log2(c_doc))/2)
          p2<-pnorm(mean(ff),sd=sd/sqrt(length(ff)))
        }
      }
      pval<-append(pval,p)
      pvalMean<-append(pvalMean,p2)
    }
  }
  tName<-names(tumor)[i]
  outF<-paste(outF2,"/",tName,sep="")
  
  tOut<-NULL
  tOut$Chr<-tumor$Chr[nChanges]
  tOut$Start<-tumor$Start[nChanges]
  tOut$End<-tumor$End[nEnds]
  tOut$gene<-tumor$ID[nChanges]
  tOut$tumorDOC<-tumorDOC
  tOut$controlDOC<-controlDOC
  tOut$log2Ratio<-ratios
  tOut$pval<-pval
  pvalMean<-p.adjust(pvalMean,method="BH")
  tOut$pvalOnMean<-pvalMean
  
 # x.start<-which(tOut$Start>=38545687 & tOut$Chr=="17 " & tOut$Start<=41600867)
 # brca<-which(tOut$Start>=41196189 & tOut$Chr=="17 " & tOut$Start<=41275667)
#   png(filename=paste(outF,"/",tName,"brca1.png",sep=""),height=500,width=800,
#       bg="white")
#   plot(tOut$Start[x.start],tOut$log2Ratio[x.start],pch=20,cex=0.8,
#        ylim=c(min(-2,min(tOut$log2Ratio[x.start])),
#               max(2,max(tOut$log2Ratio[x.start]))),
#        col=(tOut$pval[x.start]<0.05&tOut$pvalOnMean[x.start]<0.05)+1,
#        main=paste(tName,"BRCA1",sep=" "),xlab="genome cordinates (bp)",ylab="log2 ratio")
#   abline(v=41196189,col="grey",lty="dashed")
#   abline(v=41276656,col="grey",lty="dashed")
#   text(x=41196189+200,y=max(2,max(tOut$log2Ratio[x.start]))-0.1,labels="BRCA1")
#   dev.off()
#   
#   x.start<-which(tOut$Start>=30782519 & tOut$Chr=="13 " & tOut$Start<=34540107)
#   #brca<-which(tOut$Start>=41196189 & tOut$Chr=="17 " & tOut$Start<=41275667)
#   png(filename=paste(outF,"/",tName,"brca2.png",sep=""),height=500,width=800,
#       bg="white")
#   plot(tOut$Start[x.start],tOut$log2Ratio[x.start],pch=20,cex=0.8,
#        ylim=c(min(-2,min(tOut$log2Ratio[x.start])),
#               max(2,max(tOut$log2Ratio[x.start]))),
#        col=(tOut$pval[x.start]<0.05&tOut$pvalOnMean[x.start]<0.05)+1,
#        main=paste(tName,"BRCA2",sep=" "),xlab="genome cordinates (bp)",ylab="log2 ratio")
#   abline(v=32889411,col="grey",lty="dashed")
#   abline(v=32973704,col="grey",lty="dashed")
#   text(x=32889411+200,y=max(2,max(tOut$log2Ratio[x.start]))-0.1,labels="BRCA2")
#   dev.off()
  
  tOut<-data.frame(tOut)
  
  write.table(tOut,paste(outF,"/",tName,"_cnv.txt",sep=""),quote=F,sep="\t",row.names=F)
  
  print("Filtering true positive calls...")
  f<-which(tOut$pvalOnMean<=0.01)
  sdEstimate<-sd(tOut$log2Ratio[-f])
  meanEst<-mean(tOut$log2Ratio[-f])
  x<-pnorm(tOut$log2Ratio[f],sd=sdEstimate)
  pval<-2*ifelse(x>0.5,1-x,x)
  #pval<-p.adjust(pval,method="BH")
  
  dev<-NULL
  for(j in f){
    x<-which(tOut$Chr==tOut$Chr[j] & tOut$Start>=tOut$Start[j]-10000 & tOut$End<=tOut$End[j]+10000 & tOut$pvalOnMean>0.01)
        
    if(length(x)<4){
      dataset<-tOut[(tOut$pvalOnMean>0.01&tOut$Chr==tOut$Chr[j]),]
      ii<-which(dataset$Start==tOut$Start[j])
      s<-max(1,ii-2)
      e<-min(ii+2,nrow(dataset))
      r<-dataset$log2Ratio[c(s:e)]
    }
    else{
      r<-tOut$log2Ratio[x]
    }
    localSD<-sd(r)
    deviation<-abs(tOut$log2Ratio[j]-mean(r))/localSD
    dev<-append(dev,deviation)
  }
 
  result<-NULL
  result$Chr<-tOut$Chr[f]
  result$Start<-tOut$Start[f]
  result$End<-tOut$End[f]
  result$geneID<-tOut$gene[f]
  result$ratio<-tOut$log2Ratio[f]
  result$pval<-pval
  result$localDeviation<-dev
  
  result<-data.frame(result)
  write.table(result,paste(outF,"/",tName,"_cnvFiltered.txt",sep=""),quote=F,sep="\t",row.names=F)
  
#   x.start<-which(tOut$Start>=38545687 & tOut$Chr=="17 " & tOut$Start<=41600867)
#   png(filename=paste(outF,"/",tName,"brca1.png",sep=""),height=500,width=800,
#       bg="white")
#   f1<-result$Start>=38545687 & result$Chr=="17 " & result$Start<=41600867
#   f2<-result$pval<=0.05
#   col.cord<-f[which(f1&f2)]
#   col.cord<-sapply(col.cord,function(x)which(x == x.start))
#   col.start<-rep(0,length(x.start))
#   if(length(col.cord)>0){
#     col.start[unlist(col.cord)]<-1
#   }
#   plot(tOut$Start[x.start],tOut$log2Ratio[x.start],pch=20,cex=0.8,
#        ylim=c(min(-2,min(tOut$log2Ratio[x.start])),
#               max(2,max(tOut$log2Ratio[x.start]))),col=(col.start)+1,
#        main=paste(tName,"BRCA1",sep=" "),xlab="genome cordinates (bp)",ylab="log2 ratio")
#   abline(v=41196189,col="grey",lty="dashed")
#   abline(v=41276656,col="grey",lty="dashed")
#   text(x=41196189+200,y=max(2,max(tOut$log2Ratio[x.start]))-0.1,labels="BRCA1")
#   dev.off()
#   
#   x.start<-which(tOut$Start>=23400115 & tOut$Chr=="16 " & tOut$Start<=27269004)
#   png(filename=paste(outF,"/",tName,"palb2.png",sep=""),height=500,width=800,
#       bg="white")
#   f1<-result$Start>=23400115 & result$Chr=="16 " & result$Start<=27269004
#   f2<-result$pval<=0.05
#   col.cord<-f[which(f1&f2)]
#   col.cord<-sapply(col.cord,function(x)which(x == x.start))
#   col.start<-rep(0,length(x.start))
#   if(length(col.cord)>0){
#     col.start[unlist(col.cord)]<-1
#   }
#   plot(tOut$Start[x.start],tOut$log2Ratio[x.start],pch=20,cex=0.8,
#        ylim=c(min(-2,min(tOut$log2Ratio[x.start])),
#               max(2,max(tOut$log2Ratio[x.start]))),col=(col.start)+1,
#        main=paste(tName,"PALB2",sep=" "),xlab="genome cordinates (bp)",ylab="log2 ratio")
#   abline(v=23603459,col="grey",lty="dashed")
#   abline(v=23641157,col="grey",lty="dashed")
#   text(x=23603459+200,y=max(2,max(tOut$log2Ratio[x.start]))-0.1,labels="PALB2")
#   dev.off()
  
#   x.start<-which(tOut$Start>=30782519 & tOut$Chr=="13 " & tOut$Start<=34540107)
#   png(filename=paste(outF,"/",tName,"brca2.png",sep=""),height=500,width=800,
#       bg="white")
#   f1<-result$Start>=30782519 & result$Chr=="13 " & result$Start<=34540107
#   f2<-result$pval<=0.05
#   col.cord<-f[which(f1&f2)]
#   col.cord<-sapply(col.cord,function(x)which(x == x.start))
#   col.start<-rep(0,length(x.start))
#   if(length(col.cord)>0){
#     col.start[unlist(col.cord)]<-1
#   }
#   plot(tOut$Start[x.start],tOut$log2Ratio[x.start],pch=20,cex=0.8,
#        ylim=c(min(-2,min(tOut$log2Ratio[x.start])),
#               max(2,max(tOut$log2Ratio[x.start]))),col=(col.start)+1,
#        main=paste(tName,"BRCA2",sep=" "),xlab="genome cordinates (bp)",ylab="log2 ratio")
#   abline(v=32889411,col="grey",lty="dashed")
#   abline(v=32973704,col="grey",lty="dashed")
#   text(x=32889411+200,y=max(2,max(tOut$log2Ratio[x.start]))-0.1,labels="BRCA2")
#   dev.off()
  
  dfChanges<-NULL
  dfChanges$Chr<-tumor$Chr
  dfChanges$loc<-tumor$Start
  dfChanges$cnv<-rep(2,nrow(tumor))
  dfChanges<-data.frame(dfChanges)
  colnames(dfChanges)<-c("Chr","loc",tName)
  
  ####segmentation
  segments<-segProfile(resi,tName,tumor$Chr,tumor$Start)
  kk<-kmeans(segments$output$seg.mean,centers=c(meanEst,log2(1.5),log2(0.5)))
  segVar<-abs(segments$output$seg.mean-meanEst)/sdEstimate
  fSeg<-segments$output[kk$cluster!=1|segVar>=2.5,]
  #fSeg<-segments$output[kk$cluster!=1&(segments$output$seg.mean<(-0.7)|segments$output$seg.mean>0.38),]
  #fSeg<-segments$output[segments$output$seg.mean>log2(1.3)]
  pval<-NULL
  ratio<-NULL
  for(k in 1:nrow(fSeg)){
    seg.start<-fSeg$loc.start[k]
    seg.end<-fSeg$loc.end[k]
    dat<-which(tumor$Chr==fSeg$chrom[k] & tumor$Start>=seg.start & tumor$Start<=seg.end)
    x<-resi[dat]
    c<-which(is.infinite(x) | is.na(x))
    if(length(c)>0){
      x<-x[-c]
    }
    n<-length(x)
    x.bar<-mean(x)
    h<-abs(min(x))-abs(max(x))
    if(n>1 & abs(h)>1e-6){
      p<-t.test(x,mu=0)$p.value
      #x<-pnorm(x.bar,sd=sdEstimate/sqrt(n),mean=meanEst)
      #p<-2*ifelse(x>0.5,1-x,x)
    }
    else{
      x<-pnorm(x.bar,sd=sdEstimate)
      p<-2*ifelse(x>0.5,1-x,x)
    }
    pval<-append(pval,p)
    ratio<-append(ratio,x.bar)
    if(p<=0.01){
      dfChanges[dat,tName]<-ifelse(x.bar<0,1,3)
    }
  }

  #pval<-p.adjust(pval,method="BH")
  segOut<-NULL
  segOut$Chr<-fSeg$chrom
  segOut$seg.Start<-fSeg$loc.start
  segOut$seg.End<-fSeg$loc.end
  segOut$log2Ratio<-ratio
  segOut$num.mark<-fSeg$num.mark
  segOut$pval<-pval
  
  segOut<-data.frame(segOut)
  #print(segOut)
  
#   png(filename=paste(outF,"/",tName,"chr17_seg.png",sep=""),height=500,width=800,
#       bg="white")
#   f<-segments$data$chrom=="17 " & !is.infinite(segments$data[,3]) & !is.na(segments$data[,3])
#   #plot(segments$data$maploc[f],segments$data[f,3],pch=20,cex=0.8,
#    #    main=paste(tName,"Chr17",sep=" "),xlab="genome cordinates (bp)",
#     #   ylab="log2 ratio",col=rgb(20,50,200,50,maxColorValue=255))
#   plot(segments$data[f,3],pch=20,cex=0.8,
#        main=paste(tName,"Chr17",sep=" "),xlab="genome loci",
#        ylab="log2 ratio",col=rgb(20,50,200,50,maxColorValue=255))
#   ff<-segments$output$chrom=="17 "
#   fff<-segOut$Chr=="17 " & segOut$pval<0.01
#   s=0
#   col2=rep("black",sum(ff))
#   lwd2=rep(1,sum(ff))
#   for(seg in 1:sum(ff)){
#     #lines(c(segments$output$loc.start[ff][seg],segments$output$loc.end[ff][seg]),
#      #      rep(segments$output$seg.mean[ff][seg],2))
#     c<-which(segments$output$loc.start[ff][seg]==segOut$seg.Start[fff])
#     if(length(c)>0){
#       col2[seg]<-"red"
#       lwd2[seg]<-2
#     }
#     e=s+segments$output$num.mark[ff][seg]
#     lines(seq(s+1,e),col=col2[seg],lwd=lwd2[seg],
#           rep(segments$output$seg.mean[ff][seg],segments$output$num.mark[ff][seg]))
#     s=e
#   }
#   d<-which(segments$data$maploc[f]<=41196189)
#   d<-d[length(d)]
#   abline(v=d,col="dark grey",lty="dashed")
#   d<-which(segments$data$maploc[f]<=41276656)
#   d<-d[length(d)]
#   abline(v=d,col="dark grey",lty="dashed")
#   dev.off()
#   
#   png(filename=paste(outF,"/",tName,"chr16_seg.png",sep=""),height=500,width=800,
#       bg="white")
#   f<-segments$data$chrom=="16 " & !is.infinite(segments$data[,3]) & !is.na(segments$data[,3])
#   #plot(segments$data$maploc[f],segments$data[f,3],pch=20,cex=0.8,
#   #    main=paste(tName,"Chr17",sep=" "),xlab="genome cordinates (bp)",
#   #   ylab="log2 ratio",col=rgb(20,50,200,50,maxColorValue=255))
#   plot(segments$data[f,3],pch=20,cex=0.8,
#        main=paste(tName,"Chr16",sep=" "),xlab="genome loci",
#        ylab="log2 ratio",col=rgb(20,50,200,50,maxColorValue=255))
#   ff<-segments$output$chrom=="16 "
#   fff<-segOut$Chr=="16 " & segOut$pval<0.01
#   s=0
#   col2=rep("black",sum(ff))
#   lwd2=rep(1,sum(ff))
#   for(seg in 1:sum(ff)){
#     #lines(c(segments$output$loc.start[ff][seg],segments$output$loc.end[ff][seg]),
#     #      rep(segments$output$seg.mean[ff][seg],2))
#     c<-which(segments$output$loc.start[ff][seg]==segOut$seg.Start[fff])
#     if(length(c)>0){
#       col2[seg]<-"red"
#       lwd2[seg]<-2
#     }
#     e=s+segments$output$num.mark[ff][seg]
#     lines(seq(s+1,e),col=col2[seg],lwd=lwd2[seg],
#           rep(segments$output$seg.mean[ff][seg],segments$output$num.mark[ff][seg]))
#     s=e
#   }
#   d<-which(segments$data$maploc[f]<=23603459)
#   d<-d[length(d)]
#   abline(v=d,col="dark grey",lty="dashed")
#   d<-which(segments$data$maploc[f]<=23641157)
#   d<-d[length(d)]
#   abline(v=d,col="dark grey",lty="dashed")
#   dev.off()
  
#   png(filename=paste(outF,"/",tName,"chr13_seg.png",sep=""),height=500,width=800,
#       bg="white")
#   f<-segments$data$chrom=="13 " & !is.infinite(segments$data[,3]) & !is.na(segments$data[,3])
#   #plot(segments$data$maploc[f],segments$data[f,3],pch=20,cex=0.8,
#   #    main=paste(tName,"Chr17",sep=" "),xlab="genome cordinates (bp)",
#   #   ylab="log2 ratio",col=rgb(20,50,200,50,maxColorValue=255))
#   plot(segments$data[f,3],pch=20,cex=0.8,
#        main=paste(tName,"Chr13",sep=" "),xlab="genome loci",
#        ylab="log2 ratio",col=rgb(20,50,200,50,maxColorValue=255))
#   ff<-segments$output$chrom=="13 "
#   fff<-segOut$Chr=="13 " & segOut$pval<0.01
#   s=0
#   col2=rep("black",sum(ff))
#   lwd2=rep(1,sum(ff))
#   for(seg in 1:sum(ff)){
#     #lines(c(segments$output$loc.start[ff][seg],segments$output$loc.end[ff][seg]),
#     #      rep(segments$output$seg.mean[ff][seg],2))
#     c<-which(segments$output$loc.start[ff][seg]==segOut$seg.Start[fff])
#     if(length(c)>0){
#       col2[seg]<-"red"
#       lwd2[seg]<-2
#     }
#     e=s+segments$output$num.mark[ff][seg]
#     lines(seq(s+1,e),col=col2[seg],lwd=lwd2[seg],
#           rep(segments$output$seg.mean[ff][seg],segments$output$num.mark[ff][seg]))
#     s=e
#   }
#   d<-which(segments$data$maploc[f]<=32889411)
#   d<-d[length(d)]
#   abline(v=d,col="dark grey",lty="dashed")
#   d<-which(segments$data$maploc[f]<=32973704)
#   d<-d[length(d)]
#   abline(v=d,col="dark grey",lty="dashed")
#   dev.off()
  
#   chr<-unique(tumor$Chr)
#   pdf(filename=paste(dir,"Profile.png",sep=""),height=750,width=1024,bg="white")
#   ch<-table(segments$data$chrom)
#   chrom<-unlist(names(ch)," ")[ch[]>0]
#   changes<-match(chrom,segments$data$chrom)
#   chrom<-chrom[order(changes)]
#   changes<-changes[order(changes)]
#   mt<-match(chr,chrom)
#   limits<-match(chrom,segments$data$chrom)
#   limits<-c(limits[-1]-1,nrow(segments$data))
#   limits<-limits[mt]
#   tot<-sum(as.numeric(segments$data$maploc[limits]))
#   end_limit<-cumsum(as.numeric(segments$data$maploc[limits]))
#   plot(c(0,tot),c(0,max(segments$data[,3])),type="n",xlab="",
#        ylab="DOC ratio",xaxt="n",main="CNV profile",ylim=c(-4,4))
#   end_limit2<-c(0,end_limit[-length(end_limit)])[match(chrom,chr)]
#   len<-rep(end_limit2,ch[chrom][])
#   x<-segments$data$maploc + len
#   points(x,segments$data[,3],pch=20,cex=0.8,col=rgb(0,0,255,50,maxColorValue=255))
#   abline(v=0,lty=1,col="lightgrey")
#   text(end_limit[1]/2,-3.5,chr[1], pos = 1, cex = 1)
#   for (j in 2:length(end_limit)) {
#     vpos = end_limit[j-1];
#     tpos = (end_limit[j-1]+end_limit[j])/2
#     text(tpos,-3.5,chr[j], pos = 1, cex = 1)
#     abline(v=vpos,lty=1,col="lightgrey")
#   }
  
  
  #dev.off()
  write.table(segOut,paste(outF,"/",tName,"_segmentsFiltered.txt",sep=""),quote=F,
              sep="\t",row.names=F)
  write.table(segOut[segOut$pval<=0.01,],paste(outF,"/",tName,"_significantSegments.txt",sep=""),quote=F,
              sep="\t",row.names=F)
  write.table(dfChanges,paste(outF,"/",tName,"_changes.txt",sep=""),quote=F,
              sep="\t",row.names=F)
  write.table(segments$output[,2:6],paste(outF,"/",tName,"_segments.txt",sep=""),
              quote=F,sep="\t",row.names=F)
  write.table(segments$data,paste(outF,"/",tName,"_rawData.txt",sep=""),quote=F,
              sep="\t",row.names=F)
}
