options<-commandArgs(trailingOnly=T)
segments_data=options[1]
segments_output=options[2]
segOut=options[3]
outF=options[4]
tName=options[5]

segments_data=read.delim(segments_data)
segments_output=read.delim(segments_output)
segOut=read.delim(segOut)
chr=unique(segments_data$chrom)
for (i in chr){
  png(filename=paste(outF,"/",tName,"_Chr",i,".png",sep=""),height=500,width=800,
      bg="white")
  f<-segments_data$chrom==i & !is.infinite(segments_data[,3]) & !is.na(segments_data[,3])
  plot(segments_data[f,3],pch=20,cex=0.8,
       main=paste(tName,"_Chr",i,sep=" "),xlab="genome loci",
       ylab="log2 ratio",col=rgb(20,50,200,50,maxColorValue=255),ylim=c(-4,4))
  ff<-segments_output$chrom==i
  fff<-(segOut$Chr==i | segOut$Chr==sub("\\s+$", "", i)) & segOut$pval<0.01
  s=0
  col2=rep("black",sum(ff))
  lwd2=rep(1,sum(ff))
  for(seg in 1:sum(ff)){
    c<-which(segments_output$loc.start[ff][seg]==segOut$seg.Start[fff])
    if(length(c)>0){
      col2[seg]<-"red"
      lwd2[seg]<-2
    }
    e=s+segments_output$num.mark[ff][seg]
    lines(seq(s+1,e),col=col2[seg],lwd=lwd2[seg],
          rep(segments_output$seg.mean[ff][seg],segments_output$num.mark[ff][seg]),
          ylim=c(-4,4))
    s=e
  }
  dev.off()
}

