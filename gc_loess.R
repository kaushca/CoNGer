options <- commandArgs(trailingOnly = T)
gc = options[1]
doc=options[2]
outF = options[3]

gcF = read.delim(gc,header=F)
df = read.delim(doc,check.names=F) #normalized coverage file
samples<-names(df)[5:ncol(df)]
print(cat("GC normalization of",samples,sep=" "))
if(ncol(df)>5){
  t<-rowSums(df[,5:ncol(df)])==0
} else{
  t<-(df[,5])==0
}
out<-df
eps<-0.001 #small amount to add to 0 coverage
for (i in samples){
  up<-quantile(df[,i],prob=0.9,names=F)
  f<-(df[,i]>0&df[,i]<up) & (gcF[,7]>20 & gcF[,7]<80)
  gc.loess<-loess(log2(df[f,i])~gcF[f,7],control=loess.control(surface="direct"))
  print(i)
  x<-predict(gc.loess,gcF[!t,7])
  out[!t,i]<-2^(log2(df[!t,i])-x)
  out[t,i]<-rep(0,sum(t))
}

write.table(out,outF,quote=F,sep="\t",row.names=F)