### pairwise Fst between mng/sou and three non-native populations

## sliding window ngsFST
### mngsou vs non-native
### 50000 bp window; 10000 steps

library(reshape2)
rm(list=ls())
meta <- read.csv('data/meta-SNP.pops_edited.V2.csv')
meta$Site.Abb <- tolower(meta$Site.Abb)
all <- c()
### all output files (fst.txt)
a <- list.files('data/fst/',pattern="mngsou")
for(i in 1:length(a))
{
  tmp <- read.delim(paste('data/fst/',a[i],sep=""),header=T,row.names=NULL)
  tmp$pop <- substr(a[i],1,10)
  colnames(tmp) <- c("region","chr","midPos","Nsites","fst","pop")
  all <- rbind(all,tmp)
}
all$chr_pos <- paste(all$chr,"_",all$midPos,sep="")

pops <- c("mngsou.tmb","mngsou.osh","fdm.mngsou") #unique(all$pop); pops <- pops[c(3,2,1)]
chrs.cols <- data.frame(chrs=levels(all$chr),cols=c(rep(c("red","blue"),12)))

out <- c()
stats <- c()
pdf('output/Fstoutlier.mngsou.vs.non.pdf',width=5,height=8)
par(mfrow=c(3,1),mar=c(4,3,2,1))
for (j in 1:length(pops))
{
  tmp <- all[all$pop==pops[j],]
  cols.to.use <- chrs.cols$cols[match(tmp$chr,chrs.cols$chrs)]
  plot(tmp$fst,pch=20,cex=.75,col=as.character(cols.to.use),
       ylab="Fst",xlab="",xaxt="n",main=pops[j],ylim=c(0,1))
  x=c(tapply(1:nrow(tmp),tmp$chr,mean))
  minx=c(tapply(1:nrow(tmp),tmp$chr,min))-1
  mtext(levels(tmp$chr),at=x,las=3,cex=.8,side=1,line=.5,col=as.character(chrs.cols$cols))
  segments(x0 = minx,y0 = -10,x1 = minx,y1 = 10^3,col="grey")
  # quantiles
  segments(x0=-10,y0=quantile(tmp$fst,c(0.01,0.99)),x1=10^6,y1=quantile(tmp$fst,c(0.01,0.99)),col="black")
  points(tmp$fst,pch=20,cex=.75,col=as.character(cols.to.use))
  out.tmp <- data.frame(pop=pops[j],chr_pos=tmp$chr_pos[tmp$fst>quantile(tmp$fst,c(0.01,0.99))[2]])
  out <- rbind(out,out.tmp)
  stats <- rbind(stats,data.frame(pop=pops[j],mean=mean(tmp$fst),median=median(tmp$fst)))
}
dev.off()
print(stats)
