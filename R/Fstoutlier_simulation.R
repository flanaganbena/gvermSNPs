### pairwise Fst between mng/sou and three non-native populations

## sliding window ngsFST
### mngsou vs non-native
### 50000 bp window; 10000 steps

library(reshape2); library(dplyr)
rm(list=ls())

dat <- read.table("data/mast_fst_out.txt.gz", header = T)
dat <- dat %>% filter(rep != "rep")
dat$rep2 <- substr(dat$rep, 4, 10)

goi <- c(1, 25, 100)

out <- c()
stats <- c()
pdf('output/Fstoutlier_simulation.pdf',width=5,height=8, compress = T)
par(mfrow=c(3,1),mar=c(4,3,2,1))
for (j in unique(goi)){
  set.seed(1900)
  tmp <-  dat %>% filter(gen == j) %>% filter(rep2 %in% sample(unique(dat$rep2),26)) %>% arrange(rep2)
  tmp$fst <- as.numeric(as.character(tmp$fst))
  rep2.cols <- data.frame(rep2=unique(tmp$rep2),cols=c(rep(c("red","blue"),13)))
  cols.to.use <- rep2.cols$cols[match(tmp$rep2,rep2.cols$rep2)]
  plot(tmp$fst, pch=20, cex=.75, col=as.character(cols.to.use),
       ylab="Fst",xlab="",xaxt="n",main=j,ylim=c(0,1))
  x=na.omit(c(tapply(1:nrow(tmp),tmp$rep2, mean)))
  minx=c(tapply(1:nrow(tmp),tmp$rep2,min))-1
  mtext(rep2.cols$rep2, at=x, las=3, cex=.8, side=1, line=0.5, col= as.character(rep2.cols$cols))
  segments(x0 = minx,y0 = -10,x1 = minx,y1 = 10^3,col="grey")
  # quantiles
  segments(x0=-10,y0=quantile(tmp$fst,c(0.01,0.99)),x1=10^6,y1=quantile(tmp$fst,c(0.01,0.99)),col="black")
  points(tmp$fst,pch=20,cex=.75,col=as.character(cols.to.use))
  out.tmp <- data.frame(gen=j,chr_pos=tmp$pos[tmp$fst>quantile(tmp$fst,c(0.01,0.99))[2]])
  out <- rbind(out,out.tmp)
  stats <- rbind(stats,data.frame(gen=j,mean=mean(tmp$fst),median=median(tmp$fst)))
}
dev.off()
print(stats)
