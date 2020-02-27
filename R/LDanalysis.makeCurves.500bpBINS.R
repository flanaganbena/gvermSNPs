### analyze output from ngsLD

#ngsLD outputs a TSV file with LD results for all pairs of sites for which LD was calculated, where the first two columns are positions of the SNPs, the third column is the distance (in bp) between the SNPs, and the following 4 columns are the various measures of LD calculated (r^2 from pearson correlation between expected genotypes, D from EM algorithm, D' from EM algorithm, and r^2 from EM algorithm). 

# We generated mean r2 within 500 bp bins of distances using r-code.
# 3 columns: kb.bins, r2, kb.midpt

rm(list=ls())
library(ggplot2) # cut_interval()
load("data/LDanalysis.500bpBINS.rdata")
pop <- substr(names(out),1,3)
meta <- read.csv('data/meta-SNP.pops_edited.V2.csv')
nat.non <- meta$Country[match(pop,meta$Site.Abb)]

pdf('output/LDanalysis.makeCurves.500bpBINS.pdf',width=8,height=5)
#quartz(width=8,height=5)
par(mfrow=c(1,2),mar=c(2,2,2,2))
stats <- c()
for (i in 1:2)
{
  plot(1,1,xlim=c(0,3600),ylim=c(0,0.6),xlab="distance (kbp)", ylab="LD (r^2)",type="n")
  tmp <- out[nat.non==levels(nat.non)[i]]
  
for(j in 1:length(tmp))
{
  xbar <- tmp[[j]]
  Cstart <- c(C=0.1)
  CDist <- function(n,C,distance)
{
  ((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance)))
}
n=15
modelC = try(nls(r2 ~ CDist(n,C,kb.midpt), data=xbar, start=Cstart, control=nls.control(maxiter=100)))#
                 #error = function(e) modelC=print("oops"))
if(length(modelC)>1)
{
  xbar$prd <- predict(modelC)
  halfdecay = (max(xbar$prd))*0.5
  halfdecaydist <- xbar$kb.midpt[which.min(abs(xbar$prd-halfdecay))]
  dist.LD.is.10perc <- xbar$kb.midpt[which.min(abs(xbar$prd-0.1))] 
  stats <- rbind(stats,data.frame(pop=substr(names(tmp)[j],1,3),
                                      chr=substr(names(tmp)[j],5,9),
                                      halfdecay,
                                      halfdecaydist,
                                      dist.LD.is.10perc))
  lines(xbar$kb.midpt, xbar$prd, col=alpha(c("black","red")[i],.5), lwd=1)
}
else {stats <- rbind(stats,data.frame(pop=substr(names(tmp)[j],1,3),
                                          chr=substr(names(tmp)[j],5,9),
                                          halfdecaydist=NA,
                                          halfdecaydist=NA,
                                          dist.LD.is.10perc=NA))}
}
}

write.table(stats,"output/ngsLD.stats.500bpBINS.csv",sep=",",row.names=F,quote = F)

dev.off()
