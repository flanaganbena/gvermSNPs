# shift in likelihood across 10 runs
rm(list=ls())
pdf('output/ngsAdmix-lnl.pdf',width=8,height=4)
ks <- c("k02","k03","k04","k05","k06","k07","k08","k09","k10","k11","k12",
        "k13","k14","k15","k16","k17","k18","k19","k20")
filename <- list.files(path="data/NGSadmix.runs/",pattern=".log")

lnl <- c()
for (k in 1:length(ks))
{
  filename.k <- filename[grep(pattern=ks[k],filename)]
  for (i in 1:length(filename.k))
  {
    liketmp <- readLines(paste("data/NGSadmix.runs/",filename.k[i],sep=""))
    liketmp <- liketmp[grep(pattern = "best like=",liketmp)]
    #print(filename.k[i])#;print(liketmp)
    lnl <- rbind(lnl,
                 data.frame(k=ks[k],run=i,lnl=substr(liketmp,regexpr(pattern="=",liketmp)[1]+1,regexpr(pattern="=",liketmp)[1]+14)))
    
  }
}
lnl$lnl <- as.numeric(as.character(lnl$lnl))
plot(lnl~factor(k),data=lnl)
xbar <- tapply(lnl$lnl,lnl$k,mean)
std <- tapply(lnl$lnl,lnl$k,sd)
out <- data.frame(xbar,std)
# Evanno et al. 2005 ∆K = m(|L(K + 1) − 2 L(K ) + L(K − 1)|)/s[L(K )]
out$L.prime.k <- c(NA,xbar[-1]-xbar[-length(xbar)])
out$L.dblprime.k <- c(NA,out$L.prime.k[-c(1,length(xbar))]-(out$L.prime.k)[-(1:2)],NA)
out$delta <- out$L.dblprime.k/out$std
print(out)
plot(out$delta,xaxt="n",xlab="k",type="b",ylab="delta lnl");mtext(at=1:length(xbar),rownames(out),side = 1)
plot(out$L.dblprime.k,xaxt="n",xlab="k",type="b", ylab="k''");mtext(at=1:length(xbar),rownames(out),side = 1)

dev.off()