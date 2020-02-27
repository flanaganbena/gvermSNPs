### PCA

## covariance matrix from PCAngsd
### made eigenvectors and plotted
rm(list=ls())
library(seqinr) # col2alpha
meta <- read.csv('data/meta-SNP.pops_edited.V2.csv')
ind <- read.delim('data/dipsmill.351inds.txt',header=F)
pop <- substr(ind[,1],1,3)
pop2 <- meta$Site.Abb.2[match(pop,meta$Site.Abb)]
nat.non <- meta$Country[match(pop2,meta$Site.Abb.2)]; nat.non <- factor(nat.non)
dat <- read.delim('data/dipsmll.351ind.cov',header=F)
e <- eigen(dat)
pca1 <- e$vectors[,1]; pca2 <- e$vectors[,2]; pca3 <- e$vectors[,3]; pca4 <- e$vectors[,4]

### generate pop averages
xbar1 <- tapply(pca1,pop2,mean)
xbar2 <- tapply(pca2,pop2,mean)
nat.non2 <- meta$Country[match(names(xbar1),meta$Site.Abb.2)]; nat.non2 <- factor(nat.non2)
xbar3 <- tapply(pca3,pop2,mean)
xbar4 <- tapply(pca4,pop2,mean)
pdf("output/pca.pdf",width=6,height=4)
par(mar=c(3,3,2,2))
plot(pca2~pca1,pch=20,col=c(col2alpha("black",.2),col2alpha("red",.2))[nat.non])
text(x=xbar1,y=xbar2,names(xbar1),col=c("black","red")[nat.non2])
plot(pca4~pca3,pch=20,col=c(col2alpha("black",.2),col2alpha("red",.2))[nat.non])
text(x=xbar3,y=xbar4,names(xbar1),col=c("black","red")[nat.non2])
dev.off()
