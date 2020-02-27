### do assignment test based on linear discriminant analysis
### analogous to DAPC but for the output from PCAngsd
library(reshape)
library("circlize");library(RColorBrewer)
library(MASS) ### lda()
######## SNPs ###########

### 1- input the entire PCA
rm(list=ls())

meta <- read.csv('data/meta-SNP.pops_edited.V2.csv')
ind <- read.delim('data/dipsmill.351inds.txt',header=F)
pop <- substr(ind[,1],1,3)
pop <- meta$Site.Abb.2[match(pop,meta$Site.Abb)]
reg <- meta$Region[match(pop,meta$Site.Abb.2)]; reg <- factor(reg)
nat.non <- meta$Country[match(pop,meta$Site.Abb.2)]; nat.non <- factor(nat.non)
dat <- read.delim('data/dipsmll.351ind.cov',header=F)
e <- eigen(dat)
pca1 <- e$vectors[,1]; pca2 <- e$vectors[,2]; pca3 <- e$vectors[,3]; pca4 <- e$vectors[,4]

### 2- use the Native range to train the LDA

all <- data.frame(pop,pca1,pca2,pca3,pca4)
#pop.nat <- pop[nat.non=="Non-native"]
z <- lda(pop~.,all,subset=(1:351)[nat.non=="Japan"])
out <- predict(z,all[(1:351)[nat.non=="Non-native"],])
out2 <- out$posterior
out3 <- data.frame(out2)
out3$reg <- substr(reg[nat.non=="Non-native"],1,3)
md <- melt(out3)
md2 <- cast(reg~variable,data=md,sum)
md2 <- data.frame(md2)
rownames(md2) <- md2[,1]
md2 <- md2[,-1]
md3 <- data.frame(t(md2))
md3 <- cbind(md3,data.frame(matrix(rep(0,11*11),ncol=11)))
tmp <- data.frame(matrix(rep(0,14*3),nrow=3));rownames(tmp) <- c("ena","eur","wna");colnames(tmp) <- colnames(md3)
md3 <- rbind(md3,tmp)
md3 <- as.matrix(md3)

### 3-make plot
png('output/assignment-snp.png',res = 700,units = "in",height=10,width=10)
circos.clear()
par(mar = rep(0, 4), cex=0.9)
circos.par(start.degree = 90, gap.degree = 4)
cols.to.use <- c(rep("black",3),rep("red",3))
#cols.to.use[rownames(out3std)%in%c("akk","mng","sou","fut")] <- "black"
fig <- chordDiagram(x = md3, directional = 1, #order = order(toSortCols), 
             grid.col=cols.to.use,
             annotationTrack = "grid", 
             transparency = 0.25,  annotationTrackHeight = c(0.1, 0.1),
             #diffHeight  = -0.04,
             direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",
             link.arr.length =  0.15)

 dev.off()
 pdf('output/assignment-snps-details.pdf')
 fig <- chordDiagram(x = md3, directional = 1)
 dev.off()
 
