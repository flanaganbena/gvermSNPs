### get mean Fst by chr and pairwise populations
library(reshape)
library(RColorBrewer)
library(lattice)
rm(list=ls())
meta <- read.csv('data/meta-SNP.pops_edited.V2.csv')
fst <- c()
a <- list.files('data/fst/',pattern=".txt")
a.unique <- unique(a) 
for(i in 1:length(a.unique))
{
  tmp <- c()
  tmp <- read.delim(paste("data/fst/",a[i],sep=""),header=F,skip=1)
  colnames(tmp) <- c("region","chr","midPos","Nsites","fst")
  sites <- strsplit(a[i],".pbs.fst.idx.txt")[[1]][1]
  fst <- rbind(fst,data.frame(pops=sites,meanfst = mean(tmp$fst,na.rm=T)))
}
write.csv(fst,'output/meanfst.csv')

fst <- fst[-grep(pattern="mngsou",x=as.character(fst$pops)),]

#make heatmap

fst$pop1 <- substr(fst$pops,1,3)
fst$pop2 <- substr(fst$pops,5,7)
fst$meta1 <- meta$Country[match(fst$pop1,meta$Site.Abb)]
fst$meta2 <- meta$Country[match(fst$pop2,meta$Site.Abb)]
fst$meta12 <- paste(fst$meta1,fst$meta2)
fst$meta12[fst$meta12=="Japan Non-native"] <- "Non-native Japan"

md <- fst[,c("pop1","pop2","meanfst")]
md$pop1 <- factor(md$pop1); md$pop2 <- factor(md$pop2)
md <- melt(fst,id=c("pop1","pop2"),measure.vars = 'meanfst')

siteorder <- data.frame(pop = unique(c(as.character(md$pop1),as.character(md$pop2))))
siteorder$reg = meta$Region[match(siteorder$pop,meta$Site.Abb)]
siteorder$reg <- factor(siteorder$reg)
siteorder$reg <- factor(siteorder$reg,levels=levels(siteorder$reg)[c(3,4,1,2)])
siteorder$lat <- meta$Latitude[match(siteorder$pop,meta$Site.Abb)]
siteorder <- siteorder[order(siteorder$lat),]
siteorder <- siteorder[order(siteorder$reg),]
siteorder$pop2 <- c(paste("0",1:9,siteorder$pop[1:9],sep=""),paste(10:22,siteorder$pop[10:22],sep=""))

full <- c()
for (i in 1:length(siteorder$pop))
{
  tmp <- md[md$pop1==siteorder$pop[i] | md$pop2==siteorder$pop[i],]
  tmp <- rbind(tmp, data.frame(pop1=siteorder$pop[i],pop2=siteorder$pop[i],variable="meanfst",value=0))
  tmp$othersite <- ifelse(tmp$pop1==siteorder$pop[i],tmp$pop2,tmp$pop1)
  tmp <- tmp[order(match(tmp$othersite,siteorder$pop)),]
  full <- cbind(full,tmp$value)
}
### get corrected names
pop.correctednames <- meta$Site.Abb.2[match(siteorder$pop,meta$Site.Abb)]
rownames(full) <- pop.correctednames
colnames(full) <- pop.correctednames
col.l <- colorRampPalette(brewer.pal(11, "RdBu")[11:1])
f <- levelplot(full,cex=.3,col.regions=col.l)
pdf('output/heatmap.pdf',width=10,height=8)
print(f)
dev.off()

