### halfdecay ~ chromosome length X region

rm(list=ls())
library(scales) # transparency
library(ggplot2)
library(lmerTest)
chrs <- read.delim('data/Table-chr.key.txt')
meta <- read.csv('data/meta-SNP.pops_edited.V2.csv')
stats <- read.csv('output/ngsLD.stats.500bpBINS.csv')
colnames(stats) <- c("pop","chr","halfdecay","halfdecaydist","dist.LD.is.10perc")
#stats$chr <- paste("sca",substr(stats$chr,4,5),sep="")
stats$natnon <- meta$Country[match(stats$pop,meta$Site.Abb)]
stats$natnon <- factor(stats$natnon)

stats$size.kbp <- chrs$size.bp[match(stats$chr,chrs$name1)]/1000

### Pull 3 outlier points (RIS)
#stats <- stats[stats$halfdecaydist<1000,]

m <- lmer(log(halfdecaydist)~size.kbp*natnon+(1|pop),data=stats)
print(anova(m))

m2 <- lmer(log(halfdecaydist)~size.kbp+(1|pop),data=stats[stats$natnon=="Japan",])
print(anova(m2))
m3 <- lmer(log(halfdecaydist)~size.kbp+(1|pop),data=stats[stats$natnon=="Non-native",])
print(anova(m3))

#pop.unique <- unique(stats$pop)
#for(i in 1:22)
#{
#  tmp <- stats[stats$pop==pop.unique[i],]
#  abline(lm(halfdecaydist~size.kbp,data=tmp),
#         col=c("black","red")[tmp$natnon])
#}

gg <- ggplot(stats, aes(x=size.kbp, y = halfdecaydist,color=natnon)) +
  geom_point(colour=alpha(c("black","red")[stats$natnon],.5)) +
  ylim(0,200) +
  geom_smooth(method="lm") + #,color="grey",fill="red")+
  scale_colour_manual(values=alpha(c("black","red"),.5)) +
  xlab("Scaffold length (kbp)") +
  ylab("Half-decay distance (kbp)") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(.~natnon)
pdf('output/halfdecay.by.chrLength+natnon.500bpBINS.pdf',width=6,height=4)
print(gg)
dev.off()
