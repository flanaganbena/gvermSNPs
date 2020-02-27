### mean half decay dist
rm(list=ls())
library(scales) # transparency
library(ggplot2)
library(lmerTest)
meta <- read.csv('data/meta-SNP.pops_edited.V2.csv')
stats <- read.csv('output/ngsLD.stats.500bpBINS.csv')
colnames(stats) <- c("pop","chr","halfdecay","halfdecaydist","dist.LD.is.10perc")
stats$chr <- paste("sca",substr(stats$chr,4,5),sep="")
stats$natnon <- meta$Country[match(stats$pop,meta$Site.Abb)]
stats$natnon <- ifelse(stats$natnon=="Japan","nat","non")
stats$pop2 <- meta$Site.Abb.2[match(stats$pop,meta$Site.Abb)] ### corrected name

m <- lmer(log(halfdecaydist)~chr*natnon+(1|pop),data=stats)
print(anova(m))

### show populations on x-axis

lat <- meta$Latitude[match(levels(stats$pop2),meta$Site.Abb.2)]
reg <- meta$Region[match(levels(stats$pop2),meta$Site.Abb.2)]
reg <- factor(reg)
reg <- factor(reg,levels=levels(reg)[c(3,4,1,2)])
stats$pop2 <- factor(stats$pop2,levels=levels(stats$pop2)[order(reg)])
gg <- ggplot(stats, aes(x=pop2, y = halfdecaydist,fill=natnon)) +
  geom_boxplot() +#outlier.shape=NA)+
  ylim(0, 400) +
  xlab("") +
  ylab("Half-decay distance (kbp)") +
  scale_fill_manual(values=alpha(c("black","red"),.5)) +
  theme_bw() +
  theme(legend.position = "none")
pdf('output/halfdecay.by.natnon.500bpBINS-bypop.pdf',width=6,height=4)
  print(gg)
  dev.off()

