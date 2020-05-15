# combine simulation and empirical data

sink('output/ANOVA.txt')

############################
#### half-decay distance ####
############################
library(ggplot2)
library(gridExtra)
library(dplyr)
rm(list=ls())
meta <- read.csv('data/meta-SNP.pops_edited.V2.csv')

# simulation #
# 99% bottleneck; no clonality; 0 --> 200 generations #
bn_99 <- read.table("data/simulations/bottleneck_99.half.decay", header = T)
bn_99 <- bn_99[grepl("rep", bn_99$outnam),]
bn_99$outnam <- gsub("pre", "pre_0", as.character(bn_99$outnam))
tmp <- strsplit(bn_99$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",4)) ### check this is the right column with (table())
bn_99$gen <- tmp2
bn_99$gen <- sub("\\<01\\>", "1", bn_99$gen) ### Change values for plotting.
bn_99$clonality <- "99% Bottleneck + No clonality"
bn_99$halfdecaydist.kb <- bn_99$halfdecaydist/1000
bn_99.subset <- bn_99[bn_99$gen%in%c("0","1","10","100","200"),]

# 99% bottleneck; 100% clonality
bn_99_cl_100 <- read.table("data/simulations/bottleneck_99_clonality_100.half.decay", header = T)
bn_99_cl_100 <- bn_99_cl_100[grepl("rep", bn_99_cl_100$outnam),]
bn_99_cl_100$outnam <- gsub("pre", "pre_0", as.character(bn_99_cl_100$outnam))
tmp <- strsplit(bn_99_cl_100$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",6)) ### check this is the right column with (table())
bn_99_cl_100$gen <- tmp2
bn_99_cl_100$gen <- sub("\\<01\\>", "1", bn_99_cl_100$gen) ### Change values for plotting.
bn_99_cl_100$clonality <- "99% Bottleneck + 100% clonality"
bn_99_cl_100$halfdecaydist.kb <- bn_99_cl_100$halfdecaydist/1000
bn_99_cl_100.subset <- bn_99_cl_100[bn_99_cl_100$gen%in%c("0","1","10","100","200"),]

# no bottleneck; 100% clonality
bn_no_cl_100 <- read.table('data/simulations/clonality_100.half.decay',header=T)
bn_no_cl_100 <- bn_no_cl_100[grepl("rep", bn_no_cl_100$outnam),]
bn_no_cl_100$outnam <- gsub("pre", "pre_0", as.character(bn_no_cl_100$outnam))
tmp <- strsplit(bn_no_cl_100$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",4)) ### check this is the right column with (table())
bn_no_cl_100$gen <- tmp2
bn_no_cl_100$gen <- sub("\\<01\\>", "1", bn_no_cl_100$gen) ### Change values for plotting
bn_no_cl_100$clonality <- "No bottleneck + 100% clonality"
bn_no_cl_100$halfdecaydist.kb <- bn_no_cl_100$halfdecaydist/1000
bn_no_cl_100.subset <- bn_no_cl_100[bn_no_cl_100$gen%in%c("0","1","10","100","200"),]


# combine datasets #
sim1 <- rbind(bn_99.subset[,c("gen","clonality","halfdecaydist.kb")],
             bn_99_cl_100.subset[,c("gen","clonality","halfdecaydist.kb")],
             bn_no_cl_100.subset[,c("gen","clonality","halfdecaydist.kb")])

## relabel "generations"
sim1$gen[sim1$gen=="0"] <- "Before"
sim1$gen[sim1$gen=="1"] <- "0"
sim1$gen <- factor(sim1$gen,levels=c("Before","0","1","10","100","200"))

gg1a <- ggplot(sim1, aes(x=gen, y = halfdecaydist.kb, fill = clonality)) +
  geom_boxplot(outlier.shape=20,outlier.size = .8)+
  xlab("Generation") +
  ylab("") +
  ylim(0, 140) +
  scale_fill_grey(start=.5,end=1) +
  theme_bw() + 
  theme(legend.position = c(.65, .85),legend.title = element_blank(), legend.background = element_blank())

# empirical data #
emp <- read.csv('output/ngsLD.stats.500bpBINS.csv')
emp1 <- aggregate(emp$halfdecaydist,by=list(emp$pop),mean,na.rm=T)
colnames(emp1) <- c("pop","halfdecay")
emp1$natnon <- meta$Country[match(emp1$pop,meta$Site.Abb)]
gg1b <- ggplot(emp1,aes(x=natnon,y=halfdecay,fill=natnon)) +
  geom_boxplot(outlier.shape=20,outlier.size = .8,color=alpha(c("black","red"),.75))+
  geom_jitter(position=position_jitter(width=0.1),color=c("black","red")[emp1$natnon]) +
  xlab("") +
  ylab("Half-decay distance (kbp)") +
  ylim(0, 140) +
  scale_fill_manual(values=alpha(c("black","red"),.1)) +
  theme_bw() + theme(legend.position = "none")

m <- lm(log(emp1$halfdecay)~emp1$natnon)
print(anova(m))
#################################
#### Observed heterozygosity ####
################################

# 99% bottleneck; no clonality; 0 --> 1000 generations #
bn_99 <- read.table("data/simulations/bottleneck_99.het.gz", header = T)
bn_99$hobs <- (bn_99$N_SITES-bn_99$O.HOM.)/bn_99$N_SITES
bn_99 <- bn_99[grepl("rep", bn_99$outnam),]
bn_99$outnam <- gsub("pre", "pre_0", as.character(bn_99$outnam))
tmp <- strsplit(bn_99$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",4)) ### check this is the right column with (table())
bn_99$gen <- tmp2
bn_99$gen <- sub("\\<01\\>", "1", bn_99$gen) ### Change values for plotting.
bn_99$rep <- unlist(lapply(tmp,"[[",5)) # there were ~7000 sliding window 
bn_99.subset <- bn_99[bn_99$gen%in%c("0","1","10","100","200"),]
bn_99.xbar <- aggregate(bn_99.subset$hobs,by=list(bn_99.subset$rep,bn_99.subset$gen),mean,na.rm=T)
bn_99.xbar <- data.frame(bn_99.xbar,clonality="99% Bottleneck + No clonality")

# 99% bottleneck; 100% clonality
bn_99_cl_100 <- read.table("data/simulations/bottleneck_99_clonality_100.het.gz", header = T)
bn_99_cl_100$hobs <- (bn_99_cl_100$N_SITES-bn_99_cl_100$O.HOM.)/bn_99_cl_100$N_SITES
bn_99_cl_100 <- bn_99_cl_100[grepl("rep", bn_99_cl_100$outnam),]
bn_99_cl_100$outnam <- gsub("pre", "pre_0", as.character(bn_99_cl_100$outnam))
tmp <- strsplit(bn_99_cl_100$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",6)) ### check this is the right column with (table())
bn_99_cl_100$gen <- tmp2
bn_99_cl_100$gen <- sub("\\<01\\>", "1", bn_99_cl_100$gen) ### Change values for plotting.
bn_99_cl_100$rep <- unlist(lapply(tmp,"[[",7)) # there were ~7000 sliding window 
bn_99_cl_100.subset <- bn_99_cl_100[bn_99_cl_100$gen%in%c("0","1","10","100","200"),]
bn_99_cl_100.xbar <- aggregate(bn_99_cl_100.subset$hobs,by=list(bn_99_cl_100.subset$rep,bn_99_cl_100.subset$gen),mean,na.rm=T)
bn_99_cl_100.xbar <- data.frame(bn_99_cl_100.xbar,clonality="99% Bottleneck + 100% clonality")

# no bottleneck; 100% clonality
bn_no_cl_100 <- read.table('data/simulations/clonality_100.het.gz',header=T)
bn_no_cl_100$hobs <- (bn_no_cl_100$N_SITES-bn_no_cl_100$O.HOM.)/bn_no_cl_100$N_SITES
bn_no_cl_100 <- bn_no_cl_100[grepl("rep", bn_no_cl_100$outnam),]
bn_no_cl_100$outnam <- gsub("pre", "pre_0", as.character(bn_no_cl_100$outnam))
tmp <- strsplit(bn_no_cl_100$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",4)) ### check this is the right column with (table())
bn_no_cl_100$gen <- tmp2
bn_no_cl_100$gen <- sub("\\<01\\>", "1", bn_no_cl_100$gen) ### Change values for plotting.
bn_no_cl_100$rep <- unlist(lapply(tmp,"[[",5)) # there were ~7000 sliding window 
bn_no_cl_100.subset <- bn_no_cl_100[bn_no_cl_100$gen%in%c("0","1","10","100","200"),]
bn_no_cl_100.xbar <- aggregate(bn_no_cl_100.subset$hobs,by=list(bn_no_cl_100.subset$rep,bn_no_cl_100.subset$gen),mean,na.rm=T)
bn_no_cl_100.xbar <- data.frame(bn_no_cl_100.xbar,clonality="No Bottleneck + 100% clonality")


# combine datasets #
sim2 <- rbind(bn_99.xbar,bn_99_cl_100.xbar,bn_no_cl_100.xbar)
colnames(sim2) <- c("rep","gen","hobs","clonality")

## relabel "generations"
sim2$gen[sim2$gen=="0"] <- "Before"
sim2$gen[sim2$gen=="1"] <- "0"
sim2$gen <- factor(sim2$gen,levels=c("Before","0","1","10","100","200"))

gg2a <- ggplot(sim2, aes(x=gen, y = hobs, fill = clonality)) +
  geom_boxplot(outlier.shape=20,outlier.size = .8)+
  xlab("Generation") +
  ylab("") +
  #ylim(0, 20) +
  scale_fill_grey(start=.5,end=1) +
  theme_bw() +
  theme(legend.position = "none")

# empirical #
emp2 <- read.csv('data/Hobs.means.csv')
#poptouse <- c("aka","ave","cfj","cnt","dhd","dns","fdm","fut","hay","hik","lhp","mng","nyc","osh","ris","shk","sou","tmb","uho","usu","waj","wwf")
#emp2 <- emp2[emp2$pop%in%poptouse,]
gg2b <- ggplot(emp2,aes(x=natnon,y=Hobs,fill=natnon)) +
  geom_boxplot(outlier.shape=20,outlier.size = .8,color=alpha(c("black","red"),.75))+
  geom_jitter(position=position_jitter(width=0.1),color=c("black","red")[emp2$natnon]) +
  xlab("") +
  ylab("Ho") +
  #ylim(0, 0.25) +
  scale_fill_manual(values=alpha(c("black","red"),.1)) +
  theme_bw() + theme(legend.position = "none")

m <- lm(emp2$Hobs~emp2$natnon)
print(anova(m))
#################################
#### Tajima D ####
################################

# 99% bottleneck; no clonality; 0 --> 1000 generations #
bn_99 <- read.table("data/simulations/bottleneck_99.tajimasD.gz", header = T)
bn_99 <- bn_99[grepl("rep", bn_99$outnam),]
bn_99$outnam <- gsub("pre", "pre_0", as.character(bn_99$outnam))
tmp <- strsplit(bn_99$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",4)) ### check this is the right column with (table())
bn_99$gen <- tmp2
bn_99$gen <- sub("\\<01\\>", "1", bn_99$gen) ### Change values for plotting.
bn_99$rep <- unlist(lapply(tmp,"[[",5)) # there were ~7000 sliding window 
bn_99.subset <- bn_99[bn_99$gen%in%c("0","1","10","100","200"),]
bn_99.xbar <- aggregate(bn_99.subset$TajimaD,by=list(bn_99.subset$rep,bn_99.subset$gen),mean,na.rm=T)
bn_99.xbar <- data.frame(bn_99.xbar,clonality="99% Bottleneck + No clonality")

# 99% bottleneck; 100% clonality
bn_99_cl_100 <- read.table("data/simulations/bottleneck_99_clonality_100.tajimasD.gz", header = T)
bn_99_cl_100 <- bn_99_cl_100[grepl("rep", bn_99_cl_100$outnam),]
bn_99_cl_100$outnam <- gsub("pre", "pre_0", as.character(bn_99_cl_100$outnam))
tmp <- strsplit(bn_99_cl_100$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",6)) ### check this is the right column with (table())
bn_99_cl_100$gen <- tmp2
bn_99_cl_100$gen <- sub("\\<01\\>", "1", bn_99_cl_100$gen) ### Change values for plotting.
bn_99_cl_100$rep <- unlist(lapply(tmp,"[[",7)) # there were ~7000 sliding window 
bn_99_cl_100.subset <- bn_99_cl_100[bn_99_cl_100$gen%in%c("0","1","10","100","200"),]
bn_99_cl_100.xbar <- aggregate(bn_99_cl_100.subset$TajimaD,by=list(bn_99_cl_100.subset$rep,bn_99_cl_100.subset$gen),mean,na.rm=T)
bn_99_cl_100.xbar <- data.frame(bn_99_cl_100.xbar,clonality="99% Bottleneck + 100% clonality")

# no bottleneck; 100% clonality
bn_no_cl_100 <- read.table("data/simulations/clonality_100.tajimasD.gz", header = T)
bn_no_cl_100 <- bn_no_cl_100[grepl("rep", bn_no_cl_100$outnam),]
bn_no_cl_100$outnam <- gsub("pre", "pre_0", as.character(bn_no_cl_100$outnam))
tmp <- strsplit(bn_no_cl_100$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",4)) ### check this is the right column with (table())
bn_no_cl_100$gen <- tmp2
bn_no_cl_100$gen <- sub("\\<01\\>", "1", bn_no_cl_100$gen) ### Change values for plotting.
bn_no_cl_100$rep <- unlist(lapply(tmp,"[[",5)) # there were ~7000 sliding window 
bn_no_cl_100.subset <- bn_no_cl_100[bn_no_cl_100$gen%in%c("0","1","10","100","200"),]
bn_no_cl_100.xbar <- aggregate(bn_no_cl_100.subset$TajimaD,by=list(bn_no_cl_100.subset$rep,bn_no_cl_100.subset$gen),mean,na.rm=T)
bn_no_cl_100.xbar <- data.frame(bn_no_cl_100.xbar,clonality="No Bottleneck + 100% clonality")

# combine datasets #
sim3 <- rbind(bn_99.xbar,bn_99_cl_100.xbar,bn_no_cl_100.xbar)
colnames(sim3) <- c("rep","gen","TajimaD","clonality")

## relabel "generations"
sim3$gen[sim3$gen=="0"] <- "Before"
sim3$gen[sim3$gen=="1"] <- "0"
sim3$gen <- factor(sim3$gen,levels=c("Before","0","1","10","100","200"))

gg3a <- ggplot(sim3, aes(x=gen, y = TajimaD, fill = clonality)) +
  geom_boxplot(outlier.shape=20,outlier.size = .8)+
  xlab("Generation") +
  ylab("") +
  ylim(-1, 3) +
  scale_fill_grey(start=.5,end=1) +
  theme_bw() + 
  theme(legend.position = "none")

# empirical #
emp <- read.csv('data/TajimaD.means.csv')
emp3 <- aggregate(emp$Tajima,by=list(emp$pop),mean,na.rm=T)
colnames(emp3) <- c("pop","TajimaD")
emp3$natnon <- meta$Country[match(emp3$pop,meta$Site.Abb)]

gg3b <- ggplot(emp3,aes(x=natnon,y=TajimaD,fill=natnon)) +
  geom_boxplot(outlier.shape=20,outlier.size = .8,color=alpha(c("black","red"),.75))+
  geom_jitter(position=position_jitter(width=0.1),color=c("black","red")[emp3$natnon]) +
  xlab("") +
  ylab("Tajima's D") +
  ylim(-1, 3) +
  scale_fill_manual(values=alpha(c("black","red"),.1)) +
  theme_bw() + theme(legend.position = "none")

m <- lm(emp3$TajimaD~emp3$natnon)
print(anova(m))

#################################
#### Sites Pi ####
################################

# 99% bottleneck; no clonality; 0 --> 1000 generations #
bn_99 <- read.table("data/simulations/bottleneck_99.sites.pi.gz", header = T)
bn_99 <- bn_99[grepl("rep", bn_99$outnam),]
bn_99$outnam <- gsub("pre", "pre_0", as.character(bn_99$outnam))
tmp <- strsplit(bn_99$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",4)) ### check this is the right column with (table())
bn_99$gen <- tmp2
bn_99$gen <- sub("\\<01\\>", "1", bn_99$gen) ### Change values for plotting.
bn_99$rep <- unlist(lapply(tmp,"[[",5)) # there were ~7000 sliding window 
bn_99.subset <- bn_99[bn_99$gen%in%c("0","1","10","100","200"),]
bn_99.xbar <- aggregate(bn_99.subset$PI,by=list(bn_99.subset$rep,bn_99.subset$gen),mean,na.rm=T)
bn_99.xbar <- data.frame(bn_99.xbar,clonality="99% Bottleneck + no clonality")

# 99% bottleneck; 100% clonality
bn_99_cl_100 <- read.table("data/simulations/bottleneck_99_clonality_100.sites.pi.gz", header = T)
bn_99_cl_100 <- bn_99_cl_100[grepl("rep", bn_99_cl_100$outnam),]
bn_99_cl_100$outnam <- gsub("pre", "pre_0", as.character(bn_99_cl_100$outnam))
tmp <- strsplit(bn_99_cl_100$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",6)) ### check this is the right column with (table())
bn_99_cl_100$gen <- tmp2
bn_99_cl_100$gen <- sub("\\<01\\>", "1", bn_99_cl_100$gen) ### Change values for plotting.
bn_99_cl_100$rep <- unlist(lapply(tmp,"[[",7)) # there were ~7000 sliding window 
bn_99_cl_100.subset <- bn_99_cl_100[bn_99_cl_100$gen%in%c("0","1","10","100","200"),]
bn_99_cl_100.xbar <- aggregate(bn_99_cl_100.subset$PI,by=list(bn_99_cl_100.subset$rep,bn_99_cl_100.subset$gen),mean,na.rm=T)
bn_99_cl_100.xbar <- data.frame(bn_99_cl_100.xbar,clonality="99% Bottleneck + 100% clonality")

# no bottleneck; 100% clonality
bn_no_cl_100 <- read.table("data/simulations/clonality_100.sites.pi.gz", header = T)
bn_no_cl_100 <- bn_no_cl_100[grepl("rep", bn_no_cl_100$outnam),]
bn_no_cl_100$outnam <- gsub("pre", "pre_0", as.character(bn_no_cl_100$outnam))
tmp <- strsplit(bn_no_cl_100$outnam,"_")
tmp2 <- unlist(lapply(tmp,"[[",4)) ### check this is the right column with (table())
bn_no_cl_100$gen <- tmp2
bn_no_cl_100$gen <- sub("\\<01\\>", "1", bn_no_cl_100$gen) ### Change values for plotting.
bn_no_cl_100$rep <- unlist(lapply(tmp,"[[",5)) # there were ~7000 sliding window 
bn_no_cl_100.subset <- bn_no_cl_100[bn_no_cl_100$gen%in%c("0","1","10","100","200"),]
bn_no_cl_100.xbar <- aggregate(bn_no_cl_100.subset$PI,by=list(bn_no_cl_100.subset$rep,bn_no_cl_100.subset$gen),mean,na.rm=T)
bn_no_cl_100.xbar <- data.frame(bn_no_cl_100.xbar,clonality="No Bottleneck + 100% clonality")

# combine datasets #
sim4 <- rbind(bn_99.xbar,bn_99_cl_100.xbar,bn_no_cl_100.xbar)
colnames(sim4) <- c("rep","gen","PI","clonality")

## relabel "generations"
sim4$gen[sim4$gen=="0"] <- "Before"
sim4$gen[sim4$gen=="1"] <- "0"
sim4$gen <- factor(sim4$gen,levels=c("Before","0","1","10","100","200"))

gg4a <- ggplot(sim4, aes(x=gen, y = PI, fill = clonality)) +
  geom_boxplot(outlier.shape=20,outlier.size = .8)+
  xlab("Generation") +
  ylab("") +
  ylim(0, 0.30) +
  scale_fill_grey(start=.5,end=1) +
  theme_bw() + 
  theme(legend.position = "none")

# empirical #
emp4 <- read.csv('data/sites.pi.means.csv')
#poptouse <- c("aka","ave","cfj","cnt","dhd","dns","fdm","fut","hay","hik","lhp","mng","nyc","osh","ris","shk","sou","tmb","uho","usu","waj","wwf")
#emp4 <- emp4[emp4$pop%in%poptouse,]
gg4b <- ggplot(emp4,aes(x=natnon,y=sites.pi,fill=natnon)) +
  geom_boxplot(outlier.shape=20,outlier.size = .8,color=alpha(c("black","red"),.75))+
  geom_jitter(position=position_jitter(width=0.1),color=c("black","red")[emp4$natnon]) +
  xlab("") +
  ylab("PI") +
  ylim(0, 0.30) +
  scale_fill_manual(values=alpha(c("black","red"),.1)) +
  theme_bw() + theme(legend.position = "none")

m <- lm(emp4$sites.pi~emp4$natnon)
print(anova(m))


pdf('output/simulation+empirical.pdf',height=11,width=6)
grid.arrange(gg1b,gg1a,
             gg2b,gg2a,
             gg3b,gg3a,
             gg4b,gg4a,
             widths=c(1,2),nrow=4,ncol=2)
dev.off()
sink()