### PCA

## covariance matrix from PCAngsd
### made eigenvectors and plotted
### Updated by BAF 5 Dec 2020 in response to reviewr comments 
rm(list=ls())
library(seqinr); library(ggplot2); library(dplyr); library(ggrepel) # col2alpha
meta <- read.csv('data/meta-SNP.pops_edited.V2.csv')
ind <- read.delim('data/dipsmill.351inds.txt',header=F)
dat <- read.delim('data/dipsmll.351ind.cov',header=F)
e <- eigen(dat)

### generate individal PCA data frame
pca_dat <- data.frame(e$vectors[,1], e$vectors[,2], e$vectors[,3], e$vectors[,4])
colnames(pca_dat) <- c("PC1", "PC2", "PC3", "PC4")
pca_dat$pop <- substr(ind[,1],1,3)
pca_dat$pop2 <- meta$Site.Abb.2[match(pca_dat$pop,meta$Site.Abb)]
pca_dat$nat.non <- meta$Country[match(pca_dat$pop2,meta$Site.Abb.2)]

### generate pop averages data frame
xbar1 <- tapply(pca_dat$PC1, pca_dat$pop2, mean)
xbar2 <- tapply(pca_dat$PC2, pca_dat$pop2, mean)
nat.non2 <- as.factor(meta$Country[match(names(xbar1),meta$Site.Abb.2)])
xbar3 <- tapply(pca_dat$PC3, pca_dat$pop2, mean)
xbar4 <- tapply(pca_dat$PC4, pca_dat$pop2, mean)
pca_pop_means <- data.frame(xbar1, xbar2, xbar3, xbar4, nat.non2)
pca_pop_means$pops <- rownames(pca_pop_means)
rownames(pca_pop_means) <- NULL
pca_pop_means$to_include_p1 <- ifelse(pca_pop_means$nat.non2 == "Japan" | pca_pop_means$pops %in% c("moo", "eld", "ptw"), "y", "n")
pca_pop_means$to_include_p2 <- ifelse(pca_pop_means$nat.non2 == "Non-native" | pca_pop_means$pops %in% c("mng", "sou", "mou"), "y", "n")


p1 <- ggplot(pca_dat, aes(x = PC1, y = PC2, color = nat.non)) +
  geom_point(alpha = 0.08) +
  geom_text_repel(data = pca_pop_means %>% filter(to_include_p1 == "y"),
                  aes(x = xbar1, y = xbar2, label = pops, color = nat.non2),
                  force = 0.05,
                  box.padding = 0.2,
                  fontface = "bold",
                  arrow = arrow(length = unit(0.006, "npc"), type = "closed", ends = "first")) +
  scale_color_manual(values=c("black", "red")) +
  theme_minimal() +
  theme(panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none")

p2 <- ggplot(pca_dat, aes(x = PC3, y = PC4, color = nat.non)) +
  geom_point(alpha = 0.08) +
  geom_text_repel(data = pca_pop_means %>% filter(to_include_p2 == "y"),
                  aes(x = xbar3, y = xbar4, label = pops, color = nat.non2),
                  force = 0.1,
                  box.padding = 0.3,
                  fontface = "bold",
                  arrow = arrow(length = unit(0.006, "npc"), type = "closed", ends = "first")) +  scale_color_manual(values=c("black", "red")) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none")

pdf("output/pca.pdf", width = 8, heigh = 5.2)
par(mar = c(3, 3, 2, 2))
p1
p2
dev.off()
      