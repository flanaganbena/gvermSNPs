

### Plots genotypic richness after IDing repeated genotypes in simulation vcffiles ###

library(stringr)
dat <- read.csv("data/simulations/genotypic_richness/cl_100_dupe_genomes.csv")
dat$fp <- gsub("pre", "pre_0", as.character(dat$fp))
dat$fp <- gsub("_01", "_1", as.character(dat$fp))

dat$cl <- str_split(dat$fp, pattern = "_", simplify = T)[,2]
dat$gen <- str_split(dat$fp, pattern = "_", simplify = T)[,4]
dat$sim <- str_split(dat$fp, pattern = "_", simplify = T)[,5]
dat <- tibble(dat)
dat$dupe_count <- as.numeric(as.character(dat$dupe_count))
dat <- na.omit(dat)
sum <- dat %>% group_by(fp) %>% summarise(R = length(dupe_count), N = sum(dupe_count))
sum$R.2 <- sum$R / sum$N
sum$cl <- str_split(sum$fp, pattern = "_", simplify = T)[,2]
sum$gen <- str_split(sum$fp, pattern = "_", simplify = T)[,4]
sum$sim <- str_split(sum$fp, pattern = "_", simplify = T)[,5]
sum$gen <- as.numeric(as.character(sum$gen))
sum$bn <- 0


dat2 <- read.csv("data/simulations/genotypic_richness/bn_100_cl_100_dupe_genomes.csv")
dat2$fp <- gsub("pre", "pre_0", as.character(dat2$fp))
dat2$fp <- gsub("_01", "_1", as.character(dat2$fp))

dat2$bn <- str_split(dat2$fp, pattern = "_", simplify = T)[,2]
dat2$cl <- str_split(dat2$fp, pattern = "_", simplify = T)[,4]
dat2$gen <- str_split(dat2$fp, pattern = "_", simplify = T)[,6]
dat2$sim <- str_split(dat2$fp, pattern = "_", simplify = T)[,7]
dat2 <- tibble(dat2)
dat2$dupe_count <- as.numeric(as.character(dat2$dupe_count))
dat2 <- na.omit(dat2)
sum2 <- dat2 %>% group_by(fp) %>% summarise(R = length(dupe_count), N = sum(dupe_count))
sum2$R.2 <- sum2$R / sum2$N
sum2$bn <- str_split(sum2$fp, pattern = "_", simplify = T)[,2]
sum2$cl <- str_split(sum2$fp, pattern = "_", simplify = T)[,4]
sum2$gen <- str_split(sum2$fp, pattern = "_", simplify = T)[,6]
sum2$sim <- str_split(sum2$fp, pattern = "_", simplify = T)[,7]
sum2$gen <- as.numeric(as.character(sum2$gen))


dat3 <- read.csv("data/simulations/genotypic_richness/null_dupe_genomes.csv")
dat3$fp <- gsub("pre", "pre_0", as.character(dat3$fp))
dat3$fp <- gsub("_01", "_1", as.character(dat3$fp))


dat3 <- tibble(dat3)
dat3$dupe_count <- as.numeric(as.character(dat3$dupe_count))
dat3 <- na.omit(dat3)
sum3 <- dat3 %>% group_by(fp) %>% summarise(R = length(dupe_count), N = sum(dupe_count))
sum3$R.2 <- sum3$R / sum3$N
sum3$bn <- 0
sum3$cl <- 0
sum3$gen <- str_split(sum3$fp, pattern = "_", simplify = T)[,3]
sum3$sim <- str_split(sum3$fp, pattern = "_", simplify = T)[,4]
sum3$gen <- as.numeric(as.character(sum3$gen))
sum3$gen <- ifelse(sum3$gen == "NA", 200, sum3$gen)


dat4 <- read.csv("data/simulations/genotypic_richness/bn_99_dupe_genomes.csv")
dat4$fp <- gsub("pre", "pre_0", as.character(dat4$fp))
dat4$fp <- gsub("_01", "_1", as.character(dat4$fp))


dat4 <- tibble(dat4)
dat4$dupe_count <- as.numeric(as.character(dat4$dupe_count))
dat4 <- na.omit(dat4)
sum5 <- dat4 %>% group_by(fp) %>% summarise(R = length(dupe_count), N = sum(dupe_count))
sum5$R.2 <- sum5$R / sum5$N
sum5$bn <- str_split(sum5$fp, pattern = "_", simplify = T)[,2]
sum5$cl <- 0
sum5$gen <- str_split(sum5$fp, pattern = "_", simplify = T)[,4]
sum5$sim <- str_split(sum5$fp, pattern = "_", simplify = T)[,5]
sum5$gen <- as.numeric(as.character(sum5$gen))
sum5$gen <- ifelse(sum5$gen == "NA", 200, sum5$gen)






sum_join <- rbind(sum, sum2, sum3, sum5) %>% na.omit() %>% filter(gen != 200)


sum_join$bn <- ifelse(sum_join$bn == "0", "No Bottleneck", sum_join$bn)
sum_join$bn <- ifelse(sum_join$bn == "99", "0.99", sum_join$bn)
sum_join$bn <- factor(sum_join$bn, levels = c("No Bottleneck", "0.99"))


sum_join$gen <- as.numeric(as.character(sum_join$gen))
sum_join$gen_plot <- ifelse(sum_join$gen == 0, "Before", sum_join$gen)
sum_join$gen_plot <- ifelse(sum_join$gen == 1, "0", sum_join$gen_plot)
sum_join$gen_plot <- factor(sum_join$gen_plot, levels = c("Before","0" ,"10", "25","50", "100", "200"))

sum_join$cl_2 <- as.character(as.numeric(as.character(sum_join$cl)) / 100)


a <- ggplot(sum_join, aes(x = as.factor(gen_plot), y = R.2, fill = cl_2)) +
  geom_boxplot() +
  facet_wrap(~bn) +
  ylim(c(0,1)) +
  scale_fill_grey(start = 1, end = 0) +
  xlab("Generation") +
  ylab("Genotypic richness (R)") +
  labs(fill = "Clonality") +
  theme_linedraw()


pdf('output/simulation_genotypicrichness.pdf',height=4.5,width=6)
print(a)
dev.off()
