#!/usr/bin/evn Rscript

# Runs vcftools on vcf files. --TajimasD --het --freq2 --site-pi --geno-r2. Calculates Hill and Weir 1988 half-decay distance.  
# Usage: Rscript runVCFtools.R <vcf file>
# Built using VCFtools (0.1.15)

args <- commandArgs(TRUE)


setwd(".")
PATH <- getwd()
lenarg <- length(args)


for (i in 1:lenarg) {
  outnam <- strsplit(args, "[.]")[[i]][1]
  if(grepl("pre", outnam)== TRUE) {
	system(paste("vcftools --vcf ", outnam, ".vcf --kept-sites", sep = ""))
	}
}

outtajD <- NULL
outhet <- NULL
outfreq2 <- NULL
outld <- NULL
outsitespi <- NULL
for (j in 1:lenarg) {
  outnam <- strsplit(args, "[.]")[[j]][1]
  system(paste("vcftools --vcf ", outnam, ".vcf --TajimaD 2000 --positions out.kept.sites --out ", outnam, sep = ""))
  system(paste("vcftools --vcf ", outnam, ".vcf --het --positions out.kept.sites --out ", outnam, sep = ""))
  system(paste("vcftools --vcf ", outnam, ".vcf --freq2 --positions out.kept.sites --out ", outnam, sep = ""))  
  system(paste("vcftools --vcf ", outnam, ".vcf --site-pi --positions out.kept.sites --out ", outnam, sep = ""))
  system(paste("vcftools --vcf ", outnam, ".vcf --geno-r2 --positions out.kept.sites --out ", outnam, " --min-r2 0.001 --ld-window-bp 100000", sep = ""))
  
  tajD <- read.table(paste(PATH,"/", outnam,".Tajima.D", sep = ""), header = T)
  het <- read.table(paste(PATH,"/", outnam,".het", sep = ""), header = T)
  freq2 <- read.table(paste(PATH,"/", outnam,".frq", sep = ""), header = T, row.names=NULL)
  pi <- read.table(paste(PATH,"/", outnam,".sites.pi", sep = ""), header = T, row.names=NULL)
  foo <- read.table(paste(PATH,"/", outnam,".geno.ld", sep = ""), header = T)
  dist <- (abs(foo$POS1 - foo$POS2))
  nloc <- length(dist)
  rsq <- foo$R.2
  n = foo$N_INDV[1]
  Cstart <- c(C=0.1)
  modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))),
                start=Cstart,
                control=nls.control(maxiter=1000))

  rho <- summary(modelC)$parameters[1]
  newrsq <- ((10+rho*dist)/((2+rho*dist)*(11+rho*dist)))*(1+((3+rho*dist)*(12+12*rho*dist+(rho*dist)^2))/(n*(2+rho*dist)*(11+rho*dist)))
  newfile <- data.frame(dist, newrsq)

  maxld <- max(newfile$newrsq) #using max LD value from adjusted data
  halfdecay = maxld*0.5 # https://jujumaan.com/2017/07/15/linkage-disequilibrium-decay-plot/
  halfdecaydist <- newfile$dist[which.min(abs(newfile$newrsq-halfdecay))]


  outtajD <- rbind(outtajD, cbind(i, tajD, outnam))
  outhet <- rbind(outhet, cbind(i, het, outnam))
  outfreq2 <- rbind(outfreq2, cbind(i, freq2, outnam))
  outsitespi <- rbind(outsitespi, cbind(i, pi, outnam))
  outld <- rbind(outld, cbind(i, halfdecay, halfdecaydist, nloc, outnam, rho))

}


write.table(outtajD, paste(PATH, "/compiled.tajimasD", sep = ""), quote = F, row.names = F)
write.table(outhet, paste(PATH, "/compiled.het", sep = ""), quote = F, row.names = F)
write.table(outfreq2, paste(PATH, "/compiled.freq2", sep = ""), quote = F, row.names = F)
write.table(outsitespi, paste(PATH, "/compiled.sites.pi", sep = ""), quote = F, row.names = F)
write.table(outld, paste(PATH, "/compiled.half.decay", sep = ""), quote = F, row.names = F)

