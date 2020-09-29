
###getcompiled.R###

#Compiles all like files for given simulation parameters. 
#Usage: Rscript getcompiled.R


files_het <- list.files("./", pattern="compiled.het", recursive=TRUE)
PATH <- getwd()

hetout <- NULL
for ( i in files_het){
	het <- read.table(i, header = TRUE)
	hetout <- rbind(hetout, het)
}

files_taj <- list.files("./", pattern="compiled.tajimasD", recursive=TRUE)

tajout <- NULL
for ( i in files_taj){
	taj <- read.table(i, header = TRUE)
	tajout <- rbind(tajout, taj)
}

files_freq2 <- list.files("./", pattern="compiled.freq2", recursive=TRUE)

freq2out <- NULL
for ( i in files_freq2){
	freq2 <- read.table(i, header = TRUE)
	freq2out <- rbind(freq2out, freq2)
}

sitespi_files <- list.files("./", pattern="compiled.sites.pi", recursive=TRUE)

piout <- NULL
for ( i in sitespi_files){
	pi <- read.table(i, header = TRUE)
	piout <- rbind(piout, pi)
}

ldout_files <- list.files("./", pattern="compiled.half.decay", recursive=TRUE)
ldout <- NULL
for ( i in ldout_files){
	ld <- read.table(i, header = TRUE)
	ldout <- rbind(ldout, ld)
}

write.table(hetout, paste(PATH, "/compiled.het", sep = ""), quote = F, row.names = F)
write.table(tajout, paste(PATH, "/compiled.tajimasD", sep = ""), quote = F, row.names = F)
write.table(freq2out, paste(PATH, "/compiled.freq2", sep = ""), quote = F, row.names = F)
write.table(piout, paste(PATH, "/compiled.sites.pi", sep = ""), quote = F, row.names = F)
write.table(ldout, paste(PATH, "/compiled.half.decay", sep = ""), quote = F, row.names = F)



