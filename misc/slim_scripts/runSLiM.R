#!/usr/bin/evn Rscript
### runSLiMthenLD.R ###

# This takes a given .slim file then performs the simulation n number of times and produces a vcf file for each simulation sequentially and the vcfile will be over written following each simulation replicate to save space. Following each simulation the output the vcf file will be placed in directory "vcf_files" created by this script. 
# SLiM version: version 3.2.1, built Jan 29 2019 18:20:03
# vcftools version: 0.1.15

# Usage: Rscript runSLiM.R <slim file> <rep_num> <vcf_num>
# User needs to edit path in .slim as needed. All paths are relative. User needs to edit child vcf file names. Example here for 90% bottleneck. 

args <- commandArgs(TRUE)

if (length(args)!=3) {
  stop("Three arguments must be supplied", call.=FALSE)
}

setwd(".")
PATH <- getwd()
slim <- args[1]
rep_num <- as.numeric(args[2])
harvest_num <- as.numeric(args[3])
if (harvest_num>rep_num){
  stop("Number of sumulations to harvest vcf files cannot be larger than rep number", call.=FALSE)
}
outnam <- strsplit(slim, "[.]")[[1]][1]
vcf_save <- sample(rep_num, harvest_num)

for (i in 1:rep_num) {
  system(paste("slim ", slim, sep = "")) # envoke SLiM, edit file path...
  system(paste("mv bn_90_pre.vcf ", "bn_90_pre_rep", i, ".vcf", sep = ""))  
  system(paste("mv bn_90_post_01.vcf ", "bn_90_post_01_rep", i, ".vcf", sep = ""))
  system(paste("mv bn_90_post_10.vcf ", "bn_90_post_10_rep", i, ".vcf", sep = ""))
  system(paste("mv bn_90_post_25.vcf ", "bn_90_post_25_rep", i, ".vcf", sep = ""))  
  system(paste("mv bn_90_post_50.vcf ", "bn_90_post_50_rep", i, ".vcf", sep = ""))
  system(paste("mv bn_90_post_100.vcf ", "bn_90_post_100_rep", i, ".vcf", sep = ""))
  system(paste("mv bn_90_post_200.vcf ", "bn_90_post_200_rep", i, ".vcf", sep = ""))
  }
system("mkdir vcf_files")
system("mv *.vcf vcf_files")






