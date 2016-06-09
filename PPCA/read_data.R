############################################################################
# Project Title: Ancestry Mapping
# Done By: Ronak Sumbaly
# Description: read 1000 genome data 
############################################################################

library(data.table)
library(stringr)
library(flashpcaR)

setwd("~/Documents/Ancestry Mapping")

##### INITIALIZING DATA #####

# read in 1000 genome data
haploid = t(fread("1000genomes/chr-22.geno.reduced.csv", header = FALSE, data.table = TRUE))
individuals = fread("1000genomes/chr-22.ind", header = FALSE,  data.table = TRUE)
snps = fread("1000genomes/chr-22.snp", header = FALSE, data.table = TRUE)

individuals.count = dim(individuals)[1] # haploid chromosomes
snps.count = dim(haploid)[2] # number of snps 
mapping = individuals$V3[seq(1, 2184, 2)] # location 
mapping = unlist(lapply(mapping, function(k){strsplit(k, ":")[[1]][2]}))
numeric.mapping = as.numeric(factor(mapping)) # numeric locations

# construct diploid data
diploid = haploid[seq(1,individuals.count,2),] + haploid[seq(2,individuals.count,2),] # add consecutive haploid rows 

# add row names and col names to diploid data
rownames(diploid) = individuals$V1[seq(1,individuals.count,2)]
colnames(diploid) = paste("snp", str_split_fixed(snps$V1, ":", 2)[c(1:snps.count),2],sep = "")

# update individual count
individuals.count = dim(diploid)[1] # diploid chromosomes

rm(haploid)