library(data.table)
# if(!exists("infer_population_structure", mode="function")) source("variational_inference.R")

source("variational_inference.R")
n_individuals = 100
n_snp = 50
K=5
threshold = 10^-6

genome = matrix(sample.int(3, n_individuals*n_snp, TRUE), n_individuals, n_snp)
genome = genome - 1

save.diploid = as.matrix(read.table("admixture-testing.geno", header = TRUE))

variational.params = infer_population_structure(save.diploid[,1:100], 5, threshold)

Q = variational.params$Admixture
P = variational.params$AlleleFreq
Q
P

