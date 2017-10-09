# Ancestry Mapping
PPCA & Variational Inference for Ancestry Mapping - CS229 Machine Learning in Bioinformatics Project

## Background

Inferring ancestry has long been recognized as a confounding factor in genetic testing. Estimation of individual ancestry from genetic data is useful in various applications such as during analysis of disease association studies, understanding human population history, interpreting personal genomic variation and also as covariates to correct for population stratification. The project from the human genetics community point of view focuses on determining and evaluating new methods of Variational Inference and Probabilistic Principal Components towards the estimation and visualization of admixture ancestral origin for an individual.

## Contents
+ Data - 1000 genome chr-22 data
+ Plots - graphs generated for project
+ Report & Presentation - final report & presentation
+ Scripts - PPCA & Variational Inference source code

## PPCA

In this project, we apply a different statistical method for inferring ancestry that avoids the painstaking computational complexity of PCA and further assumes an underlying probabilistic model of the genetic data. We purpose application of Probabilistic Principal Component Analysis (PPCA) (Tipping, M.E et al.) for the problem of ancestry inference. Unlike standard PCA methods that yields a multitude of SNPs with non-zero contributions to the principal components, PPCA method produces a limited number of influential SNPs for ancestry inferencing, with no loss in visualization of ancestry.

## Variational Inference

For the purposes of ancestry proportion estimation well-designed sampling schemes need to generate a large number of posterior samples to resolve convergence and mixing issues and yield accurate estimates of ancestry proportions, greatly increasing the time complexity of inference for large genotype data sets. To provide faster estimation we purpose a Variational Bayesian inference framework to speed up the inference of population structure.
