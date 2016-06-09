############################################################################
# Project Title: Ancestry Mapping
# Done By: Ronak Sumbaly
# Description: perform pca and ppca
############################################################################

library(pcaMethods)
library(flashpcaR)

# pca mapping
pca.baseline = flashpca(pruned.data, do_loadings=TRUE, verbose=TRUE, stand="binom", ndim=3,
                        nextra=100)
pca.baseline = prcomp(diploid)

# plot pca with color coding
pairs(pca.baseline$x[,1:3], col = numeric.mapping, oma=c(4,4,6,12), labels = c("PCA1","PCA2","PCA3"), pch = 20)
par(xpd=FALSE)
legend(0.85, 0.6, as.vector(c("EU", "AM", "AS", "AF")),  
       fill=c("red", "green3", "blue", "black"))


plot((pca.baseline$sdev^2/sum(pca.baseline$sdev^2))[1:10],type = "o",xlab = "Principal Components",ylab = "Variance Explained", main = "Variance Associated")
plot(pca.baseline, type = "l")


individuals

pairs(pca.baseline$x[,1:3], col = numeric.mapping, labels = c("PCA1","PCA2","PCA3"), pch= 20)

# ppca mapping
ppca.baseline = ppca(pruned.data, nPcs = 3)

# plot ppca with color coding
pairs(ppca.baseline@scores, col = numeric.mapping, oma=c(4,4,6,12), labels = c("PPCA1","PPCA2","PPCA3"), pch = 20)
par(xpd=FALSE)
legend(0.85, 0.6, as.vector(c("EU", "AM", "AS", "AF")),  
       fill=c("red", "green3", "blue", "black"))

plot((ppca.baseline@sdev^2/sum(ppca.baseline@sdev^2))[1:10],type = "o",xlab = "Principal Components",ylab = "Variance Explained", main = "Variance Associated")
