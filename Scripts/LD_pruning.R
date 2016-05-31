############################################################################
# Project Title: Ancestry Mapping
# Done By: Ronak Sumbaly
# Description: perform LD pruning 
############################################################################

library(trio)

threshold = 0.7

# function to plot linkage disequilibrium between SNPs
LD.plot = function(ld.mat, color = heat.colors(50), title="Linkage Disequilibrium"){
  break.points = 0:length(color)/length(color)
  mar.orig = (par.orig = par(c("mar", "las", "mfrow")))$mar
  oma = c(2, 1, 1, 2)
  def.par = par(no.readonly = TRUE)
  m = matrix(c(1, 1, 1, 2), 2, 2)
  h = (3 + mar.orig[2]) * par("csi") * 2.54
  layout(m, heights = c(1, lcm(0.5 * h)))
  image(1-ld.mat, axes=FALSE, col=color, cex.main = 1, breaks=break.points, main=title)
  box()
  par(mar = c(2, 1, 1.5, 2))
  a = matrix(1:length(color), ncol = 1)
  image(a, axes = FALSE, col = color[length(color):1])
  box(bty = "o")
  plot.label <- as.character(seq(from = 0, to = 1, length = 11))
  axis(side = 1, at = seq(from = 0 - 1/11/2, to = 1 + 1/11/2, length = 11), labels = plot.label)
  par(def.par)
}

# get LD for window size of 1000 SNPs
start = 1 # start position
end = 1000 # end position
correlated.snps = c()
while (end <= dim(diploid)[2] & start < end) {
  linkage.disequilibrium = getLD(diploid[,start:end],asMatrix = TRUE, which = "rSquare")$rSquare
  cat(start, " ", end, "\n")
  remove.snps = which(linkage.disequilibrium > threshold, arr.ind = TRUE)
  correlated.snps = c(correlated.snps, unique(colnames(linkage.disequilibrium)[remove.snps[,2]]))
  
  start = start + 1000
  end = end + 1000
  if (end > dim(diploid)[2]) end = dim(diploid)[2]
}
correlated.snps = unique(correlated.snps)

# plot LD of last SNP set
LD.plot(linkage.disequilibrium)

pruned.data = subset.data.frame(diploid, select = -(which(colnames(diploid) %in% correlate.snp2)))
update.snp.count = dim(pruned.data)[2]

cat("Number of SNPs pruned : ", length(correlated.snps))

rm(start)
rm(end)