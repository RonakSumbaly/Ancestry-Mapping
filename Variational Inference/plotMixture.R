library(magrittr)
library(dplyr)
library(ggplot2)

Nth.delete<-function(dataframe, n) return(dataframe[-(seq(n,to=nrow(dataframe),by=n)),])

plotAdmixture <- function(sampleId, popId, k, value, colors, kOrder = 1:max(k), width = 1, alpha = 1, popLabels = TRUE, labColors = "black", rot = 0, showLegend = FALSE){
  
  ## prepare plotting dataframe
  d <- data.frame(sampleId = sampleId, popId = popId, k = k, value = value)
  idxP <- levels(popId)
  d$width <- width
  d$alpha <- factor(alpha, levels = sort(unique(alpha)))
  d$k <- factor(d$k, levels = kOrder)
  d <- d[order(d$sampleId, d$popId, d$k),]
  x <- d$width[!duplicated(d$sampleId)]
  d$xend <- rep(cumsum(x), each = max(k)) + 0.5 ## end x coordinate for each sample
  idx <- diff(c(1, x)) != 0 ## index of changes in width
  x[idx] <- x[idx] / 2 + 0.5
  x1 <- cumsum(x)
  d$x <- rep(x1, each = max(k)) ## center x coordinate for each sample
  
  ## get x & xend positions for each population
  d1 <- filter(d, k == 1) %>% group_by(popId) %>% summarise(x = mean(x), xend = max(xend))
  
  ## set up labels
  if(rot == 0){
    hjust <- 0.5
    vjust <- 1
  } else {
    hjust <- 1
    vjust <- 0.5
  }
  p1 <- geom_text(aes(x = x, y = -0.05, label = sampleId, fill = NULL), colour = "black", size = 2, angle = rot, hjust = hjust, vjust = vjust, data = d)
  if(popLabels){
    p1 <- geom_text(aes(x = x, y = -0.05, label = popId, fill = NULL), colour = labColors, size = 2, angle = rot, hjust = hjust, vjust = vjust, data = d1)
  }
  o <- theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, colour = NA), strip.text = element_text(colour = "white", face = "bold"), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.ticks.length = grid::unit(0, "cm"), plot.background = element_rect(fill = "transparent", colour = NA))
  o1 <- theme(legend.position = "top", legend.box = "horizontal")
  if(!showLegend){
    o1 <- theme(legend.position = "none")
  }
  
  ## prepare plot
  p <- ggplot(d, aes(x = x, y = value, fill = k, colour = k))
  p + geom_bar(aes(width = width, alpha = alpha), stat = "identity", size = 0.05) + p1 + geom_segment(aes(x = xend, xend = xend, y = 0, yend = 1, fill = NA), colour = "black", size = 0.1, data = d1) + scale_colour_manual(values = colors) + scale_fill_manual(values = colors) + scale_alpha_manual(values = unique(sort(alpha)), guide = FALSE) + scale_y_continuous(breaks = NULL, limits = c(-0.5, 1.01)) + theme_bw() + labs(x = "", y = "") + o + o1 + geom_rect(aes(xmin = 0.5, ymin = 0, xmax = max(xend), ymax = 1), fill = NA, colour = "black", size = 0.1)
  ggsave("admix.png",width = 8, height = 2.5)
}

## read full panel results
qMatrix <- matrix(scan("final.5.Q"), ncol = 5, byrow = TRUE)

individuals =save.diploid = as.matrix(read.table("chr-22.ind", header = FALSE))
individuals = data.frame(individuals)
individuals = Nth.delete(individuals, 2)
colnames(individuals) <- c("sampleId", "V2", "popId")

rownames(qMatrix) <- individuals$sampleId

## plot
d <- reshape2::melt(qMatrix)
colnames(d) <- c("sampleId", "k", "value")

d$sampleId <- factor(d$sampleId, levels = unique(d$sampleId))
d$popId <- individuals$popId[match(d$sampleId, individuals$sampleId)]
d$popId <- factor(d$popId, levels = unique(d$popId))

# idxCLabel <- popInfo$color[match(levels(d$popId), popInfo$popId)]
idxCCluster <- c("violetred3", "wheat3", "steelblue3","paleturquoise3", "gold2")
# pdf("trial.pdf", width = 8, height = 2.5)
plotAdmixture(sampleId = d$sampleId, popId = d$popId, k = d$k, value = d$value, colors = idxCCluster, alpha = 1, width = 1, showLegend = TRUE, rot = 90)
# dev.off()