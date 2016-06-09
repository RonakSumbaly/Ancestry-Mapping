############################################################################
# Project Title: Ancestry Mapping
# Done By: Ronak Sumbaly
# Description: plot admixture proportions
############################################################################

library(reshape2)
library(ggplot2)

plot.admixture <- function(q.matrix) {
  
  q.matrix = data.frame(q.matrix)
  q.matrix = data.frame(cbind(matrix(1:dim(q.matrix)[1], ncol = 1), q.matrix))
  colnames(q.matrix) = c("ID", paste("AP", 1:(dim(q.matrix)[2] - 1), sep = ""))
  
  q.matrix.reformed <- melt(q.matrix, id.var="ID")
  
  ggplot(q.matrix.reformed, aes(x = ID, y = value, fill = variable)) + geom_bar(stat = "identity")

}