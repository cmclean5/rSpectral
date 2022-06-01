##reset
rm(list=ls())

library(rSpectral)
library(igraph)

## load graph
gg = read.graph("PPI_Presynaptic.gml",format="gml")

## store graph's edges as data.frame
el = as.data.frame(get.edgelist(gg,names=T))

## run spectral clustering on graph
ptm <- proc.time()
spec = rSpectral::spectral(DF=el)
pet <- proc.time() - ptm
cat("spectral modularity \n")
cat(sprintf("time = %.3f \n", pet[[1]]))
