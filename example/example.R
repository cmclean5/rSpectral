##reset
rm(list=ls())

library(rSpectral)
library(igraph)

## load graph
gg = read.graph("PPI_Presynaptic.gml",format="gml")

## store graph's edges as data.frame
el = as.data.frame(get.edgelist(gg,names=T))

cat('run spectral clustering - normal...\n')

## run spectral clustering on graph
ptm <- proc.time()
spec = rSpectral::spectral(DF=el,fixNeig=0)
pet <- proc.time() - ptm
cat("spectral modularity \n")
cat(sprintf("time = %.3f \n", pet[[1]]))

cat('done.\n')

## delete model
rm(spec)

cat('run spectral clustering - fixing neighbouring nodes...\n')

## run spectral clustering on graph, fixing neighbouring nodes found
## in same community.
ptm <- proc.time()
spec = rSpectral::spectral(DF=el,fixNeig=1)
pet <- proc.time() - ptm
cat("spectral modularity \n")
cat(sprintf("time = %.3f \n", pet[[1]]))

cat('done.\n')
