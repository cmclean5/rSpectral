##reset
rm(list=ls())

library(rSpectral)
library(igraph)

## load graph
gg = read.graph("PPI_Presynaptic.gml",format="gml")

## store graph's edges as data.frame
el = as.data.frame(get.edgelist(gg,names=T))

cat('load data...')
rSpectral:::load_data(df=el)
cat('done.\n')

cat('run spectral clustering - basic operation...\n')

## run spectral clustering on graph
ptm <- proc.time()
status = rSpectral:::spectral(fix_neig=0)
spec   = rSpectral:::membership(detach_graph=0)
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
status = rSpectral:::spectral(fix_neig=1)
spec   = rSpectral:::membership(detach_graph=0)
pet <- proc.time() - ptm
cat("spectral modularity \n")
cat(sprintf("time = %.3f \n", pet[[1]]))

cat('done.\n')


## delete model
rm(spec)

cat('run spectral clustering - fixing neighbouring nodes, Cnmin=5...\n')

## run spectral clustering on graph, fixing neighbouring nodes found
## in same community.
ptm <- proc.time()
status = rSpectral:::spectral(fix_neig=1,Cn_min=5)
spec   = rSpectral:::membership(detach_graph=1)
pet <- proc.time() - ptm
cat("spectral modularity \n")
cat(sprintf("time = %.3f \n", pet[[1]]))

cat('done.\n')

#cat('free model...')
#rSpectral::freeSpace()
#cat('done.\n')
