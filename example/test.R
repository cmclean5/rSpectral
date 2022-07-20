##reset
rm(list=ls())

library(rSpectral)
library(igraph)

## load graph
gg = read.graph("PPI_Presynaptic.gml",format="gml")

## store graph's edges as data.frame
el = as.data.frame(get.edgelist(gg,names=T))

logFile<-'rSpectral.log'
for(i in 1:500){

    rSpectral::load_data(df=el)

    cat(format(Sys.time(), "%b %d %X"),i,'Started.\n',file = logFile,append = TRUE)

    status = rSpectral::spectral(fix_neig=1)
    spec   = rSpectral::membership()
    rm(spec)

    cat(format(Sys.time(), "%b %d %X"),i,'Finished.\n',file = logFile,append = TRUE)
}
cat(format(Sys.time(), "%b %d %X"),i,'Done 500 times.\n',file = logFile,append = TRUE)
