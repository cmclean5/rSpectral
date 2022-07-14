#' Spectral clustering for \code{igraph} objects
#'
#' @param g \code{igraph} object
#' @inheritParams spectral
#' @return \code{data.frame} with node names and membership information
#' @export
#' @importFrom igraph get.edgelist
#'
#' @examples
#' data(karate,package='igraphdata')
#' df.mem<-spectral_igraph_membership(karate)
spectral_igraph_membership<-function(g,Cn_min = 1L, tol = 0.00001, names = 1L, fix_neig = 0L){
  if(!inherits(g,'igraph')){
    stop('Graph should be "igraph" object.')
  }
  el = as.data.frame(igraph::get.edgelist(gg,names=T))
  rSpectral::load_data(df=df)
  status = rSpectral::spectral(fix_neig=1)
  spec   = rSpectral::membership(detach_graph=1)
  spec.df<- data.frame(names=spec$ID,membership=spec$K)
  return(spec.df)
}

#' Spectral clustering for \code{igraph} objects
#'
#' @param g \code{igraph} object
#' @inheritParams spectral
#' @return \code{communitios} object
#' @export
#' @importFrom igraph modularity vcount
#'
#' @examples
#' data(karate,package='igraphdata')
#' c<-spectral_igraph_communities(karate)
spectral_igraph_communities<-function(g,Cn_min = 1L, tol = 0.00001, names = 1L, fix_neig = 0L){
  if(!inherits(g,'igraph')){
    stop('Graph should be "igraph" object.')
  }
  df.mem<-spectral_igraph_membership(g)
  res<-list()
  res$vcount <- igraph::vcount(g)
  res$algorithm <- "spectral"
  res$membership <- df.mem$membership
  res$modularity <- igraph::modularity(karate,c$membership)
  res$names <- df.mem$names
  class(res) <- "communities"
}