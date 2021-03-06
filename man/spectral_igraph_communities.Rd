% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/igraph.R
\name{spectral_igraph_communities}
\alias{spectral_igraph_communities}
\title{Spectral clustering for \code{igraph} objects}
\usage{
spectral_igraph_communities(
  g,
  Cn_min = 1L,
  tol = 1e-05,
  names = 1L,
  fix_neig = 0L
)
}
\arguments{
\item{g}{\code{igraph} object}

\item{Cn_min}{minimum cluster size}

\item{tol}{tolerance}

\item{names}{are we dealing with alphaNumeric (1) or numeric (!1) ids}

\item{fix_neig}{wether fixing neighbouring nodes found in same community}
}
\value{
\code{\link[igraph]{communities}} object
}
\description{
Spectral clustering for \code{igraph} objects
}
\examples{
data(karate,package='igraphdata')
c<-spectral_igraph_communities(karate)
}
