\name{as.graphNEL}
\alias{as.graphNEL}
\title{
Convert igraph object to Ragraph
}
\description{
  Converts an igraph object into an Ragraph object, preserving as many
  attributes of nodes and edges as possible.
}
\usage{
as.graphNEL(G)
}
\arguments{
  \item{G}{An \code{igraph} representation of a graph.}
}
\details{
  The only tricky part of the implementation is to transfer as many of
  the properties of nodes and edges from one representation to another,
  since the \code{igraph} package already has a conversion function
  called \code{\link[igraph]{as_graphnel}}. That version keeps tge ege
  list and node list, but not the attributes. Our implementation uses
  the \code{\link[igraph]{vertex_attr}} and
  \code{\link[igraph]{edge_attr}} functions from \code{igraph} to
  extract the properties and then passes them to the
  \code{\link[Rgraphviz]{agopen}} function in \code{Rgraphviz} to map
  them to the slightly different naes use din te output.

  The function should probably be named \code{as.Ragraph}, but we didn't
  understand enough about Rgraphviz when we started.

  Also, the converson will fail if the igrph input has any multiple edges,
  which are not supported in graphNEL objects.
}
\value{
  An \code{Ragraph} representation of the same graph.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Polina Bombina \email{pbombina@augusta.edu}
}
\examples{
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
G <- GPMLtoIgraph(xmlfile)
GN <- as.graphNEL(G)
}
\keyword{ utility }
\keyword{ manip }
\keyword{ graphs }
