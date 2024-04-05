\name{GPMLtoIgraph}
\alias{GPMLtoIgraph}
\alias{nodeLegend}
\alias{edgeLegend}
\title{
Converting GPML FIles to Igraph Objects
}
\description{
  Extract entities of different type from GPML files in order to convert
  the pathway to a mathematical graph that we can compute on.
}
\usage{
GPMLtoIgraph(xmldoc, returnLists = FALSE, debug = FALSE)
nodeLegend(x, graph)
edgeLegend(x, graph)
}
\arguments{
  \item{xmldoc}{Either the name of an XML file meeting the
    specifications of the Genomic Pathway Markup Languag (GPML), or an
    object of class \code{XMLInternalDocument} obtained by running such
    a file through the \code{xmlParseDoc} function of the \code{XML}
    package. (All of the functions described here will call
    \code{xmlParseDoc} if it hasn't already been used.)}
  \item{returnLists}{A logical value string; should the retunr value
    include the node list and edge list matrices?}e
  \item{debug}{A locaigl value; should debugging progress information be
  printed? Probably best to leave it equal to FALSE.}
  \item{x}{A character string, such as "topleft" indicating where
    to place the legend.}
  \item{graph}{An \code{igraph} object as produced by the
    function\code{GPMLtoIgraph}.}
}
\value{
  The \code{GMLtoIgraph} function usually returns an \code{igraph}
  object that represents the pathway defined by the input
  \code{xmlfile}. If the argument \code{returnLists = TRUE}, then it
  returns a list containing three components; \code{graph} is the
  \code{igraph} object, \code{nodes} is a data frame containing node
  information where each row is a node, and \code{edges} is a matrix
  containing edge information where each row is an edge. The node and
  edge information can be used to reproduce the graph in any network or
  graph visualization tool that accepts such matrices to describe the
  graph.

  Both \code{nodeLegend} and \code{edgeLegend} invisibly return the same
  value that id returned by the \code{\link{legend}} function that is
  used in the implementation.
}
\details{
  \code{GPMLtoIgraph} is the main function of the
  \code{WayFindR}package. It achieves the primary goal of converting
  pathways from the GPML file format to a mathematical graph, in the
  format defined by the \code{\link{igraph}} package, at which point we
  can apply a wide variety of grph algorithms from computer science to
  "compute on biological pathways".

  The implementation of \code{GPMLtoIgraph} relies on the utility
  funcrtions described in \link{gpml-utility}. Briefly, the algorithm
  starts by collecting all nodes (\code{DataNodes} in GPML) and edges
  (\code{Interactions} in GPML) from the GPML input file. However, GPML
  includes two other structrues with (semantic) biological
  meaning. First, the GPML description includes the idea of an
  (invisible) "\code{Anchor}" that allows one edge to point to another
  edge. We expand those invisible target locations into full-fledged
  nodes in the final graph. Second, Second, GPML includes "\code{Group}"
  that represent protein complexes or sets of closely related genes.  In
  \code{WayFindR}, we represent such groups as threir own nodes in the
  final graph, and add "contained" edges linking in the group members.
  The transformations of Anchors and Groups do not change the
  fnudamaental topology (in particular, the existence of cycles) of the
  resulting graph.

  Finally, GPML includes non-semantic features (including "\code{Lables}"
  and "\code{Shapes}") that are (misued) by some pathway authors as te
  targets of edges. \code{WayFindR} converts any targeted non-semantic
  features into nodes in order to preserve as much information as
  possible from the original pathway in WikiPathways.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Polina Bombina \email{pbombina@augusta.edu}
}
\examples{
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
graf <- GPMLtoIgraph(xmlfile)
set.seed(13579)
L <- igraph::layout_with_graphopt
plot(graf, layout=L)
nodeLegend("topleft", graf)
edgeLegend("bottomright", graf)
}
\keyword{ graphs }
\keyword{ aplot }
