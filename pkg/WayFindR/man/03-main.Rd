\name{XMLtoIgraph}
\alias{GPMLtoIgraph}
\alias{KGMLtoIgraph}
\alias{nodeLegend}
\alias{edgeLegend}
\title{
Converting GPML FIles to Igraph Objects
}
\description{
  Takes an XML file, either  GPML from WikiPathways, or KGML from KEGG,
  extracts the entities therein, 
  and makes minor adjustments necessary to convert it into an
  \code{\link[igraph]{igraph}} object. Along the way, it assigns a consistent set
  of colors, line types, and shapes. 
}
\usage{
GPMLtoIgraph(xmldoc, returnLists = FALSE, debug = FALSE)
KGMLtoIgraph(xmldoc, returnLists = FALSE, debug = FALSE)
nodeLegend(x, graph)
edgeLegend(x, graph)
}
\arguments{
  \item{xmldoc}{Either the name of an XML file meeting the
    specifications of the appropriate markup language (GPML or KGML), or
    an object of class \code{XMLInternalDocument} obtained by running such
    a file through the \code{xmlParseDoc} function of the \code{XML}
    package.}
  \item{returnLists}{A logical value; should the return value
    include the node list and edge list matrices?}
  \item{debug}{A logical value; should debugging progress information be
  printed? Probably best to leave it equal to FALSE.}
  \item{x}{A character string, such as "topleft" indicating where
    to place the legend.}
  \item{graph}{An \code{igraph} object as produced by the
    function\code{GPMLtoIgraph}.}
}
\details{
  \code{GPMLtoIgraph} and \code{KGMLtoIgraph} are the main functions of
  the \code{WayFindR}package. They achieve the primary goal of
  converting pathways from one of the biological graph XML file formats
  into a mathematical graph, in the format defined by the
  \code{\link[igraph]{igraph}} package. At that point, we can apply a wide
  variety of graph algorithms from computer science in  order to
  "compute on biological pathways".

  The implementation of these functions relies on the utility functions
  described in \code{\link{gpml-utility}} or \code{\link{kgml-utility}}.

  Briefly, the first algorithm starts by collecting all nodes
  (\code{DataNodes} in GPML) and edges (\code{Interactions} in GPML)
  from the GPML input file. However, GPML includes two other structures
  with (semantic) biological meaning. First, the GPML description
  includes the idea of an (invisible) "\code{Anchor}" that allows one
  edge to point to another edge. We expand those invisible target
  locations into full-fledged nodes in the final graph. Second, GPML
  includes "\code{Group}s" that represent protein complexes or sets of
  closely related genes.  In \code{WayFindR}, we represent such groups
  as their own nodes in the final graph, and add "contained" edges
  linking in the group members. The transformations of Anchors and
  Groups do not change the fundamental topology (in particular, the
  existence of cycles) of the resulting graph.

  Further, GPML includes non-semantic features (including "\code{Labels}"
  and "\code{Shapes}") that are (mis)used by some pathway authors as the
  targets of edges. \code{WayFindR} converts any targeted non-semantic
  features into nodes in order to preserve as much information as
  possible from the original pathway in WikiPathways.

  The KGML algorithm is similar in structure, burt has to deal with the
  diferent underlying sdtructure of th KGML specification. Tghei files
  contina three kinds of entities: \code{Entry}, \code{Relation}, and
  \code{Reaction}. An \code{Entry} becoemes a vertex. It can be a gene,
  a map (a link to another pathway), a group (as above, except that the
  members of the group are stored as an entity called a
  \code{Component} within the group Entry), an ortholog (a KEGG-defined
  set of genes that are teh "same" across species), or a compound
  (subdivided into compounds, glycans, or drugs, all of which we view as
  "SmallMolecules" analogous to what GPML calls a metabolite). A
  \code{Relation} is an edge that usually connects genes, but we must
  map the terminology annotating edge types into the MIM   space
  defining biological edges. Finally, a \code{Reaction} is an edge
  between compounds, which has no real analog in the WikiPathways
  universe. The only "type" associated with a reaction is whether it is
  "reversible" or "irreversible".
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
  graph. The \code{nodes} data frame includes columns for \code{color}
  and \code{shape}, and the \code{edges} data frame includes columns for
  \code{color} and \code{lty} that are recognized and used by the
  \code{\link[igraph]{plot.igraph}} function.
  
  Both \code{nodeLegend} and \code{edgeLegend} invisibly return the same
  value that is returned by the \code{\link[graphics]{legend}} function that is
  used in the implementation.
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
