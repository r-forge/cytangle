\name{gpml-utility}
\alias{collectEdges}
\alias{collectNodes}
\alias{collectGroups}
\alias{collectAnchors}
\alias{collectLabels}
\alias{collectShapes}
\title{
Utility Functions to Parse GPML FIles
}
\description{
}
\usage{
collectEdges(xmldoc)
collectNodes(xmldoc)
collectGroups(xmldoc, allnodes)
collectAnchors(xmldoc)
collectLabels(xmldoc)
collectShapes(xmldoc)
}
\arguments{
  \item{xmldoc}{Either the name of an XML file meeting the
    specifications of the Genomic pathway Markup Languag (GPML), or an
    object of class \code{XMLInternalDocument} obtained by running such
    a file through the \code{xmlParseDoc} funciton of te \code{XML}
    package. (All of the functions described here will call
    \code{xmlParseDoc} if it hasn't already been used.)}
  \item{allnodes}{A data frame containng node information, in the format
    produced by the \code{collectNodes} function.}
}
\value{
  The \code{collectNodes} function returns a data frame with three
  columns (GraphId, label, and Type), where each row consusts of a node
  or vertex of the pathway/graph.

  
}
\details{
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Polina Bombina \email{pbombina@augusta.edu}
}
\examples{
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
edges <- collectEdges(xmlfile)
nodes <- collectNodes(xmlfile)
anchors <- collectAnchors(xmlfile)
labels <- collectLabels(xmlfile)
edges <- collectShapes(xmlfile)
groups <- collectGroups(xmlfile, nodes)
}
\keyword{ manip }
