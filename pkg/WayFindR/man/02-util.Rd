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
  The \code{collectEdges} function returns a data frame with three
  columns (Source, Target, and MIM), where each row describes one edge
  (or "Interaction" in the GPML terminology) of the pathway/graph. The
  \code{Source} and \code{Target} columns are the alphanumeric
  identifiers of items decribing nodes. The \code{MIM} column is the
  edge type in GPML, which often contains termiology based on the
  "Molecular Interaction Map" standard.

  The \code{collectNodes} function returns a data frame with three
  columns (GraphId, label, and Type), where each row describes node
  or vertex of the pathway/graph. The \code{GraphId} column is a unique
  alphanumeric identifier. The \code{label} column is a human-readbale
  name for the node, often the official gene symbol.

  The \code{collectAnchors} function returns a list containing a
  \code{nodes} element (in the same format returned by
  \code{collectNodes}) and an \code{edges} element (in the same format
  returned by \code{collectEdges}).

  The \code{collectGroups} function returns a list with \code{nodes} and
  \code{edges} components, just liek the one from \code{collectAnchors}.

  Both the \code{collectLabels} and \code{collectShapes} functions return
  the same kind of data frame that is retuedn by \code{collectNodes}.

}
\details{
  These functions are primarily intended as utility functions that
  implement processes required by the main function in the package,
  \code{link{GPMLtoIgraph}}. They have been exposed by being exported to
  the end user for use in debugging problematic GPML files or to reuse
  the GPML files in other contexts than the one we focus on in this
  package.

  While the meaning of  nodes (known as \code{DataNodes} in GPML) and
  edges (known as \code{Interactions} in GPML) should be obvious, some
  of the otrher objects are less so. An \code{Anchor} in GPML, for
  example, is an invisible object used to allow an edge to point to
  another edge instead of to a node. That structure isn't alloed in
  graphs in mathematics or computer science.  \code{WayFindR} handles
  this by creating a new node type to represent the anchor position,
  breaking the target edge into two pieces separated by the anchor, and
  adding an edge from the source of the anchored edge to the new node.

  In GPML, a \code{Label} is a text box allowing extra information to be
  placed on a pathway. A \code{Shape} is a graphical display object. The
  deefintion type document (dtd) for GPML describes both of these
  entities as non-semantic, intending them for display purposes
  only. However, some authors of pathways int he WikiPathways database
  use such objects as the (usually, final or "leaf") targets of
  interaction edges. When that happens, tghe \code{WayFindR} pacakge
  deals with it by creating actual nodes to represent such labels or
  shapes. Other labels and shapes are ignored.

  Th
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Polina Bombina \email{pbombina@augusta.edu}
}
\examples{
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
xmldoc <- XML::xmlParseDoc(xmlfile)
edges <- collectEdges(xmldoc)
nodes <- collectNodes(xmldoc)
anchors <- collectAnchors(xmldoc)
labels <- collectLabels(xmldoc)
edges <- collectShapes(xmldoc)
groups <- collectGroups(xmldoc, nodes)
}
\keyword{ manip }
