\name{gpml-utility}
\alias{gpml-utility}
\alias{collectEdges}
\alias{collectNodes}
\alias{collectAnchors}
\alias{collectLabels}
\alias{collectShapes}
\alias{collectGroups}
\title{
Utility Functions to Parse GPML FIles
}
\description{
  Extract entities of different types from GPML files in order to convert
  the pathway to a mathematical graph that we can compute on.
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
    specifications of the Genomic Pathway Markup Language (GPML), or an
    object of class \code{XMLInternalDocument} obtained by running such
    a file through the \code{\link{xmlParseDoc}} function of the \code{XML}
    package. (All of the functions described here will call
    \code{xmlParseDoc} if it hasn't already been used.)}
  \item{allnodes}{A data frame containing node information, in the format
    produced by the \code{collectNodes} function.}
}
\details{
  These functions are primarily intended as utility functions that
  implement processes required by the main function in the package,
  \code{\link{GPMLtoIgraph}}. They have been made accessible to
  the end user for use in debugging problematic GPML files or to reuse
  the GPML files in contexts other than the one we focus on in this
  package.

  While the meaning of  nodes (known as \code{DataNodes} in GPML) and
  edges (known as \code{Interactions} in GPML) should be obvious, some
  of the other objects are less so. For example, an \code{Anchor} in
  GPML is an invisible object used to allow an edge to point to 
  another edge instead of to a node. That structure isn't allowed in
  graphs in mathematics or computer science.  \code{WayFindR} handles
  this by creating a new node type to represent the anchor position,
  breaking the target edge into two pieces separated by the anchor, and
  adding an edge from the source of the anchored edge to the new node.

  In GPML, a \code{Label} is a text box allowing extra information to be
  placed on a pathway, and a \code{Shape} is a graphical display object. The
  definition type document (DTD) for GPML describes both of these
  entities as non-semantic, intending them for display purposes
  only. However, some authors of pathways in the WikiPathways database
  use such objects as the (usually, final or "leaf") targets of
  interaction edges. When that happens, the \code{WayFindR} package
  deals with it by creating actual nodes to represent such labels or
  shapes. Other labels and shapes are ignored.

  GPML also uses the idea of a \code{Group} as a first class object in
  their DTD. These are defined as  "A collection of structurally or
  functionally similar or related pathway elements." The GPML file
  subclassifies some groups as "Complexes", indicating that they
  represent physical interactions and bindings between two or more
  molecules. Other groups may simply indicate that there is a related set
  of molecules (for example, STAT2 and STA3) that play the same role at
  this point in the pathway. \code{WayFindR} deals with groups by
  creating a new node to represent the group as a whole and expanding
  the component genes into nodes with a single "contained" edge that
  points to the new group node.
}
\value{
  The \code{collectEdges} function returns a data frame with three
  columns (\code{Source}, \code{Target}, and \code{MIM}), where each row
  describes one edge (or "Interaction" in the GPML terminology) of the
  pathway/graph. The \code{Source} and \code{Target} columns are the
  alphanumeric identifiers of items decribing nodes. The \code{MIM}
  column is the edge type in GPML, which often contains terminology based
  on the "Molecular Interaction Map" standard. When creating an
  \code{\link{igraph}} object from a pathway, the first two columns are
  used as identifiers to define the nodes at each end of the edge.

  The \code{collectNodes} function returns a data frame with three
  columns (GraphId, label, and Type), where each row describes node
  or vertex of the pathway/graph. The \code{GraphId} column is a unique
  alphanumeric identifier. The \code{label} column is a human-readable
  name for the node, often the official gene symbol. When creating an
  \code{igraph} object from a pathway, the first column is used as
  identifier to define the node. Also, the \code{plot} method for
  \code{igraph}s recognizes the term \code{label} as a column that
  defines the text that should be displayed in a node.

  The \code{collectAnchors} function returns a list containing a
  \code{nodes} element (in the same format returned by
  \code{collectNodes}) and an \code{edges} element (in the same format
  returned by \code{collectEdges}). The \code{collectGroups} function
  returns a list with \code{nodes} and \code{edges} components, just
  like the one from \code{collectAnchors}. 

  Both the \code{collectLabels} and \code{collectShapes} functions return
  the same kind of data frame that is returned by \code{collectNodes}.
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
\keyword{ utility }
\keyword{ manip }
\keyword{ graphs }
