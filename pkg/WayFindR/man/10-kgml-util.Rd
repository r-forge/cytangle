\name{kgml-utility}
\alias{kgml-utility}
\alias{collectEntries}
\alias{collectRelations}
\alias{collectReactions}
\title{
Utility Functions to Parse KGML FIles
}
\description{
  Extract entities of different types from KGML files in order to convert
  the pathway to a mathematical graph that we can compute on.
}
\usage{
collectEntries(xmldoc)
collectRelations(xmldoc)
collectReactions(xmldoc)
}
\arguments{
  \item{xmldoc}{Either the name of an XML file meeting the
    specifications of the KEGG Genomic Markup Language (KGML), or an
    object of class \code{XMLInternalDocument} obtained by running such
    a file through the \code{\link{xmlParseDoc}} function of the \code{XML}
    package. (All of the functions described here will call
    \code{xmlParseDoc} if it hasn't already been used.)}
}
\details{
  These functions are primarily intended as utility functions that
  implement processes required by the main function in the package,
  \code{KGMLtoIgraph}. They have been made accessible to
  the end user for use in debugging problematic KGML files or to reuse
  the KGML files in contexts other than the one we focus on in this
  package.

}
\value{
  The \code{collectReactions} and \code{collectRelations} functions
  return a data frame with three columns (\code{Source}, \code{Target},
  and \code{MIM}), where each row describes one edge of the
  pathway/graph. In KEGG, they distingiuiish between relations (which
  usually connect genes) and reactions (which connect chemical
  compounds). The \code{Source} and \code{Target} columns are the 
  alphanumeric identifiers of items decribing nodes. The \code{MIM}
  column is the edge type in KGML.

  The \code{collectEntries} function returns a data frame with three
  columns (GraphId, label, and Type), where each row describes one node
  or vertex of the pathway/graph. The \code{GraphId} column is a unique
  alphanumeric identifier. The \code{label} column is a human-readable
  name for the node, often the official gene symbol. When creating an
  \code{igraph} object from a pathway, the first column is used as an
  identifier to define the node. Also, the \code{plot} method for
  \code{igraph}s recognizes the term \code{label} as a column that
  defines the text that should be displayed in a node.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Polina Bombina \email{pbombina@augusta.edu}
}
\examples{
xmlfile <- system.file("pathways/WP3850.kgml", package = "WayFindR")
xmldoc <- XML::xmlParseDoc(xmlfile)
}
\keyword{ utility }
\keyword{ manip }
\keyword{ graphs }
