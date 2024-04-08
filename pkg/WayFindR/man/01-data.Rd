\name{gpml-data}
\alias{gpml-data}
\alias{edgeColors}
\alias{edgeTypes}
\alias{nodeColors}
\alias{nodeShapes}
\docType{data}
\title{GPML GraphingR Data}
\description{
  Our standard colors, line types (for edges), and shapes (for nodes)
  to display graph features from WikiPathways GPML files.
}
\usage{
data(edgeColors)
data(edgeTypes)
data(nodeColors)
data(nodeShapes)

}
\format{There are four different objects.
  \describe{
    \item{\code{edgeColors}}{A named character vector defining the
      colors used to display 17 different kinds of edges.}
    \item{\code{edgeTypes}}{A named character vector defining the
      line types (solid, dashed, etc.) used to display 17 different
      kinds of edges.}
    \item{\code{nodeColors}}{A named character vector defining the
      colors used to display 16 different kinds of nodes.}
    \item{\code{nodeShapes}}{A named character vector defining the
      symbols used to display 16 different kinds of nodes.}
  }
}
\source{
  We downloaded the complete set of 889 human (Homo sapiens) pathways
  from WikiPathways, in GPML format, at the end of March 2024. We wrote
  and ran perl scripts to iterate over this set of pathways and extract
  all the node types and all the edge types used anywhere in the set. We
  added a handful of node types (Undefined, Shape, Label, EDGE) in order
  to enable almost all of the pathways to be converted into mathematical
  graphs. We then defined colors and shapes to be used to distinguish
  different types in plots.

  The situation with edges was slightly more complicated. The pathway
  editor (PathVisio; https://pathvisio.org/) typically used to create
  pathways in WikiPathways supports two different biologically
  meaningful vocabularies for edges: Molecular Interaction Maps (MIM;
  https://discover.nci.nih.gov/mim/formal_mim_spec.pdf) and System
  Biology Graphical Notation (SBGN; https://sbgn.github.io/), along with
  a very simple descrption as "Arrows" or "TBars". To make things more
  consistent, we converted both the SBGN and simple systems to match the
  MIM specification. We then defined colors and edge types so we could
  distinguish different edge types in plots.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Polina Bombina \email{pbombina@augusta.edu}
}
\examples{
data(edgeColors)
data(edgeTypes)
data(nodeColors)
data(nodeShapes)
if (require(Polychrome)) {
  swatch(edgeColors)
  swatch(nodeColors)
}
plot(0,0, type = "n", xlab="", ylab = "")
legend("center", legend = names(edgeColors),  lwd = 3,
       col = edgeColors,  lty = edgeTypes)
num <- c(rectangle = 15, circle = 16)
plot(0,0, type = "n", xlab="", ylab = "")
legend("center", legend = names(nodeColors),  cex = 1.5,
       col = nodeColors,  pch = num[nodeShapes])
}
\keyword{datasets}
\keyword{color}

