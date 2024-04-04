\name{GPMLtoIgraph}
\alias{GPMLtoIgraph}
\title{
Converting GPML FIles to Igraph Objects
}
\description{
  Extract entities of different type from GPML files in order to convert
  the pathway to a mathematical graph that we can compute on.
}
\usage{
GPMLtoIgraph(xmldoc)
}
\arguments{
  \item{xmldoc}{Either the name of an XML file meeting the
    specifications of the Genomic Pathway Markup Languag (GPML), or an
    object of class \code{XMLInternalDocument} obtained by running such
    a file through the \code{xmlParseDoc} function of the \code{XML}
    package. (All of the functions described here will call
    \code{xmlParseDoc} if it hasn't already been used.)}
}
\value{
  Pending
}
\details{
  Pending
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Polina Bombina \email{pbombina@augusta.edu}
}
\examples{
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
graf <- GPMLtoIgraph(xmlfile)
}
\keyword{ manip }
