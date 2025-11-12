\name{gpml-cycles}
\alias{gpml-cycles}
\alias{findCycles}
\alias{discoverCycles}
\alias{interpretCycle}
\alias{cycleSubgraph}
\title{
Finding Cycles in Biological Pathways
}
\description{
  Tools to find and interpret cycles in graphs derived from pathways in
  WikiPathways.
}
\usage{
findCycles(graph)
discoverCycles(graph, nsteps = 50, futility = 20)
interpretCycle(v, graph)
cycleSubgraph(graph, cycles)
}
\arguments{
  \item{graph}{An \code{igraph} object produced by
    \code{\link{GPMLtoIgraph}} from a pathway file in GPML format.}
  \item{v}{One of the cycles produced by the \code{findCycles} function.}
  \item{cycles}{The list of cycles produced by the \code{findCycles}
    function.}
  \item{nsteps}{Length of each random walk}
  \item{futility}{After this number of random walks without seeing a new
    cycle, halt.}
}
\details{
  The implementation of the \code{findCycles} function, while provably
  correct, makes no concessions to efficiency. It is likely to be
  adequate for the vast majority of biological pathways present in
  WikiPathways, but may prove incapable of handling large complex graphs
  in general. The algorithm is a sraightforward double loop. It first
  iterates over all nodes, and uses the \code{neighbors} function from
  \code{igraph} to find all nodes that are directly connected to the
  start node. It then iterates over those neighbors and uses the
  \code{all_simple_paths} function from \code{igraph} to find all
  paths that lead back to the starting node. One aspect of its lack of
  efficiency is that each cycle with N nodes is found N times, once for
  each place you could start traversing the cycle. Before returning the
  value to the user, it chooses a unique representative for each such
  cycle by finding the "earliest" node, based on its index, to start
  from.

  As documented in the \code{igraph} package, the \code{all_simple_paths}
  function scales exponentially and may run out of memory, If that
  happens (or if you suspect that will happen), you can fall back to the
  \code{discoverCycles} function. This function makes no claim to
  completeness. Instead, it simply makes a lot of random walks through
  the graph, testing each time whether the walk revealed a cycle. If a
  consecutive set of \code{futililty} random walks fails to discover a
  new cycle, it halts.
}
\value{
  The \code{findCycles} function returns a list. Each element of the
  list is a cycle, repreented by a named character vector specifying the
  nodes in the order that they can be traversed. The
  \code{discoverCycles} function forces its output into the same format.

  The \code{interpretCycle} function returns a matrix with two columns,
  \code{genes} and \code{arrows}. Each row contains the gene name (or
  label) of a node and the form of the interaction arrow relating it
  to the next node in the cycle. Likely to be most often used inside an
  \code{lapply} function in order to interpret all cycles at once.

  The \code{cycleSubgraph} function returns an \code{igraph}
  object. This value repesents the subgraph of the full graph defined by
  all nodes that are part of at least one cycle.
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Polina Bombina \email{pbombina@augusta.edu}
}
\examples{
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
graf <- GPMLtoIgraph(xmlfile)
cyc <- findCycles(graf)
cyc
CS <- cycleSubgraph(graf, cyc)
plot(CS)
}
\keyword{ graphs }
