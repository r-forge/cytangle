# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

findCycles <- function(graph) {
  ## Online reference:
  ## https://stackoverflow.com/questions/55091438/r-igraph-find-all-cycles

  ## Loop over every node, and every starting edge
  Cycles = NULL
  for(v1 in V(graph)) {
    for(v2 in neighbors(graph, v1, mode="out")) {
      Cycles = c(Cycles, 
                 lapply(all_simple_paths(graph, v2, v1, mode="out"),
                        function(p) c(v1,p)))
    }
  }
  ## current list is redundant, since a loop of length K appears K times.
  ## Want just unique versions of ewach cycle.
  UniqueCycles <- Cycles[sapply(Cycles, min) == sapply(Cycles, `[`, 1)]
  UniqueCycles
}

## V must be a cycle in the graph
interpretCycle <- function(v, graph) {
  eid <- get.edge.ids(graph, c(rbind(head(v, -1), v[-1])))
  arrows <- E(graph)$MIM[eid]
  genes <- V(graph)$label[v[-length(v)]]
  data.frame(genes, arrows)
}

cycleSubgraph <- function(graph, cycles) {
  uu <- unique(unlist(sapply(cycles, names)))
  uu <- uu[uu != ""]
  subgraph(graph, uu)
}

