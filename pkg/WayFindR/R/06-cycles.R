# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

## character vector representation
shift <- function(spl, N) {
  L <- length(spl)
  if (N > L)
    warning("Shifting cycles by more than their length is undefined.")
  spl[c((N+1):L, 1:N)]
}

## hyphern-separated string representation
slide <- function(seq, N) {
  moved <- shift(strsplit(seq, "-")[[1]], N)
  paste(moved, collapse = "-")
}

## Don't try to find all cycles. Settle for whatever you
## can discover through random walks.
discoverCycles <- function(graph, nsteps = 50, futility = 20) {
  Cycles <-  NULL    # eventually, a list/hash for cycles
  futile <- futility # number of allowed steps without seeing a new cycle
  while(futile > 0) {
    startNode <- sample(length(V(graph)), 1)
    rw <- as_ids(random_walk(graph, 1, nsteps))
    found <- 0
    N = 0
    while(any(dup <- duplicated(rw))) {
      ## Lots of switching back-and-forth between hyphen-separeted
      ## strings (cyc) and character vectors (bits)
      uu <- unique(rw[dup])
      pending <- which(rw %in% uu)[1]
      A <- rw[pending]
      w <- which(rw == A)
      bits <- rw[w[1]:(w[2]-1)]
      cyc <- paste(bits, collapse = "-")
      if (any(bdup <- duplicated(bits))) { # figure-eights?
        ell <- which(bits %in% bits[which(bdup)[1]])
        bits <- shift(bits, ell[1] - 1)
        ww <- which(duplicated(bits))[1] - 1
        bits <- bits[1:ww]
        cyc <- paste(bits, collapse = "-")
      } 
      oop <- order(bits)
      w0 <- which(oop == 1)
      if (w0 > 1) {
        cyc <- slide(paste(bits, collapse = "-"), w0-1)
      }
      if (cyc %in% names(Cycles)) {
        Cycles[[cyc]] <- Cycles[[cyc]] + 1
      } else {
        found <- 1
        Cycles[[cyc]] <- 1
        N <- length(Cycles)
        if (N %% 1000 == 0) {
          cat("N = ", N,  "\n", file = stderr())
        }
        if (N > 50000) break
      }
      rw <- rw[(w[2]+1):length(rw)]
    } # end (while(dup))
    if (found == 0) futile <- futile - 1
    if (N > 50000) break
  } # end while (futile)
  if (!is.null(Cycles)) {
    Cycles <- unlist(Cycles)
    temp <- strsplit(names(Cycles), "-")
    gack <- as_ids(V(graph))
    idx <- 1:length(gack)
    names(idx) <- gack
    dumb <- sapply(temp, function(X) idx[X])
    fc <- sapply(dumb, function(X) c(X[length(X)], X))
  } else {
    fc <- NULL
  }
}

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

