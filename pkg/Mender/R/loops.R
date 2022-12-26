## Copright (C) 2022 Kevin R. Coombes, RB McGee, and Jake Reed

## helper function to get (longest persisting) cycle
##    rips = a Rips diagramm from TDA
##    dimension = a non-negative integer (probably less than three)
##    target = an integer cycle ID, or NULL to find longest persistence
getCycle <- function(rips, dimension = 1, target = NULL) {
  ## Find cycles of the desired dimension
  RD <- rips$diagram
  whichD <- which(RD[, "dimension"] == dimension)
  if (length(whichD) == 0) {
    stop("Diagram contains no cycles of dimension ", dimension, ".")
  }
  ## Find duration of persistence
  duration <- (RD[, 3] - RD[, 2])[whichD]
  if (is.null(target)) {
    target <- which.max(duration) # id of most persistent cycle
  }
  ## Find cycles and pull out the one requested
  cyc <- rips$cycleLocation[whichD]
  targetCycle <- cyc[target][[1]]
  targetCycle
}

## How mathematicians name coordinate axes
mathNames <- letters[c(24:26, 21:23, 16:20, 1:19)]

## Given a view/layout, extract the coordinates of the point
## with index J
getCoords <- function(J, view) {
  dimn <- ncol(view)
  coord <- view[J,]
  names(coord) <- mathNames[1:dimn]
  coord
}

## Given a view/layout and a cycle, extract the coordinates
## of all points in the support
cycleSupport <- function(cycle, view) {
  foo <- unique(as.vector(cycle))
  t(sapply(foo, getCoords, view = view))
}

cycleEdges <- function(cycle, view, ...) {
  edges <- apply(cycle, 1, function(arow) {
    bdry <- sapply(arow, getCoords, view = view)
    list(x = bdry["x",], y = bdry["y",], ...)
  })
  edges
}

showCycle <- function(cycle, view, col = "black", ...) {
  pts <- cycleSupport(cycle, view)
  edges <- cycleEdges(cyc, view, col = col, ...)
  points(pts, pch = 16, col = col)
  sapply(edges, function(E) {
    do.call(lines, E)
  })
  invisible(cycle)
}
