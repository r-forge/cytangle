## disentangle.R
## Copyright (C) 2022 Kevin R. Coombes, RB McGee, and Jake Reed
## LICENSE: Perl Artistic License 2.0

indexCycles <- function(cycle, dataset) {
  L <- length(dim(cycle))
  if (L == 2) return(cycle)
  ## Now we have to convert the data-set coordinates to indicies
  N1 <- ncol(dataset)
  N2 <- dim(cycle)[3]
  if (N1 != N2) {
    stop("Mismatched sizes. Dataset = ", N1, " and cycle = ", N2, "\n")
  }
  input <- apply(dataset, 1, paste, collapse = "|")
  dex <- apply(cycle, c(1, 2), function(X) {
    check <- paste(X, collapse = "|")
    w <- which(input == check)
    if (length(w) > 1)
      w <- min(w)
    w
  })
#  indexes <- unlist(dex)
#  len <- length(indexes)
#  cols <- dim(dex)[2]
#  ifelse(cols < 2, 
#         dex <- matrix(indexes, ncol = cols, 
#                       nrow = len),
#         dex <- matrix(indexes, ncol = N2, 
#                       nrow = len/2))
  return(dex)
}

lookup <- function(loc, dataset) {
  input <- apply(dataset, 1, paste, collapse = "|")
  apply(loc, 1, function(L) {
    1.0*which(input == paste(L, collapse = "|"))
  })
}


disentangle <- function(rips, dataset) {
  cl <- rips$cycleLocation
  fixed <- lapply(cl, indexCycles, dataset = dataset)
  rips$cycleLocation <- fixed
  rips$birthLocation <- lookup(rips$birthLocation, dataset)
  rips$deathLocation <- lookup(rips$deathLocation, dataset)
  return(rips)
}
