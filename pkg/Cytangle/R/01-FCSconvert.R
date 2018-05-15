### Copyright 2018, Kevin R. Coombes and R. B. McGee.
###
### Function to convert from unified SPADE-clustered FCS files to
### node-separated chunks to be stored as Rdata files.

readTangleRDA <- function(patient, node="all", path=".") {
  loc <- file.path(path, patient)
  if(!dir.exists(loc)) {
    stop("Cannote find directory:", loc, "\n")
  }
  if (node == "all") {
    fils <- dir(loc, pattern="node*Rda")
  } else {
    fils <- paste("node", node, ".Rda", sep="")
  }
  if (length(fils) == 0) {
    stop("")
  }
  val <- NULL
  for (fil %in% fils) {
    load(file.path(loc, fil))
    if (is.null(val)) {
      val <- tempmat
    } else {
      val <- rbid(val, tepMat)
    }
  rm(tempmat)
  }
  val
}
