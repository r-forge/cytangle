### Copyright 2018, Kevin R. Coombes and R. B. McGee.
###
### Function to convert from unified SPADE-clustered FCS files to
### node-separated chunks to be stored as Rdata files.

readTangleRDA <- function(patient, node="all", path=NULL) {
  if (is.null(path)) {
    loc <- patient
  } else {
    loc <- file.path(path, patient)
  }
  if(!dir.exists(loc)) {
    stop("Cannot find directory: ", loc, "\n")
  }
  if (length(node) == 1 && node == "all") {
    fils <- dir(loc, pattern="node(.*)Rda")
  } else {
    fils <- paste("node", node, ".Rda", sep="")
  }
  if (length(fils) == 0) {
    stop("No files to process.\n")
  }
  val <- NULL
  for (fil in fils) {
    load(file.path(loc, fil))
    if (is.null(val)) {
      val <- tempMat
    } else {
      val <- rbind(val, tempMat)
    }
    rm(tempMat)
  }
  val
}
