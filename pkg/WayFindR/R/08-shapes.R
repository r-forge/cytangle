# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

## Not sure if we should export this or just keep it internal,
## only used in an overall wrapper.
collectShapes <- function(xmldoc) {
  if (inherits(xmldoc, "XMLInternalDocument")) {
    mydoc <- xmldoc
    xmldoc <- "internal"
  } else {
    if (!file.exists(xmldoc)) {
      stop("Cannot locate file '", xmldoc, "'!")
    }
    mydoc <- xmlParseDoc(xmldoc)         # read/load the file
  }
  nsp <- xmlNamespace(xmlRoot(mydoc)) # extract the namespace
  rasp <- c(sm = as.character(nsp))   # assign abbreviation to the namespace

  shapes <- getNodeSet(xmlRoot(mydoc), "/sm:Pathway/sm:Shape", rasp)
  nick <- matrix(NA, nrow = length(shapes), ncol = 3)
  colnames(nick) <- c("GraphId", "label", "Type")
  rowcount <- 0
  R <- rep(NA, length(shapes))
  for (shape in shapes) {
    nid <- xmlGetAttr(shape, "GraphId")
    if (is.null(nid)) next
    type <- "Shape"
    rowcount <- rowcount + 1
    shape <- paste("Shape", rowcount, sep = "")
    nick[rowcount, ] <- c(nid, shape, type)
    R[rowcount] <- nid
  }
  nick <- nick[!is.na(nick[,1]),, drop = FALSE]
  nick
}
