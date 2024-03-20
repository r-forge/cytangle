# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

## Not sure if we should export this or just keep it internal,
## only used in an overall wrapper.
collectLabels <- function(xmldoc) {
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

  labels <- getNodeSet(xmlRoot(mydoc), "/sm:Pathway/sm:Label", rasp)
  nick <- matrix(NA, nrow = length(labels), ncol = 3)
  colnames(nick) <- c("GraphId", "label", "Type")
  rowcount <- 0
  R <- rep(NA, length(labels))
  for (lbl in labels) {
    nid <- xmlGetAttr(lbl, "GraphId")
    if (is.null(nid)) next
    label <- xmlGetAttr(lbl, "TextLabel")
    label <- gsub("[\r\n]", "", label)
    type <- "Label"
    rowcount <- rowcount + 1
    nick[rowcount, ] <- c(nid, label, type)
    R[rowcount] <- nid
  }
  nick <- nick[!is.na(nick[,1]),, drop = FALSE]
  nick
}
