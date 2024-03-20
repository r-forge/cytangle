# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

## Not sure if we should export this or just keep it internal,
## only used in an overall wrapper.
collectNodes <- function(xmldoc) {
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

  ## We use the terms "node" and "vertex" interchangeably. In the
  ## GPML spec, these are called "DataNodes".
  nodes <- getNodeSet(xmlRoot(mydoc), "/sm:Pathway/sm:DataNode", rasp)
  ## Allocate space to hold the result.
  nodeInfo <- matrix(NA, nrow = length(nodes), ncol = 3)
  colnames(nodeInfo) <- c("GraphId", "label", "Type")
  ## Iterate through the GPML nodes. All the information we care about
  ## at this point is stores as XML attributes of the DataNode object.
  rowcount <- 0
  R <- rep(NA, length(nodes))
  for (node in nodes) {
    rowcount <- rowcount + 1
    nid <- xmlGetAttr(node, "GraphId")
    if (is.null(nid)) {
      warning("Nodes: Node ", rowcount,  " has no GraphId! Creating our own\n")
      nid <- paste("Node", rowcount, sep = "")
    }
    label <- xmlGetAttr(node, "TextLabel")
    label <- gsub("[\r\n]", "", label)
    type <- xmlGetAttr(node, "Type")
    if (is.null(type)) type <- "Undefined"
    repl <- c(nid, label, type)
    if (length(repl) != 3) {
      stop("Nodes: Bad replacement: ", paste(repl, collapse =", "))
    }
    nodeInfo[rowcount, ] <- c(nid, label, type)
    R[rowcount] <- nid
  }
  rownames(nodeInfo) <- R
  nodeInfo
}
