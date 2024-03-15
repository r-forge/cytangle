# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

## Not sure if we should export this or just keep it internal,
## only used in an overall wrapper.
collectEdges <- function(xmldoc) {
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

  ## Edges are called Interactions in the GPML spec.
  edges <- getNodeSet(xmlRoot(mydoc), "/sm:Pathway/sm:Interaction", rasp)
  ## Allocate space to store the results.
  edgeList <- matrix(NA, nrow = length(edges), ncol = 3)
  colnames(edgeList) <- c("Source", "Target", "MIM")
  ## Iterate through the GPML edges.
  rowcount <- 0;
  R <- rep(NA, length(edges))
  for (edge in edges) {
    anchors <- getNodeSet(edge, "./sm:Graphics/sm:Anchor", rasp)
    if (length(anchors) > 0) { # This edge points to another edge
      ## Skip this for now and deal with it after the basics are put together
      next
    }
    ## Each interaction should contain one Graphics object, with two Points
    ## (the source and target of the edge). The second point should contain
    ## an ArrowHead that designates the semantic meaning of the edge.
    rowcount = rowcount + 1
    eid <- xmlGetAttr(edge, "GraphId")
    pts <- getNodeSet(edge, "./sm:Graphics/sm:Point", rasp)
    if (length(pts) != 2) { # should be impossible
      stop("Got wrong number of points in an interaction in ", xmldoc, "!\n")
    }
    counter <- 0
    for (point in pts) {
      counter <- counter + 1
      node = xmlGetAttr(point, "GraphRef")
      arrow = xmlGetAttr(point, "ArrowHead")
      if (counter == 1) {
        if(!is.null(arrow)) { # may be impossible?
          stop("Source node has an arrow head!\n")
        } else {
          src <- node
        }
      } else {
        arrow <- ifelse(is.null(arrow), "none", arrow)
        tgt <- node
        edgeList[rowcount, ] <- c(src, tgt, arrow)
        R[rowcount] <- eid
      }
    }
  }
  rownames(edgeList) <- R
  ## Only have empty rows if there were anchors, indicating edges that
  ## point to other edges.
  empty <- apply(edgeList, 1, function(X) any(is.na(X)))
  edgeList <- edgeList[!empty,]
  edgeList
}