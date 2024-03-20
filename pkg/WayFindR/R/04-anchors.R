# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

## Not sure if we should export this or just keep it internal,
## only used in an overall wrapper.
##
## The fourth and final structure of interest in the GPML specification
## is the Anchor. An Anchor is used to mark edges (or interactions, or
## arrows) that are the target of another edge instead of the usual node
## (or DataNode, or vertex). We deal with this by changing the target
## arrow, A -> B, into a two-step arrow, A -> EDGE -> B, and make the
## funny arrow into a simpler, more standard thing tha tnow points to the
## new EDGE node.
collectAnchors <- function(xmldoc) {
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

  ## Don't know how many things we will get, so just define empty structures
  nodal <- data.frame(GraphID = NA, label = NA, Type = NA)[-1,]
  edger <- data.frame(Source = NA, Target = NA, MIM = NA)[-1,]

  ## Find all the anchors
  anchors <- getNodeSet(xmlRoot(mydoc),
                        "/sm:Pathway/sm:Interaction/sm:Graphics/sm:Anchor",
                        rasp)
  ## Iterate over the anchors
  acount <- 0
  for (anchor in anchors) {
    acount <- acount + 1
    lbl <- paste("Anchor", acount, sep = "")
    gfx <- xmlParent(anchor) # must be a Graphics object
    edge <- xmlParent(gfx)   # must be an Interaction (edge) object
    ## Create a node of type "EDGE" so that the edge that is pointed to
    ## gets transformed to A -> EDGE -> B so that we can make a node for
    ## the target of the one we are working on.
    gid <- xmlGetAttr(edge, "GraphId")
    aid <- xmlGetAttr(anchor, "GraphId")
    if (length(aid) == 0) {
      warning("Anchors: Skipping Anchor that has no GraphId!\n")
      next
    }
    edgenode <- data.frame(GraphId = aid,
                           label = lbl,
                           Type = "EDGE")
    pts <- getNodeSet(edge, "./sm:Graphics/sm:Point", rasp)
    src <- xmlGetAttr(pts[[1]], "GraphRef")
    tgt <- xmlGetAttr(pts[[2]], "GraphRef")
    ## Create two edges; src -> EDGE, EDGE -> tgt
    newedges <- data.frame(Source = c(src, aid),
                           Target = c(aid, tgt),
                           MIM    = c("Source",
                                      xmlGetAttr(pts[[2]], "ArrowHead")))
    nodal <- rbind(nodal, edgenode)
    edger <- rbind(edger, newedges)
  }
  list(nodes = nodal, edges = edger)
}
