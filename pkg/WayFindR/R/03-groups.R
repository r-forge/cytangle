# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

## Not sure if we should export this or just keep it internal,
## only used in an overall wrapper.
##
## We have already written functions to collect all edges (Interactions)
## and all vertices (DataNodes) from a GPML file. The third major type
## of object in a GPML file is called a "Group". The semantics of a
## Group are not completely clear. In some cases, they represent a complex
## and there is also a separate DataNode with attribute Type="Complex" to
## go along with it. In other cases, they sometimes appear to represent
## genes with interchangeable roles in the pathway. For example, the
## IGF1-AKT pathway includes an anonymous group with SMAD2 and SMAD3.
collectGroups <- function(xmldoc, allnodes) {
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
  nodal <- data.frame(GraphId = NA, label = NA, Type = NA)[-1,]
  edger <- data.frame(Source = NA, Targte = NA, MIM = NA)[-1,]

  ## Now we get the actual group objects
  ghash <- NULL # map GroupId to GraphId
  groups <- getNodeSet(xmlRoot(mydoc), "/sm:Pathway/sm:Group", rasp)
  gcounter <- 0
  for (group in groups) {
    gid <- xmlGetAttr(group, "GraphId")
    grid <- xmlGetAttr(group, "GroupId")
    if (is.null(gid)) {
      warning("Groups: GroupId has no GraphId. Using GroupId.\n")
      gid <- grid
    } else if (gid != grid) {
      warning("Groups: GroupId does not match GraphId.\n")
    }
    ghash[grid] <- gid
    ## See if there is a separate node that contains a name for the complex.
    ## Key point: XML package doesn't know how to expand variable names,
    ## so you need to do so explicitly, by hand, to build a query.
    query <- paste("/sm:Pathway/sm:DataNode[@GroupRef='", gid,
                   "' and @Type='Complex']", sep = "")
    repr <- getNodeSet(xmlRoot(mydoc), query, rasp)
    if (length(repr) > 0) {
      ## update the node that represent the complex
      lbl <- xmlGetAttr(repr[[1]], "TextLabel")
    } else {
      gcounter <- gcounter + 1
      lbl <- paste("Group", gcounter, sep = "")
    }
    sty <- xmlGetAttr(group, "Style")
    if (is.null(sty)) sty <- "Group"
    newn <- data.frame(GraphId = gid,
                       label = lbl,
                       Type = sty)
    rownames(newn) <- newn[,1]
    newn <- as.matrix(newn)
    nodal <- rbind(nodal, newn)
  }

  ## We start by finding all DataNodes that include a GroupRef attribute.
  ## Most of these consist of gene products or proteins that are part of
  ## the group (indicated by the Type attribute). In some cases, however,
  ## there is one DataNode that represents the group as a whole, having
  ## Type = "Complex".
  grefs <- getNodeSet(xmlRoot(mydoc), "/sm:Pathway/sm:DataNode[@GroupRef]", rasp)
  edgeCounter <- 0
  for (gref in grefs) {
    grf <- xmlGetAttr(gref, "GroupRef")
    if (!(grf %in% names(ghash))) {
      stop("Groups: Reference to non existent group!\n")
    }
    nam <- xmlGetAttr(gref, "TextLabel")
    gid <- xmlGetAttr(gref, "GraphId")
    if (is.null(gid)) {
      warning("Groups: Node ", nam, " has no GraphId! Need to locate one.\n")
      who <- which(allnodes[, "label"] == nam)
      if (length(who) == 0) { # should be impossible
        stop("Groups: Cannot locate node without a graphId!\n")
      } else if (length(who) > 1) { # might happen, but it's still bad
        stop("Groups: Two different nodes have the same label!\n")
      }
      gid <- allnodes[who, "GraphId"]
    }
    typ <- xmlGetAttr(gref, "Type")
    if (is.null(typ)) {
      typ <- "Complex" #??
    }
    if (typ == "Complex") { # Already handed this case
      next
    } else if (typ %in% c("GeneProduct", "Metabolite",
                          "Pathway", "Protein", "RNA", "Rna",
                          "Other", "CellularComponent")) {
      ## create a "contained" edge
      edgeCounter <- edgeCounter + 1
      nedg <- data.frame(Source = gid, Target = ghash[grf], MIM = "contained")
      rownames(nedg) <- paste("ec", edgeCounter, sep = "")
      edger <- rbind(edger, nedg)
    } else {
      stop("Groups: Unexpected node type (", typ, ") in gref = ", grf,
           " at gid = ", gid, "\n")
    }
  }

  list(nodes = nodal, edges = edger)
}
