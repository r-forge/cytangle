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
collectGroups <- function(xmldoc) {
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
  edger <- data.frame(Source = NA, Targte = NA, MIM = NA)[-1,]

  ## We start by finding all DataNodes that include a GroupRef attribute.
  ## Most of these consist of gene products or proteins that are part of
  ## the group (indicated by the Type attribute). In some cases, however,
  ## there is one DataNode that represents the group as a whole, having
  ## Type = "Complex".
  grefs <- getNodeSet(xmlRoot(mydoc), "/sm:Pathway/sm:DataNode[@GroupRef]", rasp)
  edgeCounter <- 0
  for (gref in grefs) {
    grf <- xmlGetAttr(gref, "GroupRef")
    gid <- xmlGetAttr(gref, "GraphId")
    nam <- xmlGetAttr(gref, "TextLabel")
    typ <- xmlGetAttr(gref, "Type")
    if (typ == "Complex") {
#      if (gid %in% rownames(nodeInfo)) next
      ## create an edge to represent the complex
      edgeCounter <- edgeCounter + 1
      nedg <- data.frame(Source = gid, Target = nam, MIM = "represents")
      cat("Complex!\n", file = stderr())
      rownames(nedg) <- paste("ec", edgeCounter, sep = "")
      edger <- rbind(edger, nedg)
    } else if (typ == "GeneProduct") {
      ## create a "contained" edge
      edgeCounter <- edgeCounter + 1
      nedg <- data.frame(Source = gid, Target = grf, MIM = "contained")
      rownames(nedg) <- paste("ec", edgeCounter, sep = "")
      edger <- rbind(edger, nedg)
    } else {
      stop(" Don't panic; just grab a towel.\n")
    }
  }

  ## Now we get the actual group objects
  groups <- getNodeSet(xmlRoot(mydoc), "/sm:Pathway/sm:Group", rasp)
  gcounter <- 0
  for (group in groups) {
    gcounter <- gcounter + 1
    cat(gcounter, "\n", file = stderr())
    gid <- xmlGetAttr(group, "GraphId")
    newn <- data.frame(GraphID = gid,
                       label = paste("Group", gcounter, sep = ""),
                       Type = xmlGetAttr(group, "Style"))
    rownames(newn) <- newn[,1]
    newn <- as.matrix(newn)
    nodal <- rbind(nodal, newn)
    ## Key point: XML pacakge doesn't know ho to expand variable names,
    ## so you need to do so explicitly, by hand, to build a query.
    query <- paste("/sm:Pathway/sm:DataNode[@GroupRef='", gid,
                   "' and @Type='Complex']", sep = "")
    repr <- getNodeSet(xmlRoot(mydoc), query, rasp)
    if (length(repr) > 0) {
      ## create an edge to represent the complex
      edgeCounter <- edgeCounter + 1
      nedg <- data.frame(Source = xmlGetAttr(repr[[1]], "GraphId"),
                         Target = gid,
                         MIM = "represents")
      cat("Complex!\n", file = stderr())
      rownames(nedg) <- paste("ec", edgeCounter, sep = "")
      edger <- rbind(edger, nedg)
    }
  }
  list(nodes = nodal, edges = edger)
}