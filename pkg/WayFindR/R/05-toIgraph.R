# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

simplifyArrows <- function(V) {
  remap <- c(Arrow = "mim-stimulation",
             LigandRound = NA,
             LigandSquare = NA,
             Receptor = NA,
             ReceptorRound = NA,
             ReceptorSquare = NA,
             "SBGN-Catalysis" = "mim-catalysis",
             "SBGN-Inhibition" = "mim-inhibition",
             "SBGN-Production" = "mim-stimulation",
             TBar = "mim-inhibition")
  sapply(V, function(X) {
    ifelse(X %in% names(remap), remap[X], X)
  })
}

GPMLtoIgraph <- function(xmldoc, returnLists = FALSE, debug = FALSE) {
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
  
  nodes <- collectNodes(mydoc)
  edges <- collectEdges(mydoc)
  if (nrow(edges) == 0) {
    stop("WayFindR: Pathway contains NO edges!")
    return(NA)
  }
  groups <- collectGroups(mydoc, nodes)
  if(nrow(groups$nodes) > 0) nodes <- rbind(nodes, groups$nodes)
  if(nrow(groups$edges) > 0) edges <- rbind(edges, groups$edges)
  links <- collectAnchors(mydoc)
  if(nrow(links$nodes) > 0) nodes <- rbind(nodes, links$nodes)
  if(nrow(links$edges) > 0) edges <- rbind(edges, links$edges)

  ## We may have duplicated a group/complex node. Let's check.
  cpx <- as.data.frame(nodes[nodes[, "Type"] == "Complex",])
  if (any(dup <- duplicated(cpx$label))) {
    ## figure out which duplicate to delete
    for (lbl in cpx$label[dup]) {
      X <- rownames(cpx)[cpx$label == lbl]
      for (who in X) {
        if (!(who %in% edges$Source) & !(who %in% edges$Target)) {
          ## not on any edge, so we can discard it
          nodes <- nodes[rownames(nodes) != who,]
        }
      }
    }
  }

  labels <- collectLabels(mydoc)
  if (nrow(labels) > 0) nodes <- rbind(nodes, labels)

  shapes <- collectShapes(mydoc)
  if (nrow(shapes) > 0) nodes <- rbind(nodes, shapes)

  states <- getNodeSet(xmlRoot(mydoc), "/sm:Pathway/sm:State", rasp)
  stateIds <- sapply(states, function(state) { xmlGetAttr(state, "GraphId")})
  ends <- c(edges[, "Source"], edges[, "Target"])
  if (length(states) > 0 & any(stateIds %in% ends)) {
    stop("WayFindR: Cannot yet handle pathways with 'States'!")
  }
  glines <- getNodeSet(xmlRoot(mydoc), "/sm:Pathway/sm:GraphicalLine/sm:Graphics/sm:Anchor", rasp)
  if (length(glines) > 0) {
    stop("WayFindR: Cannot yet handle pathways with 'Anchors' in 'GraphicalLines'!\n")
  }

  if (debug) {
    cat("E =", class(edges), "N = ", class(nodes), "\n")
  }
  if (any(duplicated(nodes[,"GraphId"]))) {
    stop("WayFindR: Found multiple nodes with same 'GraphId'!")
  }
  if (!all(edges[,"Source"] %in% nodes[,"GraphId"])) {
    warning("WayFindR: found an edge with unknown Source node!\n")
    # which(!(edges[,"Source"] %in% nodes[,"GraphId"]))
  }
  if (!all(edges[,"Target"] %in% nodes[,"GraphId"])) {
    warning("WayFindR: Found an edge with unknown Target node!\n")
    # which(!(edges[,"Target"] %in% nodes[,"GraphId"]))
  }

  ## Set up colors and linestyles here. Also simplify edge types.
  simpleEdges <- simplifyArrows(edges[,"MIM"])
  if(any(is.na(simpleEdges))) {
    odd <- paste(unique(names(which(is.na(simpleEdges)))),
                 collapse = ", ")
    stop("WayFindR: Bad edge type (Receptor-Ligand).\n")
  }
  edges[,"MIM"] <- simpleEdges
  edges <- as.data.frame(edges)
  nodes <- as.data.frame(nodes)
  edges[,"color"] <- edgeColors[simpleEdges]
  edges[,"lty"]   <- edgeTypes[simpleEdges]
  nodes[,"color"] <- nodeColors[as.character(nodes[,"Type"])]
  nodes[,"shape"] <- nodeShapes[as.character(nodes[,"Type"])]

  uu <- unique(c(edges$Source, edges$Target))
  needed <- nodes$GraphId %in% uu
  nodes <- nodes[needed, , drop = FALSE]

  if (any(nodes$Type == "Label")) {
    warning("Pathway uses a 'Label' as source or target of an Edge!\n")
  }
  if (any(nodes$Type == "Shape")) {
    warning("Pathway uses a 'Shape' as source or target of an Edge!\n")
  }

  mygraph <-   graph_from_data_frame(edges,
                                     directed = TRUE,
                                     vertices = nodes)

  name <- xmlGetAttr(xmlRoot(mydoc), "Name")
  version <- xmlGetAttr(xmlRoot(mydoc), "Version")
  pid <- strsplit(version, "_")[[1]][1]
  author <- xmlGetAttr(xmlRoot(mydoc), "Author")
  mygraph <- set_graph_attr(mygraph, "Name", name)
  mygraph <- set_graph_attr(mygraph, "PathID", pid)
  mygraph <- set_graph_attr(mygraph, "Authors", author)
  
  if (returnLists) {
    val <- list(graph = mygraph,
                edges = as.data.frame(edges),
                nodes = as.data.frame(nodes))
  } else {
    val <- mygraph
  }
  val
}

nodeLegend <- function(x, graph) {
  xlate <- c(circle = 16, rectangle = 15)
  daft <- data.frame(Type =V(graph)$Type,
                     color = V(graph)$color,
                     shape = V(graph)$shape)
  daft <- unique(daft)
  legend(x, legend = daft$Type, col = daft$color,
         pch = xlate[daft$shape])
}

edgeLegend <- function(x, graph) {
  daft <- data.frame(MIM = E(graph)$MIM,
                     color = E(graph)$color,
                     lty = E(graph)$lty)
  daft <- unique(daft)
  legend(x, legend = daft$MIM, col = daft$color,
         lty = daft$lty)
}
