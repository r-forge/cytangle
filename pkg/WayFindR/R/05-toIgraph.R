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
    stop("Pathway contains NO edges!")
    return(NA)
  }
  groups <- collectGroups(mydoc)
  if(nrow(groups$nodes) > 0) nodes <- rbind(nodes, groups$nodes)
  if(nrow(groups$edges) > 0) edges <- rbind(edges, groups$edges)
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
  links <- collectAnchors(mydoc)
  if(nrow(links$nodes) > 0) nodes <- rbind(nodes, links$nodes)
  if(nrow(links$edges) > 0) edges <- rbind(edges, links$edges)

  if (debug) {
    cat("E =", class(edges), "N = ", class(nodes), "\n")
  }
  
  any(duplicated(nodes[,"GraphId"]))
  all(edges[,"Source"] %in% nodes[,"GraphId"])
  all(edges[,"Target"] %in% nodes[,"GraphId"])

  ## Set up colors and linestyles here. Also simplify edge types.
  simpleEdges <- simplifyArrows(edges[,"MIM"])
  if(any(is.na(simpleEdges))) stop("Bad edge type!")
  edges[,"MIM"] <- simpleEdges
  edges[,"color"] <- edgeColors[simpleEdges]
  edges[,"lty"]   <- edgeTypes[simpleEdges]
  nodes[,"color"] <- nodeColors[as.character(nodes[,"Type"])]
  nodes[,"shape"] <- nodeShapes[as.character(nodes[,"Type"])]

  mygraph <-   graph_from_data_frame(edges,
                                     directed = TRUE,
                                     vertices = nodes)
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
  legend(x, y, legend = daft$Type, col = daft$color,
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
