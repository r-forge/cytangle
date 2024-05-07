# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

translateKGMLnodes <- function(V) {
  remap <- c(map = "Pathway",
             group = "Group",
             compound = "Metabolite", ## update this; perhaps to "SmallMolecule"
             gene = "GeneProduct",
             ortholog = "Ortholog")
  sapply(V, function(X) {
    ifelse(X %in% names(remap), remap[X], X)
  })
}

translateKGMLedges <- function(V) {
  W <- strsplit(V, "\\|") # obviously W is a spread-out V
  sapply(W, function(R) {
    if (length(R) == 1) {
      if (R %in% c("reversible", "irreversible")) {
        return(R)
      }
      subtyp <- "none"
    } else {
      subtyp <- strsplit(R[2], "\\.")[[1]]
    }
    typ <- R[1]
    if (any(subtyp %in% c("inhibition", "repression")))
      return("mim-inhibition")
    if (any(subtyp %in% c("activation", "expression")))
      return("mim-stimulation")
    if (any(subtyp %in% c("state change", "phosphorylation", "dephosphorylation",
                          "glycosylation", "ubiquitination")))
      return("mim-modification")
    if (any(subtyp %in% c("binding/association", "dissociation")))
      return("mim-binding")
    if (typ == "maplink")
      return("pathway")
    if (typ == "ECrel" | subtyp == "none")
      return("compound")
    return("none")
  })
}

KGMLtoIgraph <- function(xmldoc, returnLists = FALSE, debug = FALSE) {
  if (inherits(xmldoc, "XMLInternalDocument")) {
    mydoc <- xmldoc
    xmldoc <- "internal"
  } else {
    if (!file.exists(xmldoc)) {
      stop("Cannot locate file '", xmldoc, "'!")
    }
    mydoc <- xmlParseDoc(xmldoc)         # read/load the file
  }

  nodes <- collectEntries(mydoc)
  retype <- translateKGMLnodes(nodes$Type)
  nodes$Type <- retype
  relns <- collectRelations(mydoc)
  reacs <- collectReactions(mydoc)
  edges <- rbind(relns, reacs)
  if (nrow(edges) == 0) {
    stop("WayFindR: KGML pathway contains NO edges!")
    return(NA)
  }

  if (debug) {
    cat("E =", class(edges), "N = ", class(nodes), "\n")
  }

  ## Set up colors and linestyles here. Also simplify edge types.
  simpleEdges <- translateKGMLedges(edges$MIM)
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
  if (any(!needed)) {
    nodes <- nodes[needed, , drop = FALSE]
    warning("KGML: Some nodes have no edge connections!\n")
  }

  mygraph <-   graph_from_data_frame(edges,
                                     directed = TRUE,
                                     vertices = nodes)

  name <- xmlGetAttr(xmlRoot(mydoc), "title")
  pid <- xmlGetAttr(xmlRoot(mydoc), "name")
  mygraph <- set_graph_attr(mygraph, "Name", name)
  mygraph <- set_graph_attr(mygraph, "PathID", pid)
  mygraph <- set_graph_attr(mygraph, "Authors", "KEGG")

  if (returnLists) {
    val <- list(graph = mygraph,
                edges = as.data.frame(edges),
                nodes = as.data.frame(nodes))
  } else {
    val <- mygraph
  }
  val
}
