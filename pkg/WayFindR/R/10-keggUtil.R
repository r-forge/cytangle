# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

cpdcache <- new.env()
getIUPAC <- function(cnum) {
  if (exists(cnum, where = cpdcache)) {
    label <- cpdcache[[cnum]]
  } else {
    ans <- get_cids(cnum)
    cid <- as.data.frame(ans)[1,2]
    if (cid == "No CID") {
      label <- cnum
    } else {
      ans <- get_properties("IUPACName", cid)
      label <- ans[[1]]$IUPACName
    }
    assign(cnum, label, envir = cpdcache)
  }
  label
}

collectEntries <- function(xmldoc) {
  if (inherits(xmldoc, "XMLInternalDocument")) {
    mydoc <- xmldoc
    xmldoc <- "internal"
  } else {
    if (!file.exists(xmldoc)) {
      stop("Cannot locate file '", xmldoc, "'!")
    }
    mydoc <- xmlParseDoc(xmldoc)         # read/load the file
  }

  ## KGML uses "Entry" for what we want to call a "node" or a
  ## "vertex" in the fina grph.
  entries <- getNodeSet(xmlRoot(mydoc), "/pathway/entry")
  ## Allocate space to hold the result
  nodeInfo <- matrix(NA, nrow = length(entries), ncol = 3)
  colnames(nodeInfo) <- c("GraphId", "label", "Type")
  ## avoid flopping back and forth between matrices and data frames
  nodeInfo <- as.data.frame(nodeInfo)
  ## Iterate through the KGML entries.
  ## Need to handle different node types differently to get the info we want.
  rowcount <- 0;
  R <- rep(NA, length(entries)) # in case we couldn't process something.
  gmark <- 0
  for (entry in entries) {
    typ <- xmlGetAttr(entry, "type")
    rowcount <- rowcount + 1
    nid <- xmlGetAttr(entry, "id")
    nam <- strsplit(xmlGetAttr(entry, "name"), " ")[[1]]
    if (typ == "gene") { # figure out gene labels
      key <- gsub("hsa:", "", nam )
      sel <- suppressMessages(select(org.Hs.eg.db, keys = key,
                                     columns = c("SYMBOL", "GENETYPE")))
      sym <- paste(sel$SYMBOL, collapse = ",")
      subtyp <- paste(unique(sel$GENETYPE), collapse = ",")
      repl <- c(nid, sym, paste(typ, subtyp, sep = "|"))
    } else if (typ == "compound") {
      ctype <- strsplit(nam, ":")
      key <- ctype[[1]][2] # prefix could be 'cpd' or 'gl' or ...
      label <- getIUPAC(key)
      repl <- c(nid, label, "compound")
      Sys.sleep(1)
    } else if (typ %in% c("map", "ortholog")) {
      gent <- getNodeSet(entry, "./graphics")[[1]]
      label <- xmlGetAttr(gent, "name")
      repl <- c(nid, label, typ)
    } else if (typ == "group") {
      label <- xmlGetAttr(entry, "name")
      if (label == "undefined") { # why??
        gmark <- gmark + 1
        label <- paste("Group", gmark, sep = "")
        repl <- c(nid, label, "group")
      }
    } else {
      stop("Bad entry type", typ, "\n")
    }
    if (length(repl) != 3) {
      stop("Nodes: Bad replacement: ", paste(repl, collapse =", "))
    }
    nodeInfo[rowcount, ] <- repl
    R[rowcount] <- nid
  }
  rownames(nodeInfo) <- R
  nodeInfo
}

## There are two very different kinds of edges in KEGG:
##   relatoins (usually between genes), and
##   reactions (between compounds)

collectReactions <- function(xmldoc) {
  if (inherits(xmldoc, "XMLInternalDocument")) {
    mydoc <- xmldoc
    xmldoc <- "internal"
  } else {
    if (!file.exists(xmldoc)) {
      stop("Cannot locate file '", xmldoc, "'!")
    }
    mydoc <- xmlParseDoc(xmldoc)         # read/load the file
  }
  reactions <- getNodeSet(xmlRoot(mydoc), "/pathway/reaction")
  ## Allocate space to store the result.
  reacInfo <- matrix(NA, nrow = length(reactions), ncol = 3)
  colnames(reacInfo) <- c("Source", "Target", "MIM")
  reacInfo <- as.data.frame(reacInfo)
  ## Iterate through the KGML reactions.
  rowcount <- 0;
  for (reaction in reactions) {
    rowcount <- rowcount + 1
    typ <- xmlGetAttr(reaction, "type")
    substrate <- getNodeSet(reaction, "./substrate")[[1]]
    src <- xmlGetAttr(substrate, "id")
    product <- getNodeSet(reaction, "./product")[[1]]
    tgt <- xmlGetAttr(product, "id")
    repl <- c(src, tgt, typ)
    if (length(repl) != 3) {
      stop("Reactions: Bad replacement: ", paste(repl, collapse =", "))
    }
    reacInfo[rowcount, ] <- repl
  }
  reacInfo
}


collectRelations <- function(xmldoc) {
  if (inherits(xmldoc, "XMLInternalDocument")) {
    mydoc <- xmldoc
    xmldoc <- "internal"
  } else {
    if (!file.exists(xmldoc)) {
      stop("Cannot locate file '", xmldoc, "'!")
    }
    mydoc <- xmlParseDoc(xmldoc)         # read/load the file
  }
  relations <- getNodeSet(xmlRoot(mydoc), "/pathway/relation")
  edgeInfo <- matrix(NA, nrow = length(relations), ncol = 3)
  colnames(edgeInfo) <- c("Source", "Target", "MIM")
  edgeInfo <- as.data.frame(edgeInfo)
  ## Iterate through the KGML relations.
  rowcount <- 0;
  for (relation in relations) {
    typ <- xmlGetAttr(relation, "type")
    rowcount <- rowcount + 1
    src <- xmlGetAttr(relation, "entry1")
    tgt <- xmlGetAttr(relation, "entry2")
    sub <- getNodeSet(relation, "./subtype")[[1]]
    subtyp <- xmlGetAttr(sub, "name")
    extra <-xmlGetAttr(sub, "value")
    repl <- c(src, tgt, paste(c(typ, subtyp, extra), collapse = "/"))
    if (length(repl) != 3) {
      stop("Edges: Bad replacement: ", paste(repl, collapse =", "))
    }
    edgeInfo[rowcount, ] <- repl
  }
  edgeInfo
}
