# (C) Copyright 2024 Kevin R. Coombes and Polina Bombina

## Original annotation method are slow, using separate HTTP calls
## to a server for each entry.
WAYcache <- new.env()

prepAnno <- function(items) {
  tags <- sapply(items, function(A) {
    xmlGetAttr(A, "name")  })
  pop <- strsplit(tags, ":")
  pref <- sapply(pop, function(X) X[1])
  id <-  sapply(pop, function(X) X[2])
  cpd <- id[pref == "cpd"]
  cache <- function(G) {
    ent <- G$ENTRY
    lbl <- G$NAME[1]
    assign(ent, lbl, WAYcache)
    invisible(ent)
  }
  while(length(cpd) > 10) {
    kg <- keggGet(cpd[1:10])
    ignore <- sapply(kg, cache)
    cpd <- cpd[11:length(cpd)]
  }
  if (length(cpd) > 0) {
    kg <- keggGet(cpd)
    ignore <- sapply(kg, cache)
  }
  gly <- id[pref == "gl"]
  while(length(gly) > 10) {
    kg <- keggGet(gly[1:10])
    ignore <- sapply(kg, cache)
    gly <- gly[11:length(gly)]
  }
  if (length(gly) > 0) {
    kg <- keggGet(gly)
    ignore <- sapply(kg, cache)
  }
}

getIUPAC <- function(cnum) {
  if (exists(cnum, where = WAYcache)) {
    label <- WAYcache[[cnum]]
  } else {
    suppressMessages( ans <- CIDs(get_cids(cnum)) )
    if (ncol(ans) < 2) {
      label <- cnum
    } else {
      cid <- as.data.frame(ans)[1,2]
      ans <- get_properties("IUPACName", cid)
      R <- as.data.frame(retrieve(ans))
      label <- R$IUPACName
      }
    assign(cnum, label, envir = WAYcache)
  }
  label
}

getIUPACAll <- function(cnum) {
  if (!exists("kgc", where = WAYcache)) {
    kgc <- keggList("glycan")
    assign("kgc", kgc, envir = WAYcache)
  }
  val <- WAYcache$kgc[cnum]
  if (is.na(val) || val == "") {
    val <- cnum
  } else {
    val <- strsplit(val, "; ")[[1]][1]
  }
  val
}


getGlycan <- function(gnum) {
  kg <- keggGet(gnum)
  if (length(kg) < 1) {
    val <- gnum
  } else {
    gunk <- kg[[1]]
    val <- gunk$COMPOSITION
  }
  val
}

getGlycanAll <- function(gnum) {
  if (!exists("kgl", where = WAYcache)) {
    kgl <- keggList("glycan")
    assign("kgl", kgl, envir = WAYcache)
  }
  val <- WAYcache$kgl[gnum]
  if (is.na(val) || val == "") {
    val <- gnum
  }
  val
}

collectEntries <- function(xmldoc, anno = c("all", "one", "batch")) {
  if (inherits(xmldoc, "XMLInternalDocument")) {
    mydoc <- xmldoc
    xmldoc <- "internal"
  } else {
    if (!file.exists(xmldoc)) {
      stop("Cannot locate file '", xmldoc, "'!")
    }
    mydoc <- xmlParseDoc(xmldoc)         # read/load the file
  }
  anno <- match.arg(anno)
  glyanno <- switch(anno,
                  all = getGlycanAll,
                  batch = getGlycan, # pending
                  one = getGlycan)
  cpdanno <- switch(anno,
                  all = getIUPACAll,
                  batch = getIUPAC, # pending
                  one = getIUPAC)

  ## KGML uses "Entry" for what we want to call a "node" or a
  ## "vertex" in the fina grph.
  entries <- getNodeSet(xmlRoot(mydoc), "/pathway/entry")
  if (anno == "batch") {
    prepAnno(entries)
  }
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
      sy <- sel$SYMBOL
      if (length(sy) > 3) sy <- sy[1:3]
      sym <- paste(sy, collapse = ",")
      repl <- c(nid, sym, typ)
      self <- nam[1]
    } else if (typ == "compound") {
      ctype <- strsplit(nam, ":")
      key <- ctype[[1]][2] # prefix could be 'cpd' or 'gl' or ...
      tag <- ctype[[1]][1]
      if (tag == "gl") {
        label <- glyanno(key)
      } else { ##if (tag == "cpd")
        label <- cpdanno(key)
      }
      repl <- c(nid, label, "compound")
      self <- key
    } else if (typ %in% c("map", "ortholog")) {
      key <- getNodeSet(entry, "./graphics")[[1]]
      label <- sub("^TITLE:", "", xmlGetAttr(key, "name"))
      repl <- c(nid, label, typ)
      self <- sub("ko:", "", nam[1])
    } else if (typ == "group") {
      label <- xmlGetAttr(entry, "name")
      if (label == "undefined") { # why??
        gmark <- gmark + 1
        label <- paste("Group", gmark, sep = "")
        repl <- c(nid, label, "group")
      }
      self <- label
    } else {
      stop("Bad entry type", typ, "\n")
    }
    if (length(repl) != 3) {
      stop("Entries: Bad replacement: ", paste(repl, collapse =", "))
    }
    nodeInfo[rowcount, ] <- repl
    if (length(self) > 1) stop("Entries: Bad name")
    R[rowcount] <- self
  }
  while(any(dd <- duplicated(R))) {
    R[dd] <- paste(R[dd], ".", sep = "")
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
    sub <- getNodeSet(relation, "./subtype")
    if (length(sub) == 0) {
      val <- typ
    } else {
      subtyp <- paste(sapply(sub, function(S) xmlGetAttr(S, "name")),
                      collapse = ",")
      val <- paste(c(typ, subtyp), collapse = "|")
    }
    repl <- c(src, tgt, val)
    if (length(repl) != 3) {
      stop("Edges: Bad replacement: ", paste(repl, collapse =", "))
    }
    edgeInfo[rowcount, ] <- repl
  }
  edgeInfo
}
