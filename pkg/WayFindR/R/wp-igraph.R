library(XML)

## Common Portion
fname <- "WP3850.gpml"
foo <- xmlParseDoc(fname) # read/load the file
nsp <- xmlNamespace(xmlRoot(foo)) # extract the namesspace
rasp <- c(sm = as.character(nsp)) # assign short abbreviationn to the namespace

## Edges
edges <- getNodeSet(xmlRoot(foo), "/sm:Pathway/sm:Interaction", rasp)
length(edges) # should be 47
edgeList <- matrix(NA, nrow = length(edges), ncol = 3)
colnames(edgeList) <- c("Source", "Target", "MIM")
rowcount <- 0;
R <- rep(NA, length(edges))
for (edge in edges) {
  anchors <- getNodeSet(edge, "./sm:Graphics/sm:Anchor", rasp)
  if (length(anchors) > 0) { # This edge points to another edge
    ## Skip this for now and deal with it after the basics are put together
    next
  }
  rowcount = rowcount + 1
  eid <- xmlGetAttr(edge, "GraphId")
  pts <- getNodeSet(edge, "./sm:Graphics/sm:Point", rasp)
  if (length(pts) != 2) {
    warning("Got wrong number of points in an interaction in ", fname, "!\n")
  }
  counter <- 0
  for (point in pts) {
    counter <- counter + 1
    node = xmlGetAttr(point, "GraphRef")
    arrow = xmlGetAttr(point, "ArrowHead")
    if (counter == 1) {
      if(!is.null(arrow)) {
        warning("Source node has an arrow head!\n")
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
empty <- apply(edgeList, 1, function(X) any(is.na(X)))
edgeList <- edgeList[!empty,]


## Vertices
nodes <- getNodeSet(xmlRoot(foo), "/sm:Pathway/sm:DataNode", rasp)
length(nodes) # should be 52
nodeInfo <- matrix(NA, nrow = length(nodes), ncol = 3)
colnames(nodeInfo) <- c("GraphId", "label", "Type")
rowcount <- 0
R <- rep(NA, length(nodes))
for (node in nodes) {
  rowcount <- rowcount + 1
  nid <- xmlGetAttr(node, "GraphId")
  label <- xmlGetAttr(node, "TextLabel")
  label <- gsub("[\r\n]", "", label)
  type <- xmlGetAttr(node, "Type")
  nodeInfo[rowcount, ] <- c(nid, label, type)
  R[rowcount] <- nid
}
rownames(nodeInfo) <- R

## Groups
grefs <- getNodeSet(xmlRoot(foo), "/sm:Pathway/sm:DataNode[@GroupRef]", rasp)
length(grefs) # should be 15
edgeCounter <- 0
for (gref in grefs) {
  grf <- xmlGetAttr(gref, "GroupRef")
  gid <- xmlGetAttr(gref, "GraphId")
  nam <- xmlGetAttr(gref, "TextLabel")
  typ <- xmlGetAttr(gref, "Type")
  if (typ == "Complex") {
    if (gid %in% rownames(nodeInfo)) next
    ## create an edge to represent the complex
    edgeCounter <- edgeCounter + 1
    nedg <- data.frame(Source = gid, Target = nam, Type = "represents")
    cat("Complex!\n", file = stderr())
    print(nedg)
    rownames(nedg) <- paste("ec", edgeCounter, sep = "")
    edgeList <- rbind(edgeList, nedg)
  } else if (typ == "GeneProduct") {
    ## create a "contained" edge
    edgeCounter <- edgeCounter + 1
    nedg <- data.frame(Source = gid, Target = grf, MIM = "contained")
    rownames(nedg) <- paste("ec", edgeCounter, sep = "")
    edgeList <- rbind(edgeList, nedg)
  } else {
    stop(" Don't panic; just grab a towel.\n")
  }
}

groups <- getNodeSet(xmlRoot(foo), "/sm:Pathway/sm:Group", rasp)
length(groups) # should be 4
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
  nodeInfo <- rbind(nodeInfo, newn)
  query <- paste("/sm:Pathway/sm:DataNode[@GroupRef='", gid, "' and @Type='Complex']", sep = "")
  repr <- getNodeSet(xmlRoot(foo), query, rasp)
  if (length(repr) > 0) {
    ## create an edge to represent the complex
    edgeCounter <- edgeCounter + 1
    nedg <- data.frame(Source = xmlGetAttr(repr[[1]], "GraphId"), Target = gid, MIM = "represents")
    cat("Complex!\n", file = stderr())
    print(nedg)
    rownames(nedg) <- paste("ec", edgeCounter, sep = "")
    edgeList <- rbind(edgeList, nedg)
  }
}

## Anchors

anchors <- getNodeSet(xmlRoot(foo), "/sm:Pathway/sm:Interaction/sm:Graphics/sm:Anchor", rasp)
acount <- 0
for (anchor in anchors) {
  acount <- acount + 1
  lbl <- paste("Anchor", acount, sep = "")
  gfx <- xmlParent(anchor) # must be a Graphics object
  edge <- xmlParent(gfx)   # must be an Interactiomn (edge) object
  ## Create a node of type "edge" so that the edge that is pointed to
  ## gets transformed to A -> EDGE -> B so that we can make a node for
  ## the target of the one we are working on.
  gid <- xmlGetAttr(edge, "GraphId")
  aid <- xmlGetAttr(anchor, "GraphId")
  edgenode <- data.frame(GraphId = aid,
                         label = lbl,
                         Type = "EDGE")
  pts <- getNodeSet(edge, "./sm:Graphics/sm:Point", rasp)
  src <- xmlGetAttr(pts[[1]], "GraphRef")
  tgt <- xmlGetAttr(pts[[2]], "GraphRef")
  ## Create two edges; src -> EDGE, EDGE -> tgt
  newedges <- data.frame(Source = c(src, aid),
                         Target = c(aid, tgt),
                         MIM    = c("Source", xmlGetAttr(pts[[2]], "ArrowHead")))
  nodeInfo <- rbind(nodeInfo, edgenode)
  edgeList <- rbind(edgeList, newedges)
}
nodeInfo <- as.data.frame(nodeInfo)
nodeInfo$Type <- factor(nodeInfo$Type)
summary(nodeInfo)
edgeList <- as.data.frame(edgeList)
edgeList$MIM <- factor(edgeList$MIM)
summary(edgeList)

library(Polychrome)
data(Light24)
punk <- Light24[c(24, 4, 19, 16, 18, 6, 10)]
names(punk) <- levels(nodeInfo$Type)
swatch(punk)
nodeInfo$color <- punk[nodeInfo$Type]
## set diferent shapes?
shap <- c("circle", "circle", "rectangle", "circle",
          "rectangle", "circle", "rectangle")
names(shap) <- levels(nodeInfo$Type)
nodeInfo$shape <- shap[nodeInfo$Type]
atypes <- c(Arrow = "solid",
            contained = "dotted",
            "mlm-inhibition" = "dashed",
            represents = "twodash",
            Source = "dotdash")
acols <- c(Arrow = "forestgreen",
           contained = "black",
           "mlm-inhibition" = "red",
           represents = "darkcyan",
           Source = "purple")
data.frame(atypes, acols)
edgeList$lty <- atypes[edgeList$MIM]
edgeList$color <- acols[edgeList$MIM]

ttl <- xmlGetAttr(xmlRoot(foo), "Name")

library(igraph)
G <- graph_from_data_frame(edgeList, directed = TRUE, vertices = nodeInfo)
## might want to use set_graph_attr to associate things like the
## name/title with the graph as a whole
E(G)$width = 3
set.seed(12345)
L <- layout_with_graphopt(G)
windows(width = 12, height=12)
opar <- par(mai = c(0.1, 0.1, 1, 0.1), bg = "white")
plot(G, layout = L)
legend("bottomleft", names(atypes), lwd = 4, lty = atypes, col = acols)
legend("topleft", names(punk), pch = 15, col = punk)
title(ttl)
par(opar)
resn <- 300
dev.copy(png, file = "IGF1.png", width = 12*resn, height = 12*resn, res = resn, bg = "white")
dev.off()

save(G, punk, atypes, acols, file = "wp3850.Rda")

