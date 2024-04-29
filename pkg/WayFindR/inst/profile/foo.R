savePMNG <- FALSE
library(WayFindR)
library(igraph)
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
G <- GPMLtoIgraph(xmlfile)
wc <- which(V(G)$shape == "circle")
G <- set_vertex_attr(G, "shape", index = wc, value = "ellipse")
windows(12, 12)
opar <- par(mai = c(0.05, 0.05, 1, 0.05))
plot(0,0, type = "n")
sz <- (strwidth(V(G)$label) + strwidth("oo")) * 92
G <- set_vertex_attr(G, "size", value = sz)
G <- set_vertex_attr(G, "size2", value = strheight("I") * 2 * 92)
set.seed(12345)
L <- layout_nicely(G)
L2 <- layout_with_kk(G, coords=L)
plot(G, layout = L2)
title("Two-step layout algorithm")
edgeLegend("bottomleft", G)
nodeLegend("bottomright", G)


if (savePNG) {
  resn = 300
  png(file = "igf-layout.png", width = 14*resn, height = 14*resn,
      res = resn, bg = "white")
  plot(G, layout= L2)
  title("Two-step layout algorithm")
  edgeLegend("bottomleft", G)
  nodeLegend("bottomright", G)
  dev.off()
}

if (!requireNamespace("graph")) {
  BiocManager::install("graph", update = FALSE, ask = FALSE)
}
library(graph)
if (!requireNamespace("Rgraphviz")) {
  BiocManager::install("Rgraphviz", update = FALSE, ask = FALSE)
}

library(Rgraphviz)
GN <- as_graphnel(G)
RA <- agopenSimple(GN, "fuckoff")

source("11-graphNEL.R")
gun <- as.graphNEL(G)
all(gun@nodes == GN@nodes)
identical(gun@edgeL, GN@edgeL)
identical(gun@nodeData, GN@nodeData)
identical(gun@edgeData, GN@edgeData)
identical(gun@graphData, GN@graphData)
identical(gun@renderInfo, GN@renderInfo)
identical(gun@renderInfo@pars, GN@renderInfo@pars)
identical(gun@renderInfo@nodes, GN@renderInfo@nodes)
identical(gun@renderInfo@edges, GN@renderInfo@edges)
identical(gun@renderInfo@graph, GN@renderInfo@graph)

GN@renderInfo
length(gun@renderInfo@nodes)
length(GN@renderInfo@nodes)


plot(GN)
plot(gun)


nms <- vertex_attr(G, "name")
lbl <- vertex_attr(G, "label")
shp <- vertex_attr(G, "shape")
col <- vertex_attr(G, "color")
fix <- rep(FALSE, length(nms))
names(lbl) <- names(shp) <- names(col) <- names(fix) <- nms
nAttrs <- list(label = lbl, shape = shp,
               fixedsize = fix, fill = col)
nodeRenderInfo(GN) <- nAttrs
identical(gun@renderInfo@nodes, GN@renderInfo@nodes)

enms <- sub("\\|", "~", as_ids(E(G)))
if (!all(edgeNames(GN) %in% enms)) stop("unknown edge name.")
col <- edge_attr(G, "color")
lty <- edge_attr(G, "lty")
names(col) <- names(lty) <- enms
eAttrs <- list(color = col[edgeNames(GN)],
               style = lty[edgeNames(GN)])
edgeRenderInfo(GN) <- eAttrs
identical(gun@renderInfo@edges, GN@renderInfo@edges)

L <- layoutGraph(GN)

sapply(slotNames(L), function(S) identical(slot(L, S), slot(GN, S)))
RG <- GN@renderInfo
RL <- L@renderInfo
sapply(slotNames(RL), function(S) identical(slot(RL, S), slot(RG, S)))
sapply(names(RG@nodes), function(N) identical(RG@nodes[[N]], RL@nodes[[N]]))

renderGraph(L)

plot(GN, nodeAttrs = nAttrs, edgeAttrs = eAttrs,
     attrs = list(edge = list(lwd = 3)))

plot(gun, nodeAttrs = nAttrs, edgeAttrs = eAttrs,
     attrs = list(edge = list(lwd = 3)))


if (savePNG) {
  resn = 300
  png(file = "igf-graphviz.png", width = 14*resn, height = 14*resn,
      res = resn, bg = "white")
  plot(GN, "dot",  nodeAttrs = nAttrs, edgeAttrs = eAttrs,
       attrs = list(node = list(fixedsize = FALSE),
                    edge = list(lwd = 3)))
  title("Rgraphiz layout algorithm")
  par(lwd=3)
  edgeLegend("bottomleft", G)
  par(lwd = 1, cex = 1.3)
  nodeLegend("topright", G)
  dev.off()
}

plot(0, 0, type = "n") # strwidth doesn't work until plot has been called
set.seed(13579)
L <- igraph::layout_with_graphopt(G)
tkplot(G, canvas.width=1000, canvas.height = 1000)
bonk <- 2.5*(L +200)
bonk[,2] <- 1000 - bonk[,2]
tk_set_coords(8, bonk)
ML <- tk_coords(8)
plot(G, layout=ML, axes=TRUE)





## run test/10-* first
G <- graph_from_data_frame(rbind(reac, rela), vertices = entries)
windows(14, 14)
plot(0, 0, type = "n") # strwidth doesn't work until plot has been called
sz <- (strwidth(V(G)$label) + strwidth("oo")) * 100
G <- set_vertex_attr(G, "shape", value = "rectangle")
G <- set_vertex_attr(G, "size", value = sz)
G <- set_vertex_attr(G, "size2", value = strheight("I") * 2 * 100)
L <- layout_with_kk(G)
plot(G, layout = L)

## Need to figure out best way to truncate long compound names
x <- nchar(V(G)$label)
rev(V(G)$label[order(x)])[1:10]
plot(x ~ factor(V(G)$Type))
summary(gg <- grepl("\\[", V(G)$label))

############################################
library(WayFindR)
library(igraph)
library(RBGL)
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
G <- WayFindR:::GPMLtoIgraph(xmlfile)
fc <- findCycles(G)
GN <- as_graphnel(UG <- as.undirected(G))
boyerMyrvoldPlanarityTest(GN)
pft <- planarFaceTraversal(GN)

subg <- subgraph(G, unique(pft[[1]]))
lg <- layout_nicely(subg)
plot(subg, layout = lg)

nm <- names(V(G))
up <- unique(pft[[1]])
sum(!(nm %in% up))
barf <- nm[!(nm %in% up)]
subn <- subgraph(G, barf)
ln <- layout_nicely(subn)
plot(subn, layout = ln)

mc <- merge_coords(graphs = list(subg, subn),
                   layouts = list(lg, ln))
bigG <- disjoint_union(list(subg, subn))
plot(bigG, layout = mc)

plot(G, layout = mc)
plot(G, layout = layout_nicely)

if (FALSE) { # these both crash R with exit code 5
  pco <- planarCanonicalOrdering(GN)
  CP <- chrobakPayneStraightLineDrawing(GN)
}
