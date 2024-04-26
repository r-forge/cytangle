as.graphNEL <- function(G) {
  GN <- as_graphnel(G)
  ## prepare to copy vertex/node attributes from igraph to graphNEL,
  ## which supports shape, fixedsize, color, fillcolor, label, fontcolor, fontsize
  nms <- vertex_attr(G, "name")
  lbl <- vertex_attr(G, "label")
  shp <- vertex_attr(G, "shape")
  col <- vertex_attr(G, "color")
  names(lbl) <- names(shp) <- names(col) <- nms
  nAttrs <- list(label = lbl, shape = shp, fixedsize = FALSE, fillcolor = col)
  nodeRenderInfo(GN) <- nAttrs
  ## prepare to copy edge attributes from igraph to graphNEL,
  ## which supports color, style, lwd, and arrowhead (ha!)
  enms <- sub("\\|", "~", as_ids(E(G)))
  if (!all(edgeNames(GN) %in% enms)) {
    stop("graphNEL version contains edges not in igraph version!\n")
  }
  col <- edge_attr(G, "color")
  lty <- edge_attr(G, "lty")
  names(col) <- names(lty) <- enms
  eAttrs <- list(color = col[edgeNames(GN)],
                 style = lty[edgeNames(GN)])
  edgeRenderInfo(GN) <- eAttrs
  ## graph-wide parameters in graphNEL are:
  ## main, col.main, cex.main, sub, col.sub, cex.sub
  attrs <- list(main = paste(graph_attr(G, "PathID"),
                             graph_attr(G, "Name"),
                             sep = ": "))
  graphRenderInfo(GN) <- attrs
  
  GN
}
