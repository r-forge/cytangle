library(WayFindR)
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
gr <- GPMLtoIgraph(xmlfile)
cyc <- findCycles(gr)
cyc
lapply(cyc, interpretCycle, graph = gr)

S <- WayFindR:::cycleSubgraph(gr, cyc)
#set.seed(13579)
#L <- igraph::layout_with_graphopt(G)
#plot(G, layout=L)

disc <- discoverCycles(gr)
