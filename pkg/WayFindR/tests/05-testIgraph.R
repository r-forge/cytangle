library(WayFindR)
data(edgeColors)
edgeColors
data(edgeTypes)
edgeTypes
data(nodeColors)
nodeColors
data(nodeShapes)
nodeShapes
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
G <- WayFindR:::GPMLtoIgraph(xmlfile)
set.seed(13579)
L <- igraph::layout_with_graphopt(G)
plot(G, layout=L)
