library(WayFindR)
library(Rgraphviz)
## input directly from file
xmlfile <- system.file("pathways/kegg_hsa00510.xml", package = "WayFindR")
graf <- KGMLtoIgraph(xmlfile)
set.seed(24086)
L <- igraph::layout_with_kk(graf)
plot(graf, layout = L)
## FAILS because of multiple edges:
try( RA <- as.graphNEL(graf) )

xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
G <- GPMLtoIgraph(xmlfile)
set.seed(13579)
GN <- as.graphNEL(G)
plot(GN)
