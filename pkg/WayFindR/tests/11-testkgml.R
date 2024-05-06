library(WayFindR)
## input directly from file
xmlfile <- system.file("pathways/kegg_hsa00510.xml", package = "WayFindR")
graf <- KGMLtoIgraph(xmlfile)
set.seed(24086)
L <- igraph::layout_with_kk(graf)
plot(graf, layout = L)
## FAILS becassue of multiple edges:
## RA <- as.graphNEL(graf)
