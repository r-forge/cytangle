library(WayFindR)
library(igraph)
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
G <- GPMLtoIgraph(xmlfile)

nids <- c("c0520", "c7b3c", "c90fd", "efc0d",
          "b3356", "b2d78", "b2305")
S <- subgraph(G, nids)
set.seed (54321)
L <- layout_nicely(S)
plot(S, layout = L)
sz <- (strwidth(V(S)$label) + strwidth("oo")) * 150
S <- set_vertex_attr(S, "size", value = sz)
S <- set_vertex_attr(S, "size2", value = strheight("I") * 2 * 150)
plot(S, layout = L)
S <- set_vertex_attr(S, "shape", value = "ellipse")
plot(S, layout = L)
S <- set_vertex_attr(S, "shape", value = "hexagon")
plot(S, layout = L)
