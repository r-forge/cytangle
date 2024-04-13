## run test/10-* first
library(igraph)
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
