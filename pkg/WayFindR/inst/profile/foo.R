library(WayFindR)
library(igraph)
library(plotrix)
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
G <- GPMLtoIgraph(xmlfile)

## Function for plotting an elliptical node
myellipse <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/150 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }

  vertex.size2 <- 1/150 * params("vertex", "size2")
  if (length(vertex.size2) != 1 && !is.null(v)) {
    vertex.size2 <- vertex.size2[v]
  }

  draw.ellipse(x=coords[,1], y=coords[,2],
    a = vertex.size, b=vertex.size2, col=vertex.color)
}

## Register the shape with igraph
add_shape("ellipse", clip=shapes("rectangle")$clip,
          plot=myellipse)
wc <- which(V(G)$shape == "circle")
G <- set_vertex_attr(G, "shape", index = wc, value = "ellipse")
windows(10, 10)
opar <- par(mai = c(0.02, 0.02, 0.02, 0.02))
plot(0, 0, type = "n") # strwidth doesn't work until plot has been called
sz <- (strwidth(V(G)$label) + strwidth("oo")) * 92
G <- set_vertex_attr(G, "size", value = sz)
G <- set_vertex_attr(G, "size2", value = strheight("I") * 2 * 92)
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
