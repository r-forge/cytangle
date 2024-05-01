## Function for plotting an elliptical node
igraphEllipse <- function(coords, v=NULL, params) {
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

  DescTools::DrawEllipse(x=coords[,1], y=coords[,2],
                         radius.x = vertex.size, radius.y = vertex.size2,
                         col = vertex.color)
}

## Function for plotting a hexagonal node
igraphHexagon <- function(coords, v=NULL, params) {
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
  mapply(coords[,1], coords[,2], vertex.color, vertex.size, vertex.size2,
         FUN = function(x, y, bg, W, H) {
           dx <- W/2
           dy <- H/2
           polygon(x = c(x - dx, x + dx, x + dx + dy, x + dx, x - dx, x- dx - dy),
                   y = c(y - dy, y - dy, y, y + dy, y + dy, y),
                   col = bg)
         })
}

