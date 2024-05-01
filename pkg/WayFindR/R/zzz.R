.onLoad <- function(libname, pkgname) {
  ## Register the 'ellipse' shape with igraph
  igraph::add_shape("ellipse",
                    clip = igraph::shapes("rectangle")$clip,
                    plot = igraphEllipse)
  igraph::add_shape("hexagon",
                    clip = igraph::shapes("rectangle")$clip,
                    plot = igraphHexagon)
  invisible("shaping")
}

if (FALSE) {
  library(stringi)
  hx <- stri_unescape_unicode("\u2b22")
  plot(1, 1, pch = hx, cex = 4)
  text(1, 1, labels = hx, font = 5)
  text(0.65, 0.65, letters, font = 5, adj = 1)
}
