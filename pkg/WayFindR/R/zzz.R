.onLoad <- function(libname, pkgname) {
  ## Register the 'ellipse' shape with igraph
  igraph::add_shape("ellipse",
                    clip = igraph::shapes("rectangle")$clip,
                    plot = igraphEllipse)
  invisible("ellipse")
}
