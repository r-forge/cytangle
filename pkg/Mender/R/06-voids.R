#require(rgl)
#require(ClassDiscovery)
## needs Btarget and Atarget from "isotopes.Rda"
## also needs pcdimA and pcdimB from 
voidPlot <- function(cycle, view, feature = NULL, radius = 0.01, ...) {
  index <- cycle@index
  creak <- unique(rbind(index[, 1:2], index[, 2:3], index[, c(1,3)]))
  left <- t(sapply(creak[,1], getCoords, view = view))
  right <- t(sapply(creak[,2], getCoords, view = view))
  if (is.null(feature)) {
    colors <- "gray"
  } else {
    colors <- feature@colRamp(feature@values)
  }
  open3d(windowRect=c(100, 100, 900, 900))
  rgl.bg(color="white")
  axes3d(col='black')
  rgl.viewpoint(theta=-40, phi=35)
  spheres3d(view[,1], view[,2], view[,3], radius=radius,
            alpha=1, shininess=50, col = colors)
  for (J in 1:nrow(left)) {
    spheres3d(left[J,], color = cycle@color, radius=radius)
    spheres3d(right[J,], color = cycle@color, radius=radius)
    lines3d(rbind(left[J,], right[J,]), color = cycle@color, radius=0.5, lwd=3)
  }
  if (!is.null(feature)) {
    title3d(feature@name, line = 6, cex = 2)
  }
  invisible(cycle)
}

voidFeature <- function(feature, view, radius = 0.01, ...) {
  colors <- feature@colRamp(feature@values)
  spheres3d(view[,1], view[,2], view[,3], radius = radius,
            alpha=1, shininess=90, col = colors, ...)
  title3d(feature@name, line = 6, cex = 2)
  invisible(feature)
}
