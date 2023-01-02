## Copyright (C) 2022 Kevin R. Coombes, RB McGee, and Jake Reed

Feature <- function(values, name, colors, meaning, ...) {
  colRamp <- colorRamp2(range(values, na.rm = TRUE),
                        colors = colors, ...)
  new("Feature",
      name = name,
      values = values,
      colRamp = colRamp,
      meaning = meaning)
}

setMethod("plot", c("Feature", "matrix"), function(x, y, pch = 16, ...) {
  plot(y, ...)
  points(x, y, pch = pch)
})

setMethod("points", "Feature", function(x, view, ...) {
  points(view, col = x@colRamp(x@values), ...)
})


colorFeature <- function(cycle, view, feature, featColors, cycleColor = "black", main = "", ...) {
  plot(view, col = featColors[feature], pch = 16, main = main, ...)
  showCycle(cycle, view, cycleColor)
}

prepFeature <- function(feature, dset, degrees, od) {
  G <- c(dset[, feature], dset[, feature])
  G <- G[od]
  lo <- loess(G ~ deg, span = 0.25)
  ox <- order(lo$x)

  col <- colorRamp2(c(min(G), max(G)), c("white", col.ls[2]))
  ptcol <- col(G)
  
  crackle <- (G - min(G) + 0.01)[1:(length(G)/2)]
  probs <- crackle/sum(crackle)
  set.seed(44188)
  N <- 40000
  folderol <- sample(1:length(radians), N, replace = TRUE, prob = probs)
  weighted <- degrees[folderol]
  dd <- density(weighted)
  xy <- as.matrix(usable[folderol, 1:2])
  ng <- 6
  fit <- try(lapply(1:ng, function(K) {
    movMF(xy, k = K)
  }), silent = TRUE)
  if (inherits(fit, "try-error")) fit <- NULL
  list(feature = feature, G = G, lo = lo, ox = ox, 
       ptcol = ptcol, weighted = weighted, dd = dd,
       xy = xy, fit = fit)
}

showFeature <- function(object) {
  pal <- p34[c(1, 15, 4, 5, 24, 9)]
  opar <- par(mfrow = c(2,2))
  on.exit(par(opar))
  attach(object)
  on.exit(detach())
  plot(deg, G, main = feature)
  lines(lo$x[ox], lo$fitted[ox], col = "red", lwd = 2)
  lines(seq(0, 345, 15), GM[, feature], col = "purple", lwd = 2)
  abline(v = c(0, 360))
  
  plot(usable, col = ptcol, pch = 16, main = feature)
  for(i in 1:dim(c19)[1]) {
    lines(c19[i,,1],c19[i,,2], col="black", lwd=2)
  }
  points(centroid[1], centroid[2], col = "black", pch = 16, cex = 1.5)
  
  hist(weighted, breaks = 55, prob = TRUE, main = "Weighted Samples")
  lines(dd, col = "red", lwd = 2)
  
  plot(xy, col = "gray", pch = 16, main = paste("Local Peaks of", feature))
  for (I in 1:length(fit)) points(fit[[I]]$theta, pch = 16, cex = 1.5, col = pal[I])
}
