#require(rgl)
#require(ClassDiscovery)
## needs Btarget and Atarget from "isotopes.Rda"
## also needs pcdimA and pcdimB from 
voidPlot <- function(tube, S, N, ...) {
  fname <- paste(paste("tube", tube, sep=''), S, N, "rips.Rda", sep='-')
  fp <- file.path(ripDir, tube, fname)
  cat(fp, "\n", file=stderr())
  load(fp)

  load(file.path(paths$clusters, tube, S, paste(N, "Rda", sep='.')))
  M <- tempMat[, get(paste(tube, "target", sep=''))]
  spca <- SamplePCA(t(M))
  usable <- spca@scores[, 1:3]

  RD <- rips$diagram
  dure <- RD[,3] - RD[,2]
  dimn <- RD[,1]
  CY <- rips$cycleLocation
  ptr <- which(dure > 0.25 & dimn == 2)
  if (length(ptr > 1)) {
    i <- which.max(dure[ptr])
    ptr <- ptr[i]
  }
  dure[ptr]
  dimn[ptr]
  cyc <- CY[[ptr]]
  rgl.open()
  par3d(windowRect=c(100, 100, 900, 900))
  rgl.bg(color="white")
  axes3d(col='black')
  rgl.viewpoint(theta=-40, phi=35)
  spheres3d(usable[,1], usable[,2], usable[,3], radius=0.1,
            alpha=1, shininess=50, col="#7777ff")
  for (J in 1:dim(cyc)[1]) {
    spheres3d(cyc[J,,], color="#cc2222", radius=0.15)
    lines3d(cyc[J,,], color="#cc2222", radius=0.5, lwd=3)
  }
}

exprvoidplot <- function(tube, S, N, H, CHOP=4, ...) {
  ## load rips results with cycles
  fname <- paste(paste("tube", tube, sep=''), S, N, "rips.Rda", sep='-')
  fp <- file.path(ripDir, tube, fname)
  cat(fp, "\n", file=stderr())
  load(fp)
  ## load expression data for individual cells in this sample-node
  load(file.path(paths$clusters, tube, S, paste(N, "Rda", sep='.')))
  ## principal component analysis
  M <- tempMat[, get(paste(tube, "target", sep=''))]
  spca <- SamplePCA(t(M))
  usable <- spca@scores[, 1:3]
  ## extract the cycles
  RD <- rips$diagram
  dure <- RD[,3] - RD[,2]
  dimn <- RD[,1]
  CY <- rips$cycleLocation
  ptr <- which(dure > 0.25 & dimn == 2)
  if (length(ptr > 1)) {
    i <- which.max(dure[ptr])
    ptr <- ptr[i]
  }
  dure[ptr]
  dimn[ptr]
  cyc <- CY[[ptr]]
  ## match points in cycles to PC locations
  loco <- apply(cyc, 1:2, function(x) {
    foo <- apply(sweep(usable, 2, x, "-")^2, 1, sum)
    which.min(foo)
  })
  ## standardize expression data and convert to colors
  X <- scale(M[,H])
  MX <- max(abs(X))
  if (MX > CHOP) {
    MX <- CHOP
    X[X > CHOP] <- CHOP
    X[X < -CHOP] <- -CHOP
  }
  ## should the palette be an argument to the function?
  pal <- colorRampPalette(c("cyan", "gray70", "magenta"))(31)
  XC <- pal[cut(X, breaks = seq(-MX, MX, length=31), 
                labels=FALSE, include.lowest=TRUE)]
  ## make the 3D plot
  open3d(windowRect=c(100, 100, 900, 900))
#  rgl.bg(color="white")
  axes3d(col='black')
  rgl.viewpoint(theta=-40, phi=35)
  spheres3d(usable[,1], usable[,2], usable[,3], radius=0.1,
            alpha=1, shininess=90, col=XC, ...)
  for (J in 1:dim(cyc)[1]) {
    for(K in 1:3) {
      spheres3d(cyc[J,K,], color=pal[loco[J,K]], radius=0.2, ...)
    }
    lines3d(cyc[J,,], color="black", radius=0.5, lwd=3)
  }
  title3d(H, line=6, cex=2)
}
