arrowTypes <- read.table("arrowTypes.txt")
arrowTypes = arrowTypes$V1
simplifyArrows <- function(V) {
  remap <- c(Arrow = "mim-stimulation",
             LigandRound = NA,
             LigandSquare = NA,
             Receptor = NA,
             ReceptorRound = NA,
             ReceptorSquare = NA,
             "SBGN-Catalysis" = "mim-catalysis",
             "SBGN-Inhibition" = "mim-inhibition",
             "SBGN-Production" = "mim-stimulation",
             TBar = "mim-inhibition")
  sapply(V, function(X) {
    ifelse(X %in% names(remap), remap[X], X)
  })
}
simpleEdges <- c("contained", "represents",
                 sort(unique(simplify(arrowTypes))))
edgeTypes <- c("dotted", "twodash",
           "solid", "dotdash",
           "dotdash", "dotdash",
           "dotted", "dashed",
           "solid", "dotdash",
           "solid", "dotted",
           "dashed", "solid",
           "twodash", "solid")
edgeColors <- c("gray", "gray",
           "cornflowerblue", "purple",
           "orange", "forestgreen",
           "brown", "magenta",
           "cyan", "black",
           "red", "black",
           "forestgreen", "forestgreen",
           "cornflowerblue", "brown")
names(edgeTypes) <- names(edgeColors) <- simpleEdges
data.frame(simpleEdges, edgeTypes, edgeColors)
plot(1, 1, type = "n")
legend("center", simpleEdges, lwd=4, lty = edgeTypes, col = edgeColors, cex=2)

save(simplifyArrows, edgeTypes, edgeColors, file = "../../data/edgeStyle.rda")
