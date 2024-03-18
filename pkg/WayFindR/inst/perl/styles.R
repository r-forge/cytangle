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
simpleEdges <- c("contained", "Source",
                 sort(unique(simplifyArrows(arrowTypes))))
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


nodeTypes <- read.table("nodeTypes.txt", sep = "\t", header = FALSE)
colnames(nodeTypes) <- c("Tag", "Count")
nodeTypes

library(Polychrome)
data(Light24)
set.seed(84984)
pal <- createPalette(10, Light24[c(16,19)], range = c(70,90))
nodeColors <- pal[c(5, 7, 4, 2, 6, 1, 3, 8, 9, 10, 10, 7)]
nodeColors <- c(colorNames(nodeColors), "gray75")
names(nodeColors) <- c(nodeTypes$Tag, "Group", "EDGE")
swatch(nodeColors)

nodeTypes <- rbind(nodeTypes, c(Tag = "Group", Count = 0))
nodeTypes <- rbind(nodeTypes, c(Tag = "EDGE", Count = 0))
nodeTypes$color <- nodeColors
nodeShapes <- c("circle", "circle", "circle", "rectangle", "circle",
                "rectangle", "circle", "circle", "rectangle", "rectangle",
                "rectangle", "circle", "circle")
names(nodeShapes) <- names(nodeColors)
nodeTypes$shape <- nodeShapes

save(nodeShapes, file = "../../data/nodeShapes.rda")
save(nodeColors, file = "../../data/nodeColors.rda")
save(edgeTypes,  file = "../../data/edgeTypes.rda")
save(edgeColors, file = "../../data/edgeColors.rda")

save(nodeShapes, nodeColors,
     edgeTypes, edgeColors,
     file = "../../R/sysdata.rda")


