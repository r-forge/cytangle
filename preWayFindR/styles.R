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
simpleEdges <- c("none", "contained", "Source",
                 sort(unique(simplifyArrows(arrowTypes))))
edgeTypes <- c("solid",  "dotted", "twodash",
           "solid", "dotdash",
           "dotdash", "dotdash",
           "dotted", "dashed",
           "solid", "dotdash",
           "solid", "dotted",
           "dashed", "solid",
           "twodash", "solid")
edgeColors <- c("gray", "gray", "gray",
           "cornflowerblue", "gold3",
           "orange", "forestgreen",
           "brown", "magenta",
           "cyan", "black",
           "red", "black",
           "forestgreen", "forestgreen",
           "cornflowerblue", "brown")
names(edgeTypes) <- names(edgeColors) <- simpleEdges
library(Polychrome)
hexed <- gplots::col2hex(unique(edgeColors))
set.seed(95387)
cp <- createPalette(20, hexed, range = c(40, 70))
keggE <- data.frame(name = c("pathway", "compound", "reversible", "irreversible"),
                    type = c("dotted", "dashed", "dashed","solid"),
                    color = c("cyan", "royalblue1", "purple", "purple"))
x <- keggE$type
names(x) <- keggE$name
y <- keggE$color
names(y) <- keggE$name

edgeTypes <- c(edgeTypes, x)
edgeColors <- c(edgeColors, y)
simpleEdges <- c(simpleEdges, keggE$name)
data.frame(simpleEdges, edgeTypes, edgeColors)
plot(1, 1, type = "n")
legend("center", simpleEdges, lwd=4, lty = edgeTypes, col = edgeColors, cex=1.3)

nodeTypes <- read.table("nodeTypes.txt", sep = "\t", header = FALSE)
colnames(nodeTypes) <- c("Tag", "Count")
nodeTypes
data(Light24)
set.seed(84984)
pal <- createPalette(20, Light24[c(16,19)], range = c(70,90))
nodeColors <- pal[c(5, 7, 4, 2, 6, 1, 3, 8, 9, 10, 10, 7, 11, 12)]
nodeColors <- c(colorNames(nodeColors), "gray75", "navajowhite2", "bisque3")
names(nodeColors) <- c(nodeTypes$Tag, "Group", "Undefined",
                       "Shape", "EDGE", "Label", "Ortholog")
swatch(nodeColors)

nodeTypes <- rbind(nodeTypes, c(Tag = "Group", Count = 0))
nodeTypes <- rbind(nodeTypes, c(Tag = "Undefined", Count = 0))
nodeTypes <- rbind(nodeTypes, c(Tag = "Shape", Count = 0))
nodeTypes <- rbind(nodeTypes, c(Tag = "EDGE", Count = 0))
nodeTypes <- rbind(nodeTypes, c(Tag = "Label", Count = 0))
nodeTypes <- rbind(nodeTypes, c(Tag = "Ortholog", Count = 0))
nodeTypes$color <- nodeColors
# should we change circles to ellipses here?
nodeShapes <- c("circle", "circle", "circle", "rectangle", "circle",
                "rectangle", "circle", "circle", "rectangle", "rectangle",
                "rectangle", "circle", "rectangle", "circle", "circle",
                "rectangle", "hexagon")
names(nodeShapes) <- names(nodeColors)
nodeTypes$shape <- nodeShapes

save(nodeShapes, file = "../pkg/WayFindR/data/nodeShapes.rda")
save(nodeColors, file = "../pkg/WayFindR/data/nodeColors.rda")
save(edgeTypes,  file = "../pkg/WayFindR/data/edgeTypes.rda")
save(edgeColors, file = "../pkg/WayFindR/data/edgeColors.rda")

save(nodeShapes, nodeColors,
     edgeTypes, edgeColors,
     file = "../pkg/WayFindR/R/sysdata.rda")


