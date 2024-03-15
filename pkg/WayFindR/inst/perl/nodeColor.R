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
save(nodeShapes, nodeColors, nodeTypes, file = "../../data/nodeStyle.rda")
