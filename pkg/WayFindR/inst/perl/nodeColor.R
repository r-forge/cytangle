nodeTypes <- read.table("nodeTypes.txt", sep = "\t", header = FALSE)
colnames(nodeTypes) <- c("Tag", "Count")
nodeTypes

library(Polychrome)
data(Light24)
set.seed(84984)
pal <- createPalette(10, Light24[c(16,19)], range = c(70,90))

nodeColors <- pal[c(5, 7, 4, 2, 6, 1, 3, 8, 9, 10)]
nodeColors <- colorNames(nodeColors)
names(nodeColors) <- nodeTypes$Tag[1:10]
swatch(nodeColors)

nodeTypes$color <- c(nodeColors, nodeColors[10])
save(nodeColors, nodeTypes, file = "../../data/nodeColor.rda")
