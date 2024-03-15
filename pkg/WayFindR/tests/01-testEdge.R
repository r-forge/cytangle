library(WayFindR)
## input directly from file
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
edges <- collectEdges(xmlfile)
class(edges)
dim(edges)
head(edges)
tail(edges)
## Alternative: input first, then look for edges
doc <- XML::xmlParseDoc(xmlfile)
wedges <- collectEdges(doc)
all( wedges == edges)
## Obvious failure modes
try( collectEdges() )
try( collectEdges("bogus") )
badfile <- system.file("perl/nodeTypes.txt", package = "WayFindR")
try( collectEdges(badfile) )
