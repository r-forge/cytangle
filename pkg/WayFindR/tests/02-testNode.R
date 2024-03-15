library(WayFindR)
## input directly from file
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
nodes <- collectNodes(xmlfile)
class(nodes)
dim(nodes)
head(nodes)
tail(nodes)
## Alternative: input first, then look for nodes
doc <- XML::xmlParseDoc(xmlfile)
wnodes <- collectNodes(doc)
all( wnodes == nodes )
## Obvious failure modes
try( collectNodes() )
try( collectNodes("bogus") )
badfile <- system.file("perl/nodeTypes.txt", package = "WayFindR")
try( collectNodes(badfile) )
