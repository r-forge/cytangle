library(WayFindR)
## input directly from file
xmlfile <- system.file("pathways/kegg_hsa00510.xml", package = "WayFindR")
try( entries <- collectEntries(xmlfile) )
