library(WayFindR)
## input directly from file
xmlfile <- system.file("pathways/kegg_hsa00510.xml", package = "WayFindR")

entries <- collectEntries(xmlfile, "one")
class(entries)
table(entries$Type)
head(entries)

reac <- collectReactions(xmlfile)
dim(reac)
table(reac$MIM) # all irreversible

rela <- collectRelations(xmlfile)
class(rela)
dim(rela)
table(rela$MIM)
