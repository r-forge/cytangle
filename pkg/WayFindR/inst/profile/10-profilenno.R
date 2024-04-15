library(WayFindR)
## input directly from file
xmlfile <- system.file("pathways/kegg_hsa00510.xml", package = "WayFindR")

tic <- Sys.time()
entries0 <- collectEntries(xmlfile, anno = "one")
toc <- Sys.time()
toc - tic # ~ 54 seconds
known <- ls(WayFindR:::WAYcache)
rm(list = known, envir = WayFindR:::WAYcache)

tic <- Sys.time()
entries1 <- collectEntries(xmlfile, anno = "batch")
toc <- Sys.time()
toc - tic # ~ 53 seconds
known <- ls(WayFindR:::WAYcache)
rm(list = known, envir = WayFindR:::WAYcache)

tic <- Sys.time()
entries2 <- collectEntries(xmlfile, anno = "all")
toc <- Sys.time()
toc - tic # ~ 467seconds


if (FALSE) {
  all(entries1 == entries0)
  brow <- sapply(1:nrow(entries),
                 function(I) all(entries[I,] == entries0[I,]))
  entries[!brow,]
  entries0[!brow,]
  library(XML)
  library(KEGGREST)
  library(org.Hs.eg.db)
  source("R/10-keggUtil.R")
  xmldoc = xmlfile
  anno = "all"
}
