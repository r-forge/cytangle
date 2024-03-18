library(WayFindR)
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
links <- collectAnchors(xmlfile)
class(links)
links$nodes
links$edges
