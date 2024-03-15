library(WayFindR)
xmlfile <- system.file("pathways/WP3850.gpml", package = "WayFindR")
groups <- collectGroups(xmlfile)
class(groups)
groups$nodes
groups$edges
