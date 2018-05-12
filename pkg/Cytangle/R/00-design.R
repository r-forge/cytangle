### Copyright 2018, Kevin R. Coombes and R. B. McGee.
###
### Defines classes describing the design of experiment handled by the Cytangle package.

### Completely provisional; idea is to list the things we want to bundle
### so we can later figure out how to represent each one.
setClass("CytangleDesign",
         slots=c(
           experimentData = "MIAME",         # about the whole experiment
           phenoData = "AnnotatedDataFrame", # sample information
           assays = "list"                   # list of AnnotatedDataFrames
         ))

### assays must be a list of AnnotatedDataFrames
setValidity("CytangleDesign", function(object) {
  all(sapply(object@assays, inherits, what="AnnotatedDataFrame"))
})

CytangleDesign <- function(experimentData, phenoData, assays) {
  new("CytangleDesign",
      experimentData = experimentData,
      phenoData = phenoData,
      assays = assays)
}

setGeneric("dim")

setMethod("dim", "CytangleDesign", function(x) {
  val <- c(nrow(x@assays[[1]]), nrow(x@phenoData), length(x@assays))
  names(val) <- c("antibodies", "samples", "assays")
  val
})

setMethod("summary", "CytangleDesign", function(object, ...) {
  
})
