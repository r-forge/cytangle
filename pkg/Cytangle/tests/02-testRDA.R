library(Cytangle)

try( readTangleRDA("pat1") )  # fail; no directory
dname <- system.file("Design", package="Cytangle")
try( readTangleRDA(dname) )  # fail: no files

if (FALSE) {
  mypath <- "C:/KRC/Trainees/RBMcGee/clusters91216/A"
  aml4 <- readTangleRDA("AML4", 2, mypath) # succeed
  dim(aml4) # 169 x 56
  aml4 <- readTangleRDA("AML4", 1:10, mypath) # succeed
  dim(aml4) # 524 x 56
  aml4 <- readTangleRDA("AML4",  path=mypath) # succeed, but slowly
  dim(aml4) # 89035 x 56
}

if (FALSE) {
  library(flowCore)
  mypath <- "C:/KRC/Trainees/RBMcGee/RawFCS"
  fcs <- "TA-AMLx3_AML4_MysteryRemoved.fcs"
  aml4.fcs <- read.FCS(file.path(mypath, fcs), transformation=FALSE)
  class(aml4.fcs)

  foo <- as.numeric(apply(exprs(aml4.fcs)[,1:2], 1, paste, collapse='.'))
  goo <- as.numeric(apply(aml4[,1:2], 1, paste, collapse='.'))
  ox <- order(foo)
  rx <- rank(foo)
  oy <- order(goo)
  ry <- rank(goo)

  plot(foo[ox], goo[oy])
  summary(foo[ox] - goo[oy])
  summary(foo - foo[ox][rx])
  plot(foo, goo[oy][rx])

  plot(exprs(aml4.fcs)[,7], aml4[,7][oy][rx])

}
