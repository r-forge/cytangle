library(Cytangle)

try( readTangleRDA("pat1") )  # fail; no directory
dname <- system.file("Design", package="Cytangle")
try( readTangleRDA(dname) )  # fail: no files

if (FALSE) {
  mypath <- "D:/Box Sync/cytof"
  aml4 <- readTangleRDA("AML4", 2, mypath) # succeed
  dim(aml4) # 169 x 56
  aml4 <- readTangleRDA("AML4", 1:10, mypath) # succeed
  dim(aml4) # 1120 x 56
  aml4 <- readTangleRDA("AML4",  path=mypath) # succeed, but slowly
  dim(aml4) # 89035 x 56
}

if (FALSE) {
  library(flowCore)
  mypath <- "D:/Box Sync/cytof"
  fcs <- "TA-AMLx3_AML4_MysteryRemoved.fcs"
  aml4.fcs <- read.FCS(file.path(mypath, fcs), transformation=FALSE)
  class(aml4.fcs)

  N <- 1
  repeat {
    foo <- apply(exprs(aml4.fcs)[, 1:N, drop=FALSE], 1, paste, collapse='|')
    if (! any(duplicated(foo))) break
    N <-  N + 1
  }
  goo <- apply(aml4[,1:N], 1, paste, collapse='|')
  rx <- rank(foo)
  oy <- order(goo)


  J <- sample(51, 1)
  plot(exprs(aml4.fcs)[,J], aml4[,J][oy][rx], main=colnames(aml4)[J])

}
