load("AML10-tubeA-node287.Rda")
load("isotopes.Rda")
AML10.node287 <- tempMat[, Atarget]
rm(AMat, BMat, Atarget, Btarget, clusmark, onA, onB, tempMat)

load("AML10-tubeA-node287-rips.Rda")
library(TDA)
plot(rips[["diagram"]], barcode = TRUE)
AML10.node287.rips <- rips
rm(keeper)

amldist <- dist(AML10.node287)
library(ClassDiscovery)
spca <- SamplePCA(t(AML10.node287))

f <- "testing.rda"
if (file.exists(f)) {
  load(f)
} else {
  set.seed(54227)
  Arip <- ripsDiag(AML10.node287, 2, 5, "euclidean", "Dionysus", TRUE)
  Brip <- ripsDiag(spca@scores[, 1:3], 2, 5, "euclidean", "Dionysus", TRUE)
  Crip <- ripsDiag(amldist, 2, 5, "arbitrary", "Dionysus", TRUE)
  save(Arip, Brip, Crip, file = f)
}
rm(f)

Acyc <- Arip$cycleLocation[[3]]
Bcyc <- Brip$cycleLocation[[3]]
Ccyc <- Crip$cycleLocation[[3]]

AcycL <- Arip$cycleLocation[[length(Arip$cycleLocation)]]
BcycL <- Brip$cycleLocation[[length(Brip$cycleLocation)]]
CcycL <- Crip$cycleLocation[[length(Crip$cycleLocation)]]


disentangle <- function(cycle, dataset) {
  L <- length(dim(cycle))
  if (L == 2) return(cycle)
  ## Now we have to convert the data-set coordinates to indicies
  N1 <- ncol(dataset)
  N2 <- dim(cycle)[3]
  if (N1 != N2) {
    stop("Mismatched sizes. Dataset = ", N1, " and cycle = ", N2, "\n")
  }
  input <- apply(dataset, 1, paste, collapse = "|")
  dex <- apply(cycle, c(1, 2), function(X) {
    check <- paste(X, collapse = "|")
    which(input == check)
  })
  return(dex)
}

foo <- lapply(Arip$cycleLocation,
              disentangle,
              dataset = AML10.node287)

par(mfrow = c(1, 3))
plot(Arip[["diagram"]], barcode = TRUE)
plot(Brip[["diagram"]], barcode = TRUE)
plot(Crip[["diagram"]], barcode = TRUE)

AML10.node287.rips <- Crip

save(AML10.node287, AML10.node287.rips, Arip, file = "cytof.rda")
rm(AML10.node287, AML10.node287.rips)


load("CR6-tubeA-node108.Rda")
load("isotopes.Rda")
CR6.node108 <- tempMat[, Atarget]
rm(AMat, BMat, Atarget, Btarget, clusmark, onA, onB, tempMat)

load("CR6-Node-108-rips.Rda")
library(TDA)
plot(rips[["diagram"]], barcode = TRUE)
CR6.node108.rips <- rips
rm(rips)
save(CR6.node108, CR6.node108.rips, file = "cytof-CR6.rda")
rm(CR6.node108, CR6.node108.rips)

