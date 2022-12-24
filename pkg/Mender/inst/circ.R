## read the angle-sector data
data_dir <- file.path(getwd(), 'data')
results_dir <- file.path(getwd(), 'results')
dirLoop <- file.path(results_dir, 'Loop')

## Load test data
resfile <- file.path(dirLoop, 'angleMeans_combined_result.rds')
load(file = resfile)
sectorInfo <- sapply(angle.ls, function(X) X$Mean)
## why only 23 sectors?
interp <- (sectorInfo[1,] + sectorInfo[23,])/2
sectorInfo <- rbind(sectorInfo, interp)

##########################################
## Try working with BioCircos
if (!require("BioCircos")) {
  install.packages("BioCircos")
  library("BioCircos")
}

fakegenome <- list("chr1" = 360)

W <- 15
starts <- seq(0, 360-W, by = W)
ends <- seq(W, 360, by = W)
tails <- c("#3333DD", "#338888",  "#338833", "#DD3333", "#DD33DD")
names(tails) <- colnames(sectorInfo)
Polychrome::swatch(tails, main = "Legend", cex = 2)

plot(c(0,1), c(0,1), type = "n")
legend("center", names(tails), col = tails, lwd = 4, cex = 2)

whoa <- lapply(1:ncol(sectorInfo), function(I) {
  colset <- c("#DDDDDD", tails[I])
  BioCircosHeatmapTrack(colnames(sectorInfo)[I], c("chr1"),
                        starts, ends, sectorInfo[, I],
                        color = colset,
                        labels = c(rep("", 11), colnames(sectorInfo)[I], rep("", 12)),
                        maxRadius = 1 - 0.1*I,
                        minRadius = 0.94 - 0.1*I)
})
tracks <-  whoa[[1]] + whoa[[2]] + whoa[[3]] + whoa[[4]] + whoa[[5]]


BioCircos(tracks, genome = fakegenome,
          genomeFillColor = c("skyblue", "pink"), chrPad = 0.02)
## doesn't work: legend("topleft", colnames(sectorInfo), col = tails, lwd=4)
## everything gets sent to some mysterious temp directory as an HTML file


plot(c(0,1), c(0,1), type = "n")
legend("center", names(tails), col = tails, lwd = 4, cex = 2)
BioCircos(whoa[[1]], genome = fakegenome,
          genomeFillColor = c("skyblue", "pink"), chrPad = 0.02)


########################

if (!require("circlize")) {
  install.packages("circlize")
  library("circlize")
}



cr <- colorRamp2(seq(0, 1, 0.1), viridisLite::viridis(11))

cmat <- cbind(V1 = vals,
              V2 = c(vals[5:L], vals[1:4]),
              V3 = c(vals[9:L], vals[1:8]))

sectors <- factor(LETTERS[1:24])
circos.initialize(sectors, x = rep(15, 24), xlim = c(-1,10))
spunk <- factor(rep(c("W", "X", "Y", "Z"), each = 6))
split(1:24, spunk)
circos.heatmap(sectorInfo, split = spunk)

#############################################
## Yet anmother package
library(RCircos)

cp <- data.frame(Chromosome = "chr1", chromStart = starts, chromEnd = ends)
dim(cp)
dim(sectorInfo)
hm <- cbind(cp, sectorInfo)
idio <- data.frame(Chromosome = "chr1",
                   chromStart = seq(0, 270, 90),
                   chromEnd = seq(90, 360, 90),
                   Band = c("p12", "p11", "q11", "q12"),
                   Stain = rep(c("gpos100", "gneg"), times = 2))
chr2 <- data.frame(Chromosome = "chr2", chromStart = c(0, 5), chromEnd = c(4, 9),
                   Band = c("p11", "q11"),
                   Stain = c("gneg", "gneg"))
idio = rbind(idio, chr2)

data(UCSC.HG19.Human.CytoBandIdeogram);
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
RCircos.Set.Core.Components(cyto.info = cyto.info,
                            chr.exclude = NULL,
                            tracks.inside = 5,
                            tracks.outside = 0)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()



RCircos.Set.Core.Components(cyto.info = idio,
                            chr.exclude = NULL,
                            tracks.inside = 5,
                            tracks.outside = 0)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()


 RCircos.Label.Chromosome.Names
function (chr.name.pos = NULL) 
{
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    if (is.null(chr.name.pos)) {
        chr.name.pos <- RCircos.Par$chr.name.pos
    }
    else {
        if (chr.name.pos < RCircos.Par$track.in.start) 
            stop("Chromosome name positions overlap with inside track.\n")
        if (chr.name.pos > RCircos.Par$track.out.start) 
            message("May not plot data tracks outside of chromosomes.\n")
    }
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    chroms <- unique(RCircos.Cyto$Chromosome)
    rightSide <- nrow(RCircos.Pos)/2
    for (aChr in seq_len(length(chroms))) {
        chrRows <- which(RCircos.Cyto$Chromosome == chroms[aChr])
        chrStart <- min(RCircos.Cyto$StartPoint[chrRows])
        chrEnd <- max(RCircos.Cyto$EndPoint[chrRows])
        mid <- round((chrEnd - chrStart + 1)/2, digits = 0) + 
            chrStart
        chrName <- sub(pattern = "chr", replacement = "", chroms[aChr])
        text(RCircos.Pos[mid, 1] * RCircos.Par$chr.name.pos, 
            RCircos.Pos[mid, 2] * RCircos.Par$chr.name.pos, label = chrName, 
            pos = ifelse(mid <= rightSide, 4, 2), srt = RCircos.Pos$degree[mid])
    }
}


