source("r00-paths.R")
load(file.path(paths$clean, "groups.rda"))
load(file.path(paths$clean, "isotopes.rda"))
ls()

### Tube A
colnames(onA)[5] <- "MolWt"
x <- paste("(", rownames(onA), ")Di", sep="")
onA$Header <- x
onA <- onA[, c(6, 2, 1, 5, 4, 3)]
head(onA)

all(Atarget %in% onA$Short)
x <- rep(NA, nrow(onA))
x[onA$Short %in% Atarget] <- "target"
x[is.na(x) & grepl("^CD", onA$Short)] <- "surface"
x[onA$Short == "HLA-DR"] <- "surface"
x[is.na(x)] <- "other"
onA$Role <- x

write.csv(onA, file="Ainfo.csv")

### Tube B
colnames(onB)[5] <- "MolWt"
x <- paste("(", rownames(onB), ")Di", sep="")
onB$Header <- x
onB <- onB[, c(6, 2, 1, 5, 4, 3)]
head(onB)

all(Btarget %in% onB$Short)
x <- rep(NA, nrow(onB))
x[onB$Short %in% Btarget] <- "target"
x[is.na(x) & grepl("^CD", onB$Short)] <- "surface"
x[onB$Short == "HLA-DR"] <- "surface"
x[is.na(x)] <- "other"
onB$Role <- x

write.csv(onB, file="Binfo.csv")

