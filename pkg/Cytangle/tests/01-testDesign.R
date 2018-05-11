library(Cytangle)
library(Biobase)

### Read the experimental design data
fname <- system.file("Design", "Ainfo.csv", package="Cytangle")
ainfo <- read.csv(fname, row.names=1)
fname <- system.file("Design", "Binfo.csv", package="Cytangle")
binfo <- read.csv(fname, row.names=1)
fname <- system.file("Design", "meta.csv", package="Cytangle")
vmeta <- read.csv(fname, row.names=1)

### Set up the annotated data frames
aobj <- AnnotatedDataFrame(ainfo, varMetadata = vmeta)
bobj <- AnnotatedDataFrame(binfo, varMetadata = vmeta)
rm(fname, ainfo, binfo, vmeta)

### Read the sample information
fname <- system.file("Design", "samples.csv", package="Cytangle")
sinfo <- read.csv(fname, row.names = 1)

smeta <- data.frame(labelDescription = c(
  Identifier = "Patient identifier",
  Group.Number = "Numeric code for the group to which the sample belongs",
  AML.Subgroup = "Name of the group to which the sample belongs",
  FileName = "Sample identifier; there may be multiple samples per patient"
))

sobj <- AnnotatedDataFrame(sinfo, varMetadata = smeta)
rm(fname, sinfo, smeta)

fname <- system.file("Design", "miame.csv", package="Cytangle")
mi <- read.csv(fname, row.names = 1, as.is=TRUE)
x <- as.list(mi[,1])
names(x) <- rownames(mi)
x[["samples"]] <- list(x[["samples"]])
miame <- do.call(MIAME, x)
rm(x, mi, fname)

cyd <- new("CytangleDesign",
           experimentData = miame,
           phenoData = sobj,
           assays = list(TubeA = aobj,
                              TubeB = bobj))
rm(miame, sobj, aobj, bobj)

