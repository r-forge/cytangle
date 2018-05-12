library(Cytangle)
library(Biobase)

### Read the experimental design data
fname <- system.file("Design", "TubeA.txt", package="Cytangle")
aobj <- read.AnnotatedDataFrame(fname)
fname <- system.file("Design", "TubeB.txt", package="Cytangle")
bobj <- read.AnnotatedDataFrame(fname)
rm(fname)

### Read the sample information
fname <- system.file("Design", "samples.txt", package="Cytangle")
sobj <- read.AnnotatedDataFrame(fname)
rm(fname)


### Read the experiment (MIAME) information
fname <- system.file("Design", "miame.txt", package="Cytangle")
miame <- readMIAMEfromFile(fname)
rm(fname)

### Finally create the actual CytangleDesign object
cyd <- CytangleDesign(miame, sobj,
                      list(TubeA = aobj, TubeB = bobj))
rm(miame, sobj, aobj, bobj)

### See what we got
cyd

### Check methods
dim(cyd)

writeLines(strwrap(abstract(cyd@experimentData)))

