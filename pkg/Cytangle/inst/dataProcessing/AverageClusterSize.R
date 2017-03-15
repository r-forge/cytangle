setwd("~/Desktop/CyTOF_Analysis/code/experiments")

dataset <-'AML_experiment_51184' # Nature_viSNE_data
folder <- paste("../../data/",dataset,"/FCS_edit", sep="")

load("../results/AML_info.rda")

library("flowCore")

# Collect data files
allfiles <- list.files(path = folder, pattern = ".fcs")

# Non-MDS files
MDS_indices <- grep('MDS',allfiles)
if (length(MDS_indices) == 0){
  files <- allfiles
}else{
  files <- allfiles[-MDS_indices]
}

num_subtypes <- table(groups[1:58,"Group Number"])

# Averaging loop
Tot <- 0

countA <- matrix(0,1,num_grps)
countB <- matrix(0,1,num_grps)


for (file_idx in 1:length(groups[,1])){
  # load file
  myFrame <- read.FCS(file.path(folder, files[file_idx]), transformation=FALSE)
  
  # update overall total
  Tot <- Tot + dim(myFrame)[[1]]
  
  # determine patient subtype
  type <- as.numeric(groups[file_idx,2][[1]])
  
  # add to subtype tally for the appropriate tube
  if (groups[file_idx,4] == 'A'){
    countA[1,type] = countA[1,type] + dim(myFrame)[[1]]
  }else{
    countB[1,type] = countB[1,type] + dim(myFrame)[[1]]
  }
}

Avg <- Tot/length(files)


#avg <- ceiling(233718.953846154)