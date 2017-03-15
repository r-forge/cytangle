# AML Clustering Script
# Created 4/10/16 by R. McGee
# Fills in missing entries known to cause trouble for all .fcs in specified folder.

# Set working directory (and location of .fcs data files), load libraries

.libPaths("~/Software/RLib")
#.libPaths("~/Downloads/Software/RLib")

#setwd("./Desktop/CyTOF_Analysis/code/experiments") #Not needed since I'm not using relative paths

dataset <- 'AML_experiment_51184' # 'Nature_viSNE' #
folder <- paste("/Scratch/mcgee.278/R/CyTOF_Analysis/data/",dataset,"/FCS_orig", sep="")
outFolder <- paste("/Scratch/mcgee.278/R/CyTOF_Analysis/data/",dataset,"/FCS_edit", sep="")
testFolder <- paste("/Scratch/mcgee.278/R/CyTOF_Analysis/data/",dataset,"/FCS_test", sep="")
folder <- testFolder
outFolder <- testFolder

library("tools")
library("base")
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

# Conversion loop
for (idx in seq_along(files)) {

  base <- file_path_sans_ext(files[idx])
  data_obj <- read.FCS(file.path(folder, files[idx]), transformation=FALSE)

  # Check for empty parameters in data frame
  if (nchar(description(data_obj)$'$DATE') == 0)
  { description(data_obj)$'$DATE' = "3/22/87"
  keyword(data_obj)$'$DATE' = "3/22/87"
  }
  if (nchar(description(data_obj)$'$BTIM') == 0)
  { description(data_obj)$'$BTIM' = "8:00am"
  keyword(data_obj)$'$BTIM' = "8:00am"
  }
  if (nchar(description(data_obj)$'$ETIM') == 0)
  { description(data_obj)$'$ETIM' = "9:00am"
  keyword(data_obj)$'$ETIM' = "9:00am"
  }
  if (is.null(description(data_obj)$'$P1S'))
  { description(data_obj)$'$P1S' = "" #"empty"
  keyword(data_obj)$'$P1S' = "" #"empty"
  }
  # I think the lines below actually corrupts the files
  if (is.null(description(data_obj)$'$P2S'))
  { description(data_obj)$'$P2S' = "" #"empty"
  keyword(data_obj)$'$P2S' = "" #"empty"
  }
  
  # rewrite fcs file if necessary
  rewritten_file <- paste0(outFolder,"/",base,"_edit.fcs")
  write.FCS(data_obj,rewritten_file)
  
  rewritten_file1 <- paste0(outFolder,"/",base,"_editdate.fcs")
  write.FCS(data_obj,rewritten_file1)
  
  rewritten_file2 <- paste0(outFolder,"/",base,"_editend.fcs")
  write.FCS(data_obj,rewritten_file2)
  
  rewritten_file3 <- paste0(outFolder,"/",base,"_editbegin.fcs")
  write.FCS(data_obj,rewritten_file3)
  
  
  
  rm(base,data_obj,rewritten_file)
}

src <- rewritten_file3
src <- file.path(folder, files[idx])
data_obj <- SPADE.read.FCS(src) # 
