# AML cluster analysis script
# Created before 6/8/16 by R. McGee 
# Updated 10/5/16 by R. McGee 

# Set working directory (and location of .fcs data files), load libraries

.libPaths("~/Software/RLib")

#setwd("~/Desktop/CyTOF_Analysis/code/experiments") 
if (!grepl(paste0(getwd(),"$"),"/home/mcgee.278/Dropbox/Research/CyTOF_Analysis/code/experiments",fixed = TRUE)){
  setwd("~/Dropbox/Research/CyTOF_Analysis/code/experiments")
}

# Load AML Info
load("../results/AML_info.rda")

folder <- paste(root,"/CyTOF_Analysis/data/",dataset,"/spade",run_date, sep="")

cluster_path <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date)
if (!dir.exists(cluster_path)){
  dir.create(cluster_path)
}

library("flowCore")

# Collect data files
files <- list.files(path = folder, pattern = ".fcs.density.fcs.cluster.fcs$")

for (k in 1:2){
  if (k == 1){
    tube = 'A'
  }else{
    tube = 'B'
  }
  tube_path <- paste0(cluster_path,"/",tube)
  if (!dir.exists(tube_path)){
    dir.create(tube_path)
  }
  marker_names <- get(paste0("tube",tube,"_panel"))[featInd-2]
  
  for (file_idx in seq_along(files)){ 
    if (groups[file_idx,4] == tube){
      
      src <- file.path(folder, files[file_idx])
      cluster_obj <- read.FCS(src, transformation=FALSE) #SPADE.read.FCS(src)
      temp <- exprs(cluster_obj) # read.csv(outFile)
      # normalize
      cluster_frame <- cbind(temp[,1:2],asinh(0.2*temp[,3:51]),temp[,52:56])
      rm(src,cluster_obj,temp)
      
      cur_path <- paste0(tube_path,"/",groups[file_idx,1])
      if (!dir.exists(cur_path)){
        dir.create(cur_path)
      }
      
      for (clust_idx in 1:num_clusts){
        cur_file <- paste0(cur_path,"/node",clust_idx,".Rda")
        curr_clustInd <- (cluster_frame[,"cluster"] == clust_idx)
        tempMat <- cluster_frame[curr_clustInd,]
        colnames(tempMat) <- marker_names 
        save(tempMat,file=cur_file)    
        rm(tempMat,curr_clustInd,cur_file)
      }
      rm(cur_path,cluster_frame)
    }
  }
}

# Questions
  
# What are the files again?

# Which things are far apart
# Which clusters are those - find indices
# Compare their "cms"


#Get total numbers of cells, other stats
#Doing patient by patient analysis:
#  Consider if the maximums across patient are nearby and find out if there's a region'
#Do the maximums cluster in a region
#  FInd max in one patient, choose one point as a base point and vary the correlation matrices at that point in other patients
