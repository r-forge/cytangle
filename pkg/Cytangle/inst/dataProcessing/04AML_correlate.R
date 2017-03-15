# AML cluster correlation script
# Created before 6/8/16 by R. McGee 
# Updated 10/4/16 by R. McGee 

# Set working directory (and location of .fcs data files), load libraries

.libPaths("~/Software/RLib")

#setwd("~/Desktop/CyTOF_Analysis/code/experiments") 
if (!grepl(paste0(getwd(),"$"),"/home/mcgee.278/Dropbox/Research/CyTOF_Analysis/code/experiments",fixed = TRUE)){
  setwd("~/Dropbox/Research/CyTOF_Analysis/code/experiments")
}

load("../results/AML_info.rda")

source('../util/CyTOF_fcns.R')

markerInd <- featInd 

## Compute correlation and distance matrices

# initialize storage
num_feat <- length(markerInd)

eye26 <- diag(num_feat)

count <- matrix(0,1,num_grps)

# computation loop

for (k in 1){
  if (k == 1){
    tube = 'A'
  }else{
    tube = 'B'
  }
  
  tube_path <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube)
  
  all_path <- paste0(tube_path,"/all_cms")
  if (!dir.exists(all_path)){
    dir.create(all_path)
  }
  
  grp_path <- paste0(tube_path,"/grp_cms")
  if (!dir.exists(grp_path)){
    dir.create(grp_path)
  }
  
  for (g in 1:num_grps){
    indiv_grp_path <- paste0(grp_path,'/group',g)
    if (!dir.exists(indiv_grp_path)){
      dir.create(indiv_grp_path)
    } 
  }
  
  for (clust_idx in 1:num_clusts){
    for (file_idx in 1:length(groups[,1])){
      if (groups[file_idx,4] == tube){
        #analyze data frames patient by patient 
        cm_path <- paste0(tube_path,"/",groups[file_idx,1],"/cms")
        if (!dir.exists(cm_path)){
          dir.create(cm_path)
        }
        cur_file <- paste0(tube_path,"/",groups[file_idx,1],"/node",clust_idx,".Rda")
        load(cur_file) #tempMat (patient "file_idx" info at node "clust_idx")
        if (!is.matrix(tempMat)){ #for size 1 x num_feat frames where dim() won't work
          cur_frame <- t(as.matrix(tempMat))[,markerInd,drop=FALSE] # drop forces this to stay a matrix
        }else{
          cur_frame <- tempMat[,markerInd]
        }
        rm(tempMat)
        
        assign(paste0('cm_',clust_idx),correlate(cur_frame,num_feat)) #clust 2, file 12
        save(list=paste0('cm_',clust_idx),file=paste0(cm_path,'/cm_',clust_idx,".Rda"))
        
        # Build all-patient data frame (this gathers all of the patient data at a given node)
        if (!exists('all_temp')){ #(file_idx == 1){
          all_temp <- cur_frame
        }else{
          all_temp <- rbind(all_temp,cur_frame)
        }
        
        # build group by group data frame (this builds a frame with all of the data for each group at a given node)
        
        if(groups[file_idx,2] == 1){
          if (count[1,1] == 0){
            grp_temp1 <- cur_frame 
            count[1,1] <- 1
          }else{
            grp_temp1 <- rbind(grp_temp1,cur_frame)
          }
        }else if(groups[file_idx,2] == 2){
          if (count[1,2] == 0){
            grp_temp2 <- cur_frame
            count[1,2] <- 1
          }else{
            grp_temp2 <- rbind(grp_temp2,cur_frame)
          }
        }else if(groups[file_idx,2] == 3){
          if (count[1,3] == 0){
            grp_temp3 <- cur_frame
            count[1,3] <- 1
          }else{
            grp_temp3 <- rbind(grp_temp3,cur_frame)
          }
        }else if(groups[file_idx,2] == 4){
          if (count[1,4] == 0){
            grp_temp4 <- cur_frame
            count[1,4] <- 1
          }else{
            grp_temp4 <- rbind(grp_temp4,cur_frame)
          }
        }else if(groups[file_idx,2] == 5){
          if (count[1,5] == 0){
            grp_temp5 <- cur_frame
            count[1,5] <- 1
          }else{
            grp_temp5 <- rbind(grp_temp5,cur_frame)
          }
        }else if(groups[file_idx,2] == 6){
          if (count[1,6] == 0){
            grp_temp6 <- cur_frame
            count[1,6] <- 1
          }else{
            grp_temp6 <- rbind(grp_temp6,cur_frame)
          }
        }else if(groups[file_idx,2] == 7){
          if (count[1,7] == 0){
            grp_temp7 <- cur_frame
            count[1,7] <- 1
          }else{
            grp_temp7 <- rbind(grp_temp7,cur_frame)
          }
        }else if(groups[file_idx,2] == 8){
          if (count[1,8] == 0){
            grp_temp8 <- cur_frame
            count[1,8] <- 1
          }else{
            grp_temp8 <- rbind(grp_temp8,cur_frame)
          }
        }else if(groups[file_idx,2] == 9){
          if (count[1,9] == 0){
            grp_temp9 <- cur_frame
            count[1,9] <- 1
          }else{
            grp_temp9 <- rbind(grp_temp9,cur_frame)
          }
        }else if(groups[file_idx,2] == 10){
          if (count[1,10] == 0){
            count[1,10] <- 1
            grp_temp10 <- cur_frame
          }else{
            grp_temp10 <- rbind(grp_temp10,cur_frame)
          }
        }
        rm(cur_frame) 
      } # tube check
    } # patient by patient for loop
    
    # analyze all-patient data frame for node "clust_idx"
    assign(paste0('all_cm_',clust_idx),correlate(all_temp,num_feat))
    save(list=paste0('all_cm_',clust_idx),file=paste0(all_path,'/all_cm_',clust_idx,".Rda"))
    rm(all_temp)
    
    # analyze group data frame for node "clust_idx"
    for (grp_idx in 1:num_grps){ #move was done since every group tempMat may not be created at the time of clust_idx
      indiv_grp_path <- paste0(grp_path,'/group',grp_idx) #paste0(grp_path,grp_idx)
      if (count[1,grp_idx]){
        assign(paste0('grp',grp_idx,'_cm_',clust_idx),correlate(get(paste0("grp_temp",grp_idx)),num_feat))
      }else{
        assign(paste0('grp',grp_idx,'_cm_',clust_idx),eye26)
      }
      save(list=paste0('grp',grp_idx,'_cm_',clust_idx),file=paste0(indiv_grp_path,'/grp',grp_idx,'_cm_',clust_idx,".Rda"))
    } #group analysis for loop
    #reset group info/storage for the next node
    rm(list = ls(pattern = 'grp_temp'))
    count <- matrix(0,1,num_grps)
  } # cluster by cluster for loop
}
rm(list = ls(pattern = '(*|)cm_'))

  
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
