# AML script to create difference matrices
# Created before 9/7/16 by R. McGee 
# Updated 10/4/16 by R. McGee 

# Set working directory (and location of .fcs data files), load libraries

.libPaths("~/Software/RLib")

#setwd("~/Desktop/CyTOF_Analysis/code/experiments") 
if (!grepl(paste0(getwd(),"$"),"/home/mcgee.278/Dropbox/Research/CyTOF_Analysis/code/experiments",fixed = TRUE)){
  setwd("~/Dropbox/Research/CyTOF_Analysis/code/experiments")
}

load("../results/AML_info.rda")

source('../util/CyTOF_fcns.R')
source('../util/matrix_metrics.R')

for (k in 1:2){
  if (k == 1){
    tube = 'A'
  }else{
    tube = 'B'
  }
  
  tube_path <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube)
  
  all_path <- paste0(tube_path,"/all_cms")
  grp_path <- paste0(tube_path,"/grp_cms")
  
  Diff_CMD  <- array(0,c(num_clusts,num_clusts,length(groups[,1])))
  Diff_Chol <- array(0,c(num_clusts,num_clusts,length(groups[,1]))) #Diff_CMD 
  Diff_Cond <- array(0,c(num_clusts,num_clusts,length(groups[,1]))) #Diff_CMD 
  All_Diff_CMD  <- matrix(0,num_clusts,num_clusts)
  All_Diff_Chol <- matrix(0,num_clusts,num_clusts)
  All_Diff_Cond <- matrix(0,num_clusts,num_clusts)
  Grp_Diff_CMD  <- array(0,c(num_clusts,num_clusts,num_grps))
  Grp_Diff_Chol <- array(0,c(num_clusts,num_clusts,num_grps)) #Diff_CMD 
  Grp_Diff_Cond <- array(0,c(num_clusts,num_clusts,num_grps)) #Diff_CMD 
  
  for (clust_idx in 1:num_clusts){
    for (file_idx in 1:length(groups[,1])){
      if (groups[file_idx,4] == tube){
        #analyze data frames patient by patient 
        cm_path <- paste0(tube_path,"/",groups[file_idx,1],"/cms")
        
        load(file=paste0(cm_path,'/cm_',clust_idx,".Rda"))
        # compute difference arrays for each patient
        cm_now <- get(paste0("cm_",clust_idx))
        for (j in 1:clust_idx){
          cm_cur <- get(paste0("cm_",j))
          Diff_CMD[j,clust_idx,file_idx]  <- CMD(cm_cur,cm_now) 
          Diff_Chol[j,clust_idx,file_idx] <- CholD(cm_cur,cm_now)
          Diff_Cond[j,clust_idx,file_idx] <- CondD(cm_cur,cm_now)
        }
        rm(cm_now,cm_cur)
        
        load(file=paste0(all_path,'/all_cm_',clust_idx,".Rda"))
        all_cm_now <- get(paste0("all_cm_",clust_idx))
        for (j in 1:clust_idx){
          all_cm_cur <- get(paste0("all_cm_",j))
          All_Diff_CMD[j,clust_idx]  <- CMD(all_cm_cur,all_cm_now) 
          All_Diff_Chol[j,clust_idx] <- CholD(all_cm_cur,all_cm_now)
          All_Diff_Cond[j,clust_idx] <- CondD(all_cm_cur,all_cm_now)
        }
        rm(all_cm_now,all_cm_cur)  
        
        for (grp_idx in 1:num_grps){ #move was done since every group tempMat may not be created at the time of clust_idx
          indiv_grp_path <- paste0(grp_path,'/group',grp_idx) #paste0(grp_path,grp_idx)
          load(file=paste0(indiv_grp_path,'/grp',grp_idx,'_cm_',clust_idx,".Rda"))
          grp_cm_now <- get(paste0("grp",grp_idx,"_cm_",clust_idx))
          for (j in 1:clust_idx){
            grp_cm_cur <- get(paste0("grp",grp_idx,"_cm_",j))
            Grp_Diff_CMD[j,clust_idx,grp_idx]  <- CMD(grp_cm_cur,grp_cm_now)
            Grp_Diff_Chol[j,clust_idx,grp_idx] <- CholD(grp_cm_cur,grp_cm_now)
            Grp_Diff_Cond[j,clust_idx,grp_idx] <- CondD(grp_cm_cur,grp_cm_now)
          }
          rm(grp_cm_now,grp_cm_cur)
        } #group analysis for loop
      }
    }
  }
  
  Diff_CMD <- sym_array(Diff_CMD)
  Diff_Chol <- sym_array(Diff_Chol)
  Diff_Cond <- sym_array(Diff_Cond) 
  save(Diff_Cond,Diff_Chol,Diff_CMD, file = paste0("../results/IndivPatient_Diff_Mats_",tube,run_date,".Rda"))
  
  All_Diff_CMD  <- symmetrize(All_Diff_CMD)
  All_Diff_Chol <- symmetrize(All_Diff_Chol)
  All_Diff_Cond  <- symmetrize(All_Diff_Cond)
  save(All_Diff_Cond,All_Diff_Chol,All_Diff_CMD, file = paste0("../results/AllPatient_Diff_Mats_",tube,run_date,".Rda"))
  
  Grp_Diff_CMD  <- sym_array(Grp_Diff_CMD)
  Grp_Diff_Chol <- sym_array(Grp_Diff_Chol)
  Grp_Diff_Cond  <- sym_array(Grp_Diff_Cond)
  save(Grp_Diff_Cond,Grp_Diff_Chol,Grp_Diff_CMD, file = paste0("../results/Grp_Diff_Mats_",tube,run_date,".Rda"))
  
  CMD_Mats = ls(pattern = 'Diff_CMD')  
  save(list = CMD_Mats, file = paste0("../results/CMD_Diff_Mats_",tube,run_date,".Rda"))
}

