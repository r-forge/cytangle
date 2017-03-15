# AML script to find extreme indices in difference matrices
# Created before 9/6/16 by R. McGee 
# Updated 9/28/16 by R. McGee 

# Set working directory (and location of .fcs data files), load libraries

.libPaths("~/Software/RLib")

#setwd("~/Desktop/CyTOF_Analysis/code/experiments") 
if (!grepl(paste0(getwd(),"$"),"/home/mcgee.278/Dropbox/Research/CyTOF_Analysis/code/experiments",fixed = TRUE)){
  setwd("~/Dropbox/Research/CyTOF_Analysis/code/experiments")
}

load("../results/AML_info.rda")

source('../util/CyTOF_fcns.R')

folder <- "../results"

for (k in 1){
  if (k == 1){
    tube = 'A'
  }else{
    tube = 'B'
  }
  
  clust_path <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube)
  all_path <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,"/all_cms")
  grp_path <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,"/grp_cms")
  
  marker_names <- get(paste0("tube",tube,"_panel"))[featInd-2]
  
  load(paste0(folder,"/CMD_Diff_Mats_",tube,run_date,".Rda"))
  
  Indiv_Array <- Diff_CMD
  Grp_Array <- Grp_Diff_CMD
  All_Matrix <- All_Diff_CMD
  
  cor_coeffs <- matrix(-2,1+num_grps,num_clusts)
  
  # find indices where extrema occur
  Indiv_CMD_inds_orig <- apply(Indiv_Array,3,find_inds)
  Grp_CMD_inds_orig <- apply(Grp_Array,3,find_inds)
  All_CMD_inds <- find_inds(All_Matrix)
  
  Indiv_CMD_inds <- list()
  Grp_CMD_inds <- list()
  
  # check for multiple extrema
  
  mult_extrema <- matrix("-",length(Indiv_CMD_inds_orig),2)
  colnames(mult_extrema) <- c("Individuals","Subgroups")
  
  for(i in 1:length(Indiv_CMD_inds_orig)){
    if(dim(Indiv_CMD_inds_orig[[i]])[1]>2){
      mult_extrema[i,1] = "TRUE"
      avoid <- Mode(Indiv_CMD_inds_orig[[i]])
      aug_inds <- Indiv_CMD_inds_orig[[i]][(Indiv_CMD_inds_orig[[i]][,1]!=avoid),1]
      assign(paste0("augIndiv",i,"_Diff_CMD"),Diff_CMD[-aug_inds,-aug_inds,i])
      Indiv_CMD_inds[[i]] <- find_inds(Diff_CMD[-aug_inds,-aug_inds,i])
    }else{
      mult_extrema[i,1] = "FALSE"
      Indiv_CMD_inds[[i]] <- Indiv_CMD_inds_orig[[i]]
    }
  }
  rm(i)
  
  for(g in 1:length(Grp_CMD_inds_orig)){
    if(dim(Grp_CMD_inds_orig[[g]])[1]>2){
      mult_extrema[g,2] = "TRUE"
      avoid <- Mode(Grp_CMD_inds_orig[[g]])
      aug_inds <- Grp_CMD_inds_orig[[g]][(Grp_CMD_inds_orig[[g]][,1]!=avoid),1]
      assign(paste0("augGrp",g,"_Diff_CMD"),Grp_Diff_CMD[-aug_inds,-aug_inds,g])
      Grp_CMD_inds[[g]] <- find_inds(Grp_Diff_CMD[-aug_inds,-aug_inds,g])
    }else{
      mult_extrema[g,2] = "FALSE"
      Grp_CMD_inds[[g]] <- Grp_CMD_inds_orig[[g]]
    }
  }
  rm(g)
  
  # find distance statistics [1 x length(files)]
  mu_distIndiv <- apply(Indiv_Array,3,mean)
  sd_distIndiv <- apply(Indiv_Array,3,sd)
  coeff_varIndiv <- mu_distIndiv/sd_distIndiv
  
  mu_distGrp <- apply(Grp_Array,3,mean)
  sd_distGrp <- apply(Grp_Array,3,sd)
  coeff_varGrp <- mu_distGrp/sd_distGrp
  
  mu_distAll <- mean(All_Matrix)
  sd_distAll <- sd(All_Matrix)
  coeff_varAll <- mu_distAll/sd_distAll
  
  # extract the maximum change in features
  assign(paste0("delta_",tube,"all"),extract_max(All_CMD_inds,all_path,marker_names))
  
  for (grp_idx in 1:num_grps){ #move was done since every group tempMat may not be created at the time of clust_idx
    indiv_grp_path <- paste0(grp_path,'/group',grp_idx) #paste0(grp_path,grp_idx)
    assign(paste0("delta_",tube,"grp",grp_idx),extract_max(Grp_CMD_inds[[grp_idx]],indiv_grp_path,marker_names))
  }  
  
  save(list = ls(pattern = "aug(G|I)"),file=paste0(folder,"/Augmented_CMD_Diff_Mats_",tube,run_date,".Rda"))
  rm(list = ls(pattern = "aug(G|I)"))
  
  # for (grp_idx in 1:num_grps) {
  #   #cat("Group", g, "\n", file=stderr())
  #   loc <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,"/grp_cms/group",g)
  #   stuff <- array(NA, dim=c(length(featInd),length(featInd),num_clusts))
  #   for (clust_idx in 1:num_clusts) {
  #     bollocks <- paste("grp", g, "_cm_", n, sep='')
  #     fil <- file.path(loc, paste(bollocks, ".Rda", sep=''))
  #     load(fil)
  #     stuff[,,clust_idx] <- get(bollocks)
  #     rm(list=bollocks)
  #   }
  #   dimnames(stuff) <- list(marker_names,marker_names,1:483)
  #   assign(paste("group", g, sep=""), stuff)
  #   rm(stuff, clust_idx, bollocks, fil)
  # }

  # store all correlation coefficients for previously found features
  delts_all <- get(paste0("delta_",tube,"all"))
  delts_all_inds <- which(marker_names %in% delts_all)
    for (clust_idx in 1:num_clusts){
      load(file=paste0(all_path,'/all_cm_',clust_idx,".Rda"))
      cor_coeffs[1,clust_idx] <- get(paste0("all_cm_",clust_idx))[delts_all_inds[1],delts_all_inds[2]]
      for (grp_idx in 1:num_grps){ 
        delts_grp <- get(paste0("delta_",tube,"grp",grp_idx))
        delts_grp_inds <- which(marker_names %in% delts_grp)
        indiv_grp_path <- paste0(grp_path,'/group',grp_idx) 
        load(file=paste0(indiv_grp_path,'/grp',grp_idx,'_cm_',clust_idx,".Rda"))
        cor_coeffs[1+grp_idx,clust_idx] <- get(paste0("grp",grp_idx,"_cm_",clust_idx))[delts_grp_inds[1],delts_grp_inds[2]]
      }
      rm(list = ls(pattern = '_cm'))
    }
  rownames(cor_coeffs) <- c("All-patient","Group 1","Group 2","Group 3","Group 4","Group 5","Group 6","Group 7","Group 8","Group 8","Group 10")
  
  save(list = ls(pattern = "dist|coeff|mult|CMD_inds|delta"),file=paste0(folder,"/AML_stats_",tube,run_date,".Rda"))
  rm(list = ls(pattern = "dist|coeff|mult|CMD_inds|delta"))
}



# for (i in 1:dim(mult_extrema)[1]){
#   if (mult_extrema[i,1]){
#     if(which(Indiv_CMD_inds[[i]] == All_CMD_inds[1,])){  #if/any
#       inds <- All_CMD_inds[1,]
#     }else{
#       inds <- Indiv_CMD_inds[[i]][,]
#     } 
#   }
# }
# 
# for (i in 1:num_grp){
#   if (mult_extrema[i,2]){
#     if(which(Grp_CMD_inds == All_CMD_inds[1,])){  #if/any
#       inds = All_CMD_inds[1,]
#     } 
#   }
# }

