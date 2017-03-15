
load("../results/AML_info.Rda")
 
markerInd <- featInd

for (k in 1:2) {
  if (k == 1) {
    tube = 'A'
  } else{
    tube = 'B'
  }
  
  count <- matrix(0,1,num_grps)
  
  #temp1 <- matrix("-",num_clusts,length(groups[, 1])/2)
  #temp2 <- matrix(0,num_clusts,num_grps)
  tempCounts <- cbind(matrix(NA,num_clusts,length(groups[, 1])/2),matrix(0,num_clusts,num_grps),matrix(0,num_clusts,num_grps))
  
  tube_path <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube)
  
  for (file_idx in 1:length(groups[, 1])){
    if (groups[file_idx, 4] == tube) {
      for (clust_idx in 1:num_clusts) {
        
        # load patient info at node and count number of cells
        cur_file <- paste0(tube_path, "/", groups[file_idx, 1], "/node", clust_idx, ".Rda")
        load(cur_file) #tempMat (patient "file_idx" info at node "clust_idx")
        if (!is.matrix(tempMat)) { #for size 1 x num_feat frames where dim() won't work
          cur_count <- dim(t(as.matrix(tempMat))[, markerInd, drop = FALSE])[1] # drop forces this to stay a matrix
        } else{
          cur_count <- dim(tempMat[, markerInd])[1]
        }
        rm(tempMat,cur_file)
        
        # check to make sure that numbers are stored in right spot of tempCounts, basically file_idx %% 59
        if (tube == 'B'){
          tempCounts[clust_idx,file_idx - 58] <- cur_count
        }else{
          tempCounts[clust_idx,file_idx] <- cur_count
        }
        
        #determine what group a patient is in
        k <- as.numeric(groups[file_idx,2][[1]])
        
        # efficient counting
        tempCounts[clust_idx,length(groups[, 1])/2+k] <- tempCounts[clust_idx,length(groups[, 1])/2+k] + cur_count
        
        # 04 correlation counting
        if(groups[file_idx,2] == 1){
          if (count[1,1] == 0){
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- cur_count 
            count[1,1] <- 1
          }else{
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- tempCounts[clust_idx,length(groups[, 1])/2+10+k] + cur_count
          }
        }else if(groups[file_idx,2] == 2){
          if (count[1,2] == 0){
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- cur_count 
            count[1,2] <- 1
          }else{
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- tempCounts[clust_idx,length(groups[, 1])/2+10+k] + cur_count
          }
        }else if(groups[file_idx,2] == 3){
          if (count[1,3] == 0){
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- cur_count 
            count[1,3] <- 1
          }else{
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- tempCounts[clust_idx,length(groups[, 1])/2+10+k] + cur_count
          }
        }else if(groups[file_idx,2] == 4){
          if (count[1,4] == 0){
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- cur_count 
            count[1,4] <- 1
          }else{
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- tempCounts[clust_idx,length(groups[, 1])/2+10+k] + cur_count
          }
        }else if(groups[file_idx,2] == 5){
          if (count[1,5] == 0){
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- cur_count 
            count[1,5] <- 1
          }else{
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- tempCounts[clust_idx,length(groups[, 1])/2+10+k] + cur_count
          }
        }else if(groups[file_idx,2] == 6){
          if (count[1,6] == 0){
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- cur_count 
            count[1,6] <- 1
          }else{
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- tempCounts[clust_idx,length(groups[, 1])/2+10+k] + cur_count
          }
        }else if(groups[file_idx,2] == 7){
          if (count[1,7] == 0){
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- cur_count 
            count[1,7] <- 1
          }else{
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- tempCounts[clust_idx,length(groups[, 1])/2+10+k] + cur_count
          }
        }else if(groups[file_idx,2] == 8){
          if (count[1,8] == 0){
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- cur_count 
            count[1,8] <- 1
          }else{
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- tempCounts[clust_idx,length(groups[, 1])/2+10+k] + cur_count
          }
        }else if(groups[file_idx,2] == 9){
          if (count[1,9] == 0){
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- cur_count 
            count[1,9] <- 1
          }else{
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- tempCounts[clust_idx,length(groups[, 1])/2+10+k] + cur_count
          }
        }else if(groups[file_idx,2] == 10){
          if (count[1,10] == 0){
            count[1,10] <- 1
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- cur_count 
          }else{
            tempCounts[clust_idx,length(groups[, 1])/2+10+k] <- tempCounts[clust_idx,length(groups[, 1])/2+10+k] + cur_count
          }
        }
        
        rm(cur_count)
      } # cluster by cluster for loop
    }  # tube check
  } # patient by patient for loop
  colnames(tempCounts) <- 
    c(groups[1:58,1],"Group 1","Group 2", "Group 3","Group 4", "Group 5", "Group 6", "Group 7", "Group 8", "Group 9", "Group 10","Group 1 old","Group 2  old", "Group 3 old","Group 4 old", "Group 5 old", "Group 6 old", "Group 7 old", "Group 8 old", "Group 9 old", "Group 10 old")
  assign(paste0("counts", tube), tempCounts)
} #tube by tube for loop

#storage to compare differences between pooling and direct computations of numbers per group
diff_grp <- matrix(NA,2,10)

for (grp_idx in 1:num_grps){
  
  #subset the matrix of counts to find groups
  if (grp_idx == 1){
    vars <- grepl('AML26|AML35|AML5|AML10|AML32',colnames(get(paste0("counts",tube))),ignore.case = TRUE)
  }else if (grp_idx == 2){
    vars <- grepl('APL',colnames(get(paste0("counts",tube))),ignore.case = TRUE)
  }else if (grp_idx == 3){
    vars <- grepl('AML27|AML7|AML9|AML20|AML21|AML23|AML30|AML39|AML37|AML38|AML40',colnames(get(paste0("counts",tube))),ignore.case = TRUE)
  }else if (grp_idx == 4){
    vars <- grepl('AML13|AML36|AML33|AML42|AML14|AML15',colnames(get(paste0("counts",tube))),ignore.case = TRUE)
  }else if (grp_idx == 5){
    vars <- grepl('AML18|AML4$|AML19|AML29|AML41|AML25',colnames(get(paste0("counts",tube))),ignore.case = TRUE)
  }else if (grp_idx == 6){
    vars <- grepl('AML8|AML6|AML31',colnames(get(paste0("counts",tube))),ignore.case = TRUE)
  }else if (grp_idx == 7){
    vars <- grepl('AML12',colnames(get(paste0("counts",tube))),ignore.case = TRUE)
  }else if (grp_idx == 8){
    vars <- grepl('AML22',colnames(get(paste0("counts",tube))),ignore.case = TRUE)
  }else if (grp_idx == 9){
    vars <- grepl('CR',colnames(get(paste0("counts",tube))),ignore.case = TRUE)
  }else if (grp_idx == 10){
    vars <- grepl('Nl',colnames(get(paste0("counts",tube))),ignore.case = TRUE)
  }
  
  # count the cells per cluster in for group (grp_idx)
  subCounts <- get(paste0("counts",tube))[,vars]
  if (is.matrix(subCounts)){
    diff_grp[1,grp_idx] <- max(rowSums(subCounts) - get(paste0("counts",tube))[,58+grp_idx])
    diff_grp[2,grp_idx] <- max(rowSums(subCounts) - get(paste0("counts",tube))[,58+10+grp_idx])
  }else if(is.vector(subCounts)){
    diff_grp[1,grp_idx] <- max(subCounts - get(paste0("counts",tube))[,58+grp_idx])
    diff_grp[2,grp_idx] <- max(subCounts - get(paste0("counts",tube))[,58+10+grp_idx])
  }
  
}

save(countsA,countsB, file = paste0("../results/CellPerNode",run_date,".Rda"))
