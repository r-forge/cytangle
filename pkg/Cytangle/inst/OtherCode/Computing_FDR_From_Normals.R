###############################################################################################################
###############################################################################################################
init_thresh <- 0.99
tol <- 1e-2
tube = "B"
s0 <- 0.03267334
###############################################################################################################
###############################################################################################################

load("~/Desktop/CyTOF_Analysis/code/results/AML_info.Rda")
load("~/Desktop/CyTOF_Analysis/code/results/CellPerNode91216.Rda")
source('~/Desktop/CyTOF_Analysis/code/util/CyTOF_fcns.R')

num_subtypes <- table(groups[1:58,"Group Number"])
marker_names <- get(paste0("tube",tube,"_panel"))[featInd-2]

#bonferroni correction of threshold to account for the number of biomarker pairs
num_feat <- length(featInd)
num_pairs <- choose(num_feat,2)
desired_thresh <- init_thresh + (1-1/num_pairs)*(1-init_thresh)

# find statistic values for normal comparisons
normals <- chi_delta_patients(group1 = 10, tube = tube,s0_val=s0)
seven_normals <- chi_delta_patients(group1 = 10,compar_size = 7,tube=tube,s0_val=s0)
#hist(as.vector(normals),breaks = 100)
#hist(as.vector(seven_normals),breaks = 100)

# determine cutoff for proposed null distribution defined by normals
init_cutoff <- as.numeric(quantile(as.vector(normals),desired_thresh))
# check significant findings for groups of seven to show indenpendence from groupings
seven_on_seven_thresh <- mean(seven_normals > init_cutoff)

#if thresholds for normals and seven vs. seven normals agree, proceed to normals vs. AML
if (abs((1-desired_thresh) - seven_on_seven_thresh) < tol){ # proceed if significance is indep. of group size
  group_chi_delta_stats <-chi_delta_groups_vs_normals(tube = tube,s0_val=s0)
  inds <- which(group_chi_delta_stats > init_cutoff, arr.ind = TRUE)
  stat_table <- cbind(rownames(group_chi_delta_stats)[inds[,1]],colnames(group_chi_delta_stats)[inds[,2]],group_chi_delta_stats[group_chi_delta_stats > init_cutoff])
  colnames(stat_table) <- c("Group Comparison","Markers","Chi Delta statistic")
  
  TP_and_FP <-length(stat_table[,1])/(num_clusts*(num_grps-1))
  TP <- 1-init_thresh
  FDR <- (TP_and_FP-TP)/TP_and_FP
} # tolerance check
  


