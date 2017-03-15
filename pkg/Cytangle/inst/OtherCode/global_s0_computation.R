# AML correlation coefficient comparison script
# Created before 1/31/17 by R. McGee 
# Updated 2/16/17 by R. McGee 

# Set working directory (and location of .fcs data files), load libraries

print(Sys.time())
print('Initializing')

.libPaths("~/Software/RLib")

setwd("~/Desktop/CyTOF_Analysis/code/experiments") 
#if (!grepl(paste0(getwd(),"$"),"/home/mcgee.278/Dropbox/Research/CyTOF_Analysis/code/experiments",fixed = TRUE)){
#  setwd("~/Dropbox/Research/CyTOF_Analysis/code/experiments")
#}

###############################################################################################################
###############################################################################################################
tube <- "B"
###############################################################################################################
###############################################################################################################


load("../results/AML_info.Rda")
load("../results/CellPerNode91216.Rda")
source('../util/CyTOF_fcns.R')

marker_names <- get(paste0("tube",tube,"_panel"))[featInd-2]

markerInd <- featInd 
num_feat <- length(markerInd)
num_pairs <- choose(num_feat,2)
num_comparisons <- choose(length(groups[,1]),2)

## function definitions

# fisher z transform
fisherz <- function(x) {
  0.5*(log(1+x) - log(1-x))
}

# s0/exchangeability factor estimation from samr
s0_comp <- function (resid, sd, s0.perc = seq(0, 1, by = 0.05)){
  br = unique(quantile(sd, seq(0, 1, len = 101)))
  nbr = length(br)
  a <- cut(sd, br, labels = F)
  a[is.na(a)] <- 1
  cv.sd <- rep(0, length(s0.perc))
  print(s0.perc)
  for (j in 1:length(s0.perc)) {
    w <- quantile(sd, s0.perc[j])
    w[j == 1] <- 0
    w_vec <- w*rep(1,length(sd))
    tt2 <- resid * 1/(sd + w_vec)#tt * sd/(sd + w)
    tt2[tt2 == Inf] = NA
    sds <- rep(0, nbr - 1)
    for (i in 1:(nbr - 1)) {
      sds[i] <- mad(tt2[a == i], na.rm = TRUE)
    }
    cv.sd[j] <- sqrt(var(sds))/mean(sds)
  }
  o = (1:length(s0.perc))[cv.sd == min(cv.sd)]
  s0.hat = quantile(sd[sd != 0], s0.perc[o])
  return(list(s0.perc = s0.perc, cv.sd = cv.sd, s0.hat = s0.hat))
}

print('Building arrays')
print(Sys.time())

#build bundle arrays
for (file_idx in 1:length(groups[,"Identifier"])){
  #storage and file path initialization
  loc <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,"/",groups[file_idx,"Identifier"],"/cms")
  stuff <- array(NA, dim=c(num_feat,num_feat,num_clusts))
  zstuff <- array(NA, dim=c(num_feat,num_feat,num_clusts))
  
  # grab cells/node for the current patient
  if (groups[file_idx,"Tube"] == "A"){
    N <- countsA[,groups[file_idx,"Identifier"]]
  }else{
    N <- countsB[,groups[file_idx,"Identifier"]]
  }
  
  # ensure each node has at least four cells
  N[N<4] = 4
  
  
  for (n in 1:num_clusts) {
    # create array of correlation matrices
    bollocks <- paste("cm_", n, sep='')
    fil <- file.path(loc, paste(bollocks, ".Rda", sep=''))
    load(fil)
    stuff[,,n] <- get(bollocks)
    # fisher z transform correlation matrices
    ztemp <- get(bollocks)
    ztemp[row(ztemp)<col(ztemp)] <- fisherz(ztemp[row(ztemp)<col(ztemp)])
    ztemp[row(ztemp)>=col(ztemp)] <- NA
    #ztemp[row(ztemp)!=col(ztemp)] <- fisherz(ztemp[row(ztemp)!=col(ztemp)])
    #ztemp[row(ztemp) == col(ztemp)] <- 3 # as rho -> 1, z-> 3/inf
    zstuff[,,n] <- ztemp
    rm(list=bollocks)
  }
  #save results
  dimnames(stuff) <- list(marker_names,marker_names,1:num_clusts)
  assign(paste("group", file_idx, sep=""), list(stuff,N))
  dimnames(zstuff) <- list(marker_names,marker_names,1:num_clusts)
  assign(paste("zgroup", file_idx, sep=""), list(zstuff,N))
  rm(stuff, zstuff, ztemp, n, bollocks, fil) 
}

print('Computing residuals and deviations')
print(Sys.time())

all_rs <- matrix(NA,num_pairs*num_clusts,num_comparisons)
all_sigmas <- all_rs
count <-1

for (first_grp in 1:(length(groups[,"Identifier"])-1)){
  for (sec_grp in (first_grp+1):length(groups[,"Identifier"])){
    z1 <- get(paste0("zgroup",first_grp))[[1]]
    z2 <- get(paste0("zgroup",sec_grp))[[1]]
    N1 <- get(paste0("zgroup",first_grp))[[2]]
    N2 <- get(paste0("zgroup",sec_grp))[[2]]
    
    #find residuals
    rtemp <- z1 - z2
    r <- as.vector(rtemp)
    all_rs[,count] <- r[!is.na(r)]
    #inefficient way of storing the residuals in a vector by biomarker pair
    #r <- numeric()
    # for (idx1 in 1:(num_feat-1)){
    #   for (idx2 in (idx1+1):num_feat){
    #     r <- c(r,rtemp[idx1,idx2,])
    #   }
    # }
    #replicate the sds for each biomarker pair
    all_sigmas[,count] <- rep(sqrt(1/(N1-3)+1/(N2-3)),eacb=num_pairs)

    count <- count+1
    print(count)
    #concatenate with previous results
    # if (first_grp == 1){
    #   all_rs <- r
    #   all_sigmas <- sigma
    # }else{
    #   all_rs <- c(all_rs,r)
    #   all_sigmas <- c(all_sigmas,sigma)
    # }
  }
  # should remove zgroup first_grp here # OTHER QUESTION
}
rm(list = ls(pattern = "zgroup"))
rm(z1,z2,N1,N2)

all_s0s <- sapply(1:ncol(all_rs),function(i){
  s0_comp(all_rs[,i],all_sigmas[,i],s0.perc = seq(0,1,by = 0.05))
})

print(length(all_rs))
print(length(all_sigmas))

print('Computing s0')
print(Sys.time())

#des_output <- s0_comp(all_rs,all_sigmas,s0.perc = seq(0,1,by = 0.05))

print('Finished')
print(Sys.time())

