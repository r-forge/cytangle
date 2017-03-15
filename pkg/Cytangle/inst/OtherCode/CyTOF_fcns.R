## Function declarations
# Updated 3-5-17

be_positive <-function(x){
  # Prevents approximately zero, negative values from causing issues
  x[(x<0) & (x<1e-3)] <- -x[(x<0) & (x<1e-3)]
  # if ((x<0) && (x<1e-3)){  
  #   xplus <- -x
  # }else{
  #   xplus <- x
  # }
  return(x)
}

remove_hyphen <- function(x,base){
  temp <- x
  sharedInd <- grepl('-',temp) & !(temp %in% c('DNA-1','DNA-2','Ki-67','HLA-DR'))
  if (grepl('TA',base)){  #If it's Tube A
    temp[sharedInd] <- gsub('-[A-z0-9]*',"",temp[sharedInd]) 
  }else{ # If it's Tube B
    temp[sharedInd] <- gsub(".*-","",temp[sharedInd])
  }  
  return(temp)
}

symmetrize <- function(x){
  len <- dim(x)[1] # x should be square
  x+t(x)-diag(diag(x),len,len) 
}

sym_array <- function(x){
  l <- dim(x)[1]
  w <- dim(x)[2]
  h <- dim(x)[3]
  SymData <- apply(x,3,symmetrize)  
  array(SymData,c(l,w,h))
}

getsub <- function(x){
  badrow <- apply(x,1,function(x) all(is.na(x)))
  x[!badrow,!badrow]
}

find_inds = function(x){
  which(x == max(x), arr.ind = TRUE)
}

cor_est <- function(x,a_0,b_0,m_0,B_0){
  # computes estimated correlation matrix even for low number of samples
  # inspired by: http://www.fil.ion.ucl.ac.uk/~wpenny/publications/bmn.pdf
  # Created 5/25/16 by R. McGee 
  # Updated 6/8/16 by R. McGee 
  force(x)
  force(a_0)
  force(b_0)
  force(m_0)
  force(B_0)
  
  x <- as.matrix(x)
  N <- dim(x)[1] # number of samples
  d <- dim(x)[2] # dimension of ambient space
  x_bars = colMeans(x)
  factor <- x - rep(x_bars, rep.int(nrow(x), ncol(x))) 
  #factor <-scale(x,center=TRUE,scale=FALSE)
  #attr(factor,"scaled:center") <- NULL
  Sigma_N <- (1/N)*t(factor) %*% factor # factor is named, must fix that
  
  a_N <- a_0 + N/2
  b_N <- b_0 + N
  B_N <- B_0 + (N/2)*(Sigma_N + (b_0/b_N)*( (x_bars - m_0) %*% t(x_bars - m_0) ))
  
  tempC <- (1/a_N)*B_N
  sd_fac <- diag(1/sqrt(diag(tempC)))
  C_N <- sd_fac %*% tempC %*% sd_fac
  return(C_N)
}

correlate <- function(tempMat,num_feat){
  if(dim(tempMat)[1] == 0){
    eye26  
  }else{
    cor_est(tempMat,num_feat,1,0,eye26)
  }
}

process <- function(frame,clust_idx,indices){
  #when you're passing a whole frame
  curr_clustInd <- (frame[,"cluster"] == clust_idx)
  tempMat <- frame[curr_clustInd,indices] 
  y <- correlate(tempMat,length(indices))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

cm_plots <- function(inds,cm_path,out_path,marker_names,withbar){
  require(colorspace) 
  require(fields)
  famu3 <- colorRampPalette(c("#2AA416" , "gray90", "#FF7F00" ))(128)
  #image(matrix(1:64, nrow=1), col=crp)
  famu2 <- diverge_hcl(128, h = c(130, 43), c = 250, l = c(70, 90))
  
  #grab node pair and patient number where maximum occurred
  clust_idx1 <- inds[1,1]
  clust_idx2 <- inds[1,2]
  
  load(paste0(cm_path,"/",list.files(path = cm_path,pattern = paste0('cm_',clust_idx1,".Rda"))))
  load(paste0(cm_path,"/",list.files(path = cm_path,pattern = paste0('cm_',clust_idx2,".Rda"))))
  
  file1<-gsub(".Rda","",list.files(path = cm_path,pattern = paste0('cm_',clust_idx1,".Rda")))
  file2<-gsub(".Rda","",list.files(path = cm_path,pattern = paste0('cm_',clust_idx2,".Rda")))
  
  num_feat <- length(marker_names)
  
  opar <- par()
  if (withbar){
    dev.new(width=11,height=5)
    layout(matrix(c(1,2,3),nrow = 1),heights = c(1,1,1),widths = c(1, 5,5))
    opar <- par(mai=c(1, 0.3, 1, 0.1))
    image(1, seq(-1, 1, length=128), matrix(1:128, nrow=1), col=famu3, ylab='', xlab='', xaxt='n')
    par(opar)
    #par(mfrow=c(1,2)) 
    image(1:num_feat,1:num_feat,get(file1),
          col = famu3 ,zlim = c(-1,1), xaxt = "n", yaxt = "n",xlab="",ylab="")
    #image(-32:32, matrix(-32:32, ncol=1), col=mypalette)
    mtext(marker_names,side=c(2,2),at=1:num_feat,las=2)
    mtext(marker_names,side=c(1,1),at=1:num_feat,las=2)
    #dev.new(width=8,height=8)
    image(1:num_feat,1:num_feat,get(file2),
          col = famu3 ,zlim = c(-1,1), xaxt = "n", yaxt = "n",xlab="",ylab="")
  }else{
    dev.new(width=10,height=5)
    layout(matrix(c(1,2),nrow = 1),heights = c(1,1),widths = c(5,5))
    #opar <- par(mai=c(1, 0.3, 1, 0.1))
    #image(1, seq(-1, 1, length=128), matrix(1:128, nrow=1), col=famu3, ylab='', xlab='', xaxt='n')
    #par(opar)
    #par(mfrow=c(1,2)) 
    image(1:num_feat,1:num_feat,get(file1),
          col = famu3 ,zlim = c(-1,1), xaxt = "n", yaxt = "n",xlab="",ylab="")
    #image(-32:32, matrix(-32:32, ncol=1), col=mypalette)
    mtext(marker_names,side=c(2,2),at=1:num_feat,las=2)
    mtext(marker_names,side=c(1,1),at=1:num_feat,las=2)
    #dev.new(width=8,height=8)
    image(1:num_feat,1:num_feat,get(file2),
          col = famu3 ,zlim = c(-1,1), xaxt = "n", yaxt = "n",xlab="",ylab="")
  }
  mtext(marker_names,side=c(2,2),at=1:num_feat,las=2)
  mtext(marker_names,side=c(1,1),at=1:num_feat,las=2)
  title(paste0(file1," vs. ",file2))
 
  #dev.print(pdf, width=32,file = paste0(out_path,"/",file1,file2,".png"))
  
  x = get(file1) - get(file2)
  max_cor_marker_inds <- find_inds(x)[1,]
  #max_cor_feat_inds <- which(featInd %in% marker_names[max_cor_marker_inds])
  #max_cor_feature_change <- setdiff(tubeA_markers,surface_markers)[max_cor_feat_inds]
  # dev.new(width=8,height=8)
  # par(mfrow=c(1,1))
  # image(1:num_feat,1:num_feat,x,col = famu ,zlim = c(-1,1), xaxt = "n", yaxt = "n",xlab="",ylab="")
  # mtext(marker_names,side=c(2,2),at=1:num_feat,las=2)
  # mtext(marker_names,side=c(1,1),at=1:num_feat,las=2)
  return(marker_names[max_cor_marker_inds])
}

extract_max <- function(inds,cm_path,marker_names){
  clust_idx1 <- inds[1,1]
  clust_idx2 <- inds[1,2]
  
  load(paste0(cm_path,"/",list.files(path = cm_path,pattern = paste0('cm_',clust_idx1,".Rda"))))
  load(paste0(cm_path,"/",list.files(path = cm_path,pattern = paste0('cm_',clust_idx2,".Rda"))))
  
  file1<-gsub(".Rda","",list.files(path = cm_path,pattern = paste0('cm_',clust_idx1,".Rda")))
  file2<-gsub(".Rda","",list.files(path = cm_path,pattern = paste0('cm_',clust_idx2,".Rda")))
  
  x = get(file1) - get(file2)
  max_cor_marker_inds <- find_inds(x)[1,]
  
  # dev.new(width=8,height=8)
  # image(1:num_feat,1:num_feat,get(paste0('cm_',clust_idx1)),
  #       col = famu ,zlim = c(-1,1), xaxt = "n", yaxt = "n",xlab="",ylab="")
  # #image(-32:32, matrix(-32:32, ncol=1), col=mypalette)
  # mtext(marker_names,side=c(1,2),at=1:num_feat,las=2)
  # dev.new(width=8,height=8)
  # image(1:num_feat,1:num_feat,get(paste0('cm_',clust_idx2)),
  #       col = famu ,zlim = c(-1,1), xaxt = "n", yaxt = "n",xlab="",ylab="")
  # mtext(marker_names,side=c(1,2),at=1:num_feat,las=2)
  # 
  # x = get(paste0('cm_',clust_idx1)) - get(paste0('cm_',clust_idx2))
  max_cor_marker_inds <- find_inds(x)[1,]
  #max_cor_feat_inds <- which(featInd %in% marker_names[max_cor_marker_inds])
  max_cor_feature_change <- marker_names[max_cor_marker_inds] #setdiff(tubeA_markers,surface_markers)[max_cor_feat_inds]
  return(max_cor_feature_change)
}

overlay <- function(ucoords,lwd=2) {
  if (!require("plotrix")) { # to draw ellipses
    install.packages("plotrix")
    library(plotrix)
  }

  # this function draws the annotation overlay for kinds of cells
  opar <- par(lwd=lwd)
  on.exit( par(opar) )
  text(0.8, 0.035, "CD34+ CD38-")
  draw.ellipse(0.53, 0.17, 0.25, 0.09, -17, lwd=lwd)
  
  text(0.8, 0.48, "CD34+")
  polyvert <- ucoords[c(3, 291, 238, 443, 104, 329, 479, 449),] +
    0.02*matrix(c(-1, 0, 0, 1, 1,  0,  0,  -1,
                  1, 1, 1, 1, 0, -1, -1, -1), ncol=2)
  polygon(polyvert)
  
  text(-0.95, -0.9,  "* Promyelocytes", adj=0)
  text(-0.14, 0.45, "*")
  draw.ellipse(-0.095, 0.45, 0.03, 0.03, lwd=lwd)
  
  text(0.05, 0, "Platelets", adj=0.5)
  draw.ellipse(0.04, 0.14, 0.04, 0.13, 20, lwd=lwd)
  
  text(-0.53, 0.85, "Mature", adj=0.5)
  text(-0.53, 0.80, "granulocytes", adj=0.5)
  draw.ellipse(-0.18, 0.56, 0.045, 0.045, lwd=lwd)
  arrows(-0.53, 0.77, -0.23, 0.57, length=0.15, angle=20)
  
  text(-0.05, 0.95, "Early", adj=0.5)
  text(-0.05, 0.90, "granulocytes", adj=0.5)
  polyvert <- ucoords[c(197, 271, 221, 188, 48, 428, 170, 257, 314, 37),] +
    0.02*matrix(c(-1, -1, -1, -1, -1, 0, 1, 1, 1, 1,
                  -1, 0, 0, 0, 1, 1, 1, 0, 0, -1), ncol=2)
  polygon(polyvert)
  arrows(-0.05, 0.87, -0.05, 0.80, length=0.15, angle=20)
  
  text(0.55, -0.95, "Plasma cells", adj=0)
  draw.ellipse(0.4424, -0.95, 0.1, 0.1, lwd=lwd)
  
  text(0.55, -0.76, "Early B cells", adj=0)
  draw.ellipse(0.45, -0.76, 0.11, 0.085, 135, lwd=lwd)
  
  text(0.52, -0.56, "Mature B cells", adj=0)
  draw.ellipse(0.35, -0.58, 0.20, 0.11, 35, lwd=lwd)
  
  text(-0.72, -0.15, "Mature CD14+", adj=1)
  text(-0.72, -0.2, "Monocytes", adj=1)
  draw.ellipse(-0.4, -0.10, 0.32, 0.15, 7, lwd=lwd)
  
  text(-0.3, -0.45, "pDCs", adj=1)
  draw.ellipse(-0.21, -0.45, 0.12, 0.05, 32, lwd=lwd)
  
  text(-0.9, 0.77, "NK cells", adj=0.5)
  draw.ellipse(-0.84, 0.59, 0.22, 0.14, 160, lwd=lwd)
  
  text(-0.7, 0.2, "T cells", adj=0.5)
  polyvert <- ucoords[c(11, 418, 355, 409, 460, 340, 234, 78, 453),] +
    0.03*matrix(c(1, 1, 0, 0, -1, -1, -1, -1, 1,
                  1, -1, -1, -1, -1, 1, 0, 1, 1), ncol=2)
  polygon(polyvert)
  
  text(-0.35,-0.8, "Early monocytes", adj=0.5)
  polyvert <- ucoords[c(311, 349, 413, 184, 445, 266, 6, 416, 303, 96,
                        165, 478, 467, 356),] +
    0.02*matrix(c(-2, -1, -1, -1, -1, -1, -1, 1, 2, 2,  1,  1,  1, 0,
                  0,  1,  1,  0,  0,  0,  1, 1.5, 0, 0, -1, -1, -1, -1), ncol=2)
  polygon(polyvert)
  
  text(0.1, 0.85, "Erythroid cells", adj=0)
  polyvert <- ucoords[c(258, 415, 92, 403, 395, 53, 286, 354, 204),] +
    0.02*matrix(c( 0, -1, -1, -1, -1, -1, 1, 1, 1,
                   -1,  0,  0,  0,  0,  1, 1, 0, -1), ncol=2)
  polygon(polyvert)
  
  text(0.2, 0.06, "Basophils")
  polyvert <- ucoords[c(57, 334, 160, 160, 335, 339, 281, 307, 112),] +
    0.02*matrix(c(-1, -1, -1,  2, -1, 1, 1, 1, 1,
                  1,  0, -1, -1,  0, 0, 0, 0, 1), ncol=2)
  polygon(polyvert)
}

mst_cor0 <- function(ccs,I,J,num_clusts, crp=famu2,title) {
  idxs <- 1:num_clusts #setdiff(1:num_clusts,c(I,J))
  vals <- ccs#[idxs]
  delta <- diff(range(vals))/128 
  foo <- seq(min(vals)-delta, max(vals)+delta, length=128) 
  disco <- cut(vals, breaks=foo, labels=FALSE)
  shades <- crp[disco]
  mst3 <- set_vertex_attr(mst2, "size", value=4)
  mst4 <- set_vertex_attr(mst3, "color", value=shades, index = idxs)
  #mst5 <- set_vertex_attr(mst4, "color", value="purple", index=I)
  mst5 <- mst4
  mst6 <- set_vertex_attr(mst5,"label",value="")
  mst7 <- mst6 #set_vertex_attr(mst6, "color", value="blue", index=J)
  plot(mst7, main=title, layout=lay)
  invisible(mst7)
}


mst_cor <- function(folder,ccs,I,J,num_clusts,crp,title) {
  require(igraph)
  # main spade plotting
  mst <- read_graph(paste0(folder,"/mst.gml"), format="gml")
  lay <- as.matrix(read.table(paste0(folder,"/layout.table")))
  idxs <- setdiff(1:num_clusts,c(I,J))
  vals <- ccs[idxs]
  delta <- diff(range(vals))/128 
  foo <- seq(min(vals)-delta, max(vals)+delta, length=128) 
  disco <- cut(vals, breaks=foo, labels=FALSE)
  shades <- crp[disco]
  mst3 <- set_vertex_attr(mst, "size", value=4)
  mst4 <- set_vertex_attr(mst3, "color", value=shades, index = idxs)
  #mst5 <- set_vertex_attr(mst4, "color", value="purple", index=I)
  #mst5 <- mst4
  mst6 <- set_vertex_attr(mst4,"label",value="")
  #mst7 <- set_vertex_attr(mst6, "color", value="blue", index=J)
  plot(mst6, main=title, layout=lay)
  # needed if using overlay
  xform <- function(v) {
    X <- max(v)
    N <- min(v)
    a <- 2/(X-N)
    b <- (N+X)/(N-X)
    function(y) b + a*y
  }
  f1 <- xform(lay[,1])
  u1 <- f1(lay[,1])
  f2 <- xform(lay[,2])
  u2 <- f2(lay[,2])
  ucoords <- cbind(u1, u2)
  overlay(ucoords = ucoords) #make this an if check
  invisible(mst6)
}

chi_delta_patients <- function(group1,group2=NULL,compar_size=NULL,tube,s0_val=0.03268146){
  
  marker_names <- get(paste0("tube",tube,"_panel"))[featInd-2]
  
  group_list <- c('CBF-AML',
                  'APL',
                  'NK-AML FLT3-ITD',
                  'NK-AML FLT3wt',
                  'Adverse-risk',
                  'Normal Karotype AML FLT3 tyrosine kinase domain mutation',
                  'Granulylocytic sarcoma',
                  'MLL rearrangement',
                  'Complete Recovery',
                  'Normal')
  
  markerInd <- featInd 
  num_feat <- length(markerInd)
  num_pairs <- choose(num_feat,2)
  num_grps <- length(group_list)
  num_patients <- length(groups[,"Identifier"])/2
  num_subtypes <- table(groups[1:num_patients,"Group Number"])
  
  #mesh grid
  # npairs <- data.frame(rep(n1,each = length(n2)),rep(n2,times = length(n1)))
  if (is.null(group2)){
    if (is.null(compar_size)){
      subgrp1 <- unique(groups[groups[,"Group Number"] == group1,"Identifier"])
      subgrp2 <- subgrp1
      combos <- choose(num_subtypes[[toString(group1)]],2)
    }else{
      subgrp1 <- sample(unique(groups[groups[,"Group Number"] == group1,"Identifier"]),compar_size)
      subgrp2 <- setdiff(unique(groups[groups[,"Group Number"] == group1,"Identifier"]),subgrp1)
      combos <- length(subgrp1)*length(subgrp2)
    }
  }else{
    subgrp1 <- unique(groups[groups[,"Group Number"] == group1,"Identifier"])
    subgrp2 <- unique(groups[groups[,"Group Number"] == group2,"Identifier"])
    combos <- length(subgrp1)*length(subgrp2)
  }
  cur_grps <- c(subgrp1,subgrp2)
  
  #dimnames of output vector
  chi_deltas <- matrix(0,combos,choose(num_feat,2))
  temp1 <- unlist(sapply(17:1, function(i) { rep(marker_names[18-i],each=i)}))
  temp2 <- unlist(sapply(2:18, function(i) { rep(marker_names[i:num_feat],each=1)}))
  #rep(marker_names[2:num_feat],times=choose(18,2)/(num_feat-1))
  bio_pairs <- paste0(temp1,";",temp2)
  rm(temp1,temp2)
  if (is.null(group2)&(is.null(compar_size))){
    grp_pairs <- combn(subgrp1,2, FUN = paste, collapse = ";")
  }else{
    grp_pairs <- apply(expand.grid(subgrp1, subgrp2), 1, paste, collapse=";") #paste0(temp1,";",temp2)
  }
  dimnames(chi_deltas) <- list(grp_pairs,bio_pairs)
  
  #store exchangeability factor
  S0 <- rep(s0_val,num_clusts)
  
  
  #fisher z transform
  fisherz <- function(x) {
    0.5*(log(1+x) - log(1-x))
  }
  
  #count which patient in each group we're on
  count <- rep(0,10)
  
  #build bundle arrays
  for (file_idx in 1:length(groups[,"Identifier"])){
    if (groups[file_idx,"Identifier"] %in% cur_grps){
      if (groups[file_idx,"Tube"] == tube){
        count[as.numeric(groups[file_idx,"Group Number"])] <- count[as.numeric(groups[file_idx,"Group Number"])] + 1
        loc <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,"/",groups[file_idx,"Identifier"],"/cms")
        stuff <- array(NA, dim=c(num_feat,num_feat,num_clusts))
        zstuff <- array(NA, dim=c(num_feat,num_feat,num_clusts))
        
        #store number of cells per node. make sure it's at least 4
        N <- get(paste0("counts",tube))[,groups[file_idx,"Identifier"]]
        N[N<4] = 4
        
        for (n in 1:num_clusts) {
          # create array of correlation coefficients
          bollocks <- paste("cm_", n, sep='')
          fil <- file.path(loc, paste(bollocks, ".Rda", sep=''))
          load(fil)
          stuff[,,n] <- get(bollocks)
          # fisher z transform correlation coefficients
          ztemp <- get(bollocks)
          ztemp[row(ztemp)>col(ztemp)] <- fisherz(ztemp[row(ztemp)>col(ztemp)])
          ztemp[row(ztemp) <= col(ztemp)] <- NA # as rho -> 1, z-> 3/inf
          zstuff[,,n] <- ztemp
          rm(list=bollocks)
        }
        dimnames(stuff) <- list(marker_names,marker_names,1:num_clusts)
        assign(paste("group_",groups[file_idx,"Identifier"], sep=""), list(stuff,N))
        dimnames(zstuff) <- list(marker_names,marker_names,1:num_clusts)
        assign(paste("zgroup_", groups[file_idx,"Identifier"], sep=""), list(zstuff,N))
        rm(stuff, zstuff, ztemp, n, bollocks, fil) 
      }
      
    }
  }
  
  if (is.null(group2)&(is.null(compar_size))){
    for (first_idx in 1:(length(subgrp1)-1)){
      for (sec_idx in (first_idx+1):length(subgrp1)){
        
        first_grp = subgrp1[first_idx]
        sec_grp = subgrp1[sec_idx]
        
        z1 <- get(paste0("zgroup_",first_grp))[[1]]
        z2 <- get(paste0("zgroup_",sec_grp))[[1]]
        N1 <- get(paste0("zgroup_",first_grp))[[2]]
        N2 <- get(paste0("zgroup_",sec_grp))[[2]]
        
        cur_row <- paste0(first_grp,";",sec_grp)
        
        # compute chi squared values
        temp <- sweep((z1 - z2)^2,MARGIN = 3,(sqrt(1/(N1-3)+1/(N2-3))+S0)^2,"/")
        chi_sq <- apply(temp, 1:2, sum)
        chi_deltas[cur_row,]<-chi_sq[lower.tri(chi_sq)]
      }
    }
  }else{
    for (first_grp in subgrp1){
      for (sec_grp in subgrp2){
        z1 <- get(paste0("zgroup_",first_grp))[[1]]
        z2 <- get(paste0("zgroup_",sec_grp))[[1]]
        N1 <- get(paste0("zgroup_",first_grp))[[2]]
        N2 <- get(paste0("zgroup_",sec_grp))[[2]]
        
        cur_row <- paste0(first_grp,";",sec_grp)
        
        # compute chi squared values
        temp <- sweep((z1 - z2)^2,MARGIN = 3,(sqrt(1/(N1-3)+1/(N2-3))+S0)^2,"/")
        chi_sq <- apply(temp, 1:2, sum)
        chi_deltas[cur_row,]<-chi_sq[lower.tri(chi_sq)]
      }
    }
  }
  return(chi_deltas)
}

chi_delta_groups_vs_normals <- function(tube,s0_val=0.03268146){
  
  marker_names <- get(paste0("tube",tube,"_panel"))[featInd-2]
  
  group_list <- c('CBF-AML',
                  'APL',
                  'NK-AML FLT3-ITD',
                  'NK-AML FLT3wt',
                  'Adverse-risk',
                  'Normal Karotype AML FLT3 tyrosine kinase domain mutation',
                  'Granulylocytic sarcoma',
                  'MLL rearrangement',
                  'Complete Recovery',
                  'Normal')
  
  markerInd <- featInd 
  num_feat <- length(markerInd)
  num_pairs <- choose(num_feat,2)
  num_grps <- length(group_list)
  num_patients <- length(groups[,"Identifier"])/2
  num_subtypes <- table(groups[1:num_patients,"Group Number"])
  
  #dimnames of output vector
  chi_deltas <- matrix(0,num_grps-1,choose(num_feat,2))
  temp1 <- unlist(sapply(17:1, function(i) { rep(marker_names[18-i],each=i)}))
  temp2 <- unlist(sapply(2:18, function(i) { rep(marker_names[i:num_feat],each=1)}))
 #rep(marker_names[2:num_feat],times=choose(18,2)/(num_feat-1))
  bio_pairs <- paste0(temp1,";",temp2)
  rm(temp1,temp2)
  grp_pairs <- paste0(c('CBF-AML', 
                        'APL',
                        'FLT3-ITD',
                        'FLT3wt',
                        'AdverseRisk',
                        'FLT3-TKD',
                        'Granusarcoma',
                        'MLL',
                        'CompleteRecovery'),';Normal')
  dimnames(chi_deltas) <- list(grp_pairs,bio_pairs)
  
  #store exchangeability factor
  S0 <- rep(s0_val,num_clusts)
  
  
  #fisher z transform
  fisherz <- function(x) {
    0.5*(log(1+x) - log(1-x))
  }
  
  #count which patient in each group we're on
  count <- rep(0,10)
  
  #build bundle arrays
  for (g in 1:num_grps) {
    loc <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,"/grp_cms/group",g)
    stuff <- array(NA, dim=c(num_feat,num_feat,num_clusts))
    zstuff <- array(NA, dim=c(num_feat,num_feat,num_clusts))
    
    #store number of cells per node. make sure it's at least 4
    N <- get(paste0("counts",tube))[,paste0("Group ",g)]
    N[N<4] = 4
    
    for (n in 1:num_clusts) {
      # create array of correlation coefficients
      bollocks <- paste("grp", g, "_cm_", n, sep='')
      fil <- file.path(loc, paste(bollocks, ".Rda", sep=''))
      load(fil)
      stuff[,,n] <- get(bollocks)
      # fisher z transform correlation coefficients
      ztemp <- get(bollocks)
      ztemp[row(ztemp)>col(ztemp)] <- fisherz(ztemp[row(ztemp)>col(ztemp)])
      ztemp[row(ztemp) <= col(ztemp)] <- NA # as rho -> 1, z-> 3/inf
      zstuff[,,n] <- ztemp
      rm(list=bollocks)
    }
    dimnames(stuff) <- list(marker_names,marker_names,1:num_clusts)
    assign(paste("group_", g, sep=""), list(stuff,N))
    dimnames(zstuff) <- list(marker_names,marker_names,1:num_clusts)
    assign(paste("zgroup_", g, sep=""), list(zstuff,N))
    rm(stuff, zstuff, ztemp, n, bollocks, fil) 
  }
  
  for (grp_idx in 1:(num_grps-1)){
    
    z1 <- get(paste0("zgroup_",grp_idx))[[1]]
    z2 <- zgroup_10[[1]]
    N1 <- get(paste0("zgroup_",grp_idx))[[2]]
    N2 <- zgroup_10[[2]]
    
    cur_row <- grp_pairs[grp_idx]
    
    # compute chi squared values
    temp <- sweep((z1 - z2)^2,MARGIN = 3,(sqrt(1/(N1-3)+1/(N2-3))+S0)^2,"/")
    chi_sq <- apply(temp, 1:2, sum)
    chi_deltas[cur_row,]<-chi_sq[lower.tri(chi_sq)]
  }
  return(chi_deltas)
}
