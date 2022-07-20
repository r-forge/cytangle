## Function declarations
# Updated 7-27-17

gm <- function(x) {
  exp(mean(log(x)))
}

hm <- function(x) {
  1/mean(1/x)
}

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

fisherz <- function(x) {
  0.5*(log(1+x) - log(1-x))
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

mst_cor0 <- function(ccs,num_clusts,crp,title,scale = 'relative') {
  force(folder)
  # main spade plotting
  mst <- read_graph(paste0(folder,"/mst.gml"), format="gml")
  lay <- as.matrix(read.table(paste0(folder,"/layout.table")))
  idxs <- 1:num_clusts
  vals <- ccs
  if (scale %in% 'absolute'){
    delta <- 1/128
    foo <- seq(-1-delta, 1+delta, length=128) 
  }else{
    delta <- diff(range(vals))/128 
    foo <- seq(min(vals)-delta, max(vals)+delta, length=128) 
  }
  disco <- cut(vals, breaks=foo, labels=FALSE)
  shades <- crp[disco]
  mst3 <- set_vertex_attr(mst, "size", value=4)
  mst4 <- set_vertex_attr(mst3, "color", value=shades, index = idxs)
  #mst5 <- set_vertex_attr(mst4, "color", value="purple", index=I)
  #mst5 <- mst4
  mst6 <- set_vertex_attr(mst4,"label",value="")
  #mst7 <- set_vertex_attr(mst6, "color", value="blue", index=J)
  plot(mst6, main=title, layout=lay)
  invisible(mst7)
}


mst_cor <- function(folder,ccs,num_clusts,crp,title,scale = 'relative') {
  force(folder)
  # main spade plotting
  mst <- read_graph(paste0(folder,"/mst.gml"), format="gml")
  lay <- as.matrix(read.table(paste0(folder,"/layout.table")))
  idxs <- 1:num_clusts
  vals <- ccs
  if (scale %in% 'absolute'){
    delta <- 1/128
    foo <- seq(-1-delta, 1+delta, length=128) 
  }else{
    delta <- diff(range(vals))/128 
    foo <- seq(min(vals)-delta, max(vals)+delta, length=128) 
  }
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
  overlay() #make this an if check
  invisible(mst6)
}

mst_compare <- function(folder,ccs,foo,num_clusts,crp,title,which = 2) {
  force(folder)
  # main spade plotting
  mst <- read_graph(paste0(folder,"/mst.gml"), format="gml")
  lay <- as.matrix(read.table(paste0(folder,"/layout.table")))
  idxs <- 1:num_clusts
  vals <- ccs
  disco <- cut(vals, breaks=foo, labels=FALSE)
  shades <- crp[disco]
  mst3 <- set_vertex_attr(mst, "size", value=4)
  mst4 <- set_vertex_attr(mst3, "color", value=shades, index = idxs)
  #mst5 <- set_vertex_attr(mst4, "color", value="purple", index=I)
  #mst5 <- mst4
  mst6 <- set_vertex_attr(mst4,"label",value="")
  #mst7 <- set_vertex_attr(mst6, "color", value="blue", index=J)
  plot(mst6, main=title, layout=lay)
  if (which == 1){
    fields::image.plot( zlim=c(min(foo), max(foo)), col = famu3,horizontal = TRUE,legend.only=TRUE,legend.shrink = 0.25,xlab = 'Correlation scale')
    #image(1, seq(min(foo), max(foo), length=128), matrix(1:128, nrow=1), col=famu3, ylab='', xlab='', xaxt='n')
    # scalebar(foo, xy = NULL, type = "bar", divs = 2, below = 'Correlation scale',
  #           lonlat = NULL, label = c(toString(foo[1]),toString((foo[1]+foo[length(foo)])/2),toString(foo[length(foo)])), adj=c(0.5, -0.5), lwd = 2)
  }
  overlay() #make this an if check
  invisible(mst6)
}


min_span_cor <- function(folder,both_ccs,num_clusts,crp,titles,scale = 'relative'){
  par(mfrow = c(1,2),mar=c(1,1,3,1)+0.1) # mar: c(bottom, left, top, right)
  if (scale %in% 'absolute'){
    delta <- 1/128
    foo <- seq(-1-delta, 1+delta, length=128) 
  }else{
    delta <- diff(range(both_ccs))/128 
    foo <- seq(min(both_ccs)-delta, max(both_ccs)+delta, length=128) 
  }
 # layout(matrix(c(1,2,3),nrow = 1),heights = c(1,1,1),widths = c(1, 5,5))
 # opar <- par(mai=c(1, 0.3, 1, 0.1))
#  image(1, seq(-1, 1, length=128), matrix(1:128, nrow=1), col=famu3, ylab='', xlab='', xaxt='n')
  
  mst_compare(folder,both_ccs[1,],foo,num_clusts,crp,titles[1],which = 1)
  mst_compare(folder,both_ccs[2,],foo,num_clusts,crp,titles[2],which = 2) 
  #image(1, seq(min(foo), max(foo), length=128), matrix(1:128, nrow=1), col=famu3, ylab='', xlab='Correlation scale', xaxt='n')
}


prev_mst_cor0 <- function(ccs,I,J,num_clusts, crp=famu2,title) {
  idxs <- 1:num_clusts #setdiff(1:num_clusts,c(I,J))
  vals <- ccs#[idxs]
  delta <- 1/128 #diff(range(vals))/128 
  foo <- seq(-1-delta, 1+delta, length=128) 
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


prev_mst_cor <- function(folder,ccs,I,J,num_clusts,crp,title) {
  force(folder)
  # main spade plotting
  mst <- read_graph(paste0(folder,"/mst.gml"), format="gml")
  lay <- as.matrix(read.table(paste0(folder,"/layout.table")))
  idxs <- setdiff(1:num_clusts,c(I,J))
  vals <- ccs[idxs]
  delta <- 1/128 #diff(range(vals))/128 
  foo <- seq(-1-delta, 1+delta, length=128) 
  disco <- cut(vals, breaks=foo, labels=FALSE)
  shades <- crp[disco]
  mst3 <- set_vertex_attr(mst, "size", value=4)
  mst4 <- set_vertex_attr(mst3, "color", value=shades, index = idxs)
  #mst5 <- set_vertex_attr(mst4, "color", value="purple", index=I)
  #mst5 <- mst4
  mst6 <- set_vertex_attr(mst4,"label",value="")
  #mst7 <- set_vertex_attr(mst6, "color", value="blue", index=J)
  plot(mst6, main=title, layout=lay)
  overlay() #make this an if check
  invisible(mst6)
}

node_finder <- function(I,num_clusts,title) {
  folder <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/spade",run_date)
  # main spade plotting
  mst <- read_graph(paste0(folder,"/mst.gml"), format="gml")
  lay <- as.matrix(read.table(paste0(folder,"/layout.table")))
  idxs <- setdiff(1:num_clusts,c(I))
  mst3 <- set_vertex_attr(mst, "size", value=4)
  mst4 <- set_vertex_attr(mst3, "color", value='blue', index = idxs)
  mst6 <- set_vertex_attr(mst4,"label",value="")
  mst7 <- set_vertex_attr(mst6, "color", value="red", index=I)
  plot(mst7, main=title, layout=lay)
  overlay() #make this an if check
  invisible(mst7)
}

chi_delta_patients <- function(group1,group2=NULL,compar_size=NULL,tube,which,s0_val=0){
  
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
        loc <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,which,"/",groups[file_idx,"Identifier"],"/cms")
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

chi_delta_pooled_subsample <- function(group1,group2=NULL,compar_size=NULL,tube,which,s0_val=0){
  
  marker_names <- get(paste0("tube",tube,"_panel"))[featInd-2]
  
  group_list <- c('CBF-AML',
                  'APL',
                  'NK-AML FLT3-ITD',
                  'NK-AML FLT3wt',
                  'Adverse-risk',
                  'NK-AML FLT3-TKD',
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
  if (is.null(group2)){
    if (is.null(compar_size)){
      stop("Need a comparison size for now")
      #subgrp1 <- unique(groups[groups[,"Group Number"] == group1,"Identifier"])
      #subgrp2 <- subgrp1
      #combos <- choose(num_subtypes[[toString(group1)]],2)
    }else{
      subgrp1 <- sample(unique(groups[groups[,"Group Number"] == group1,"Identifier"]),compar_size)
      subgrp2 <- setdiff(unique(groups[groups[,"Group Number"] == group1,"Identifier"]),subgrp1)
      combos <- 1 #length(subgrp1)*length(subgrp2)
    }
  }else{
    stop("Only one group for now")
  }
  #   subgrp1 <- unique(groups[groups[,"Group Number"] == group1,"Identifier"])
  #   subgrp2 <- unique(groups[groups[,"Group Number"] == group2,"Identifier"])
  #   combos <- length(subgrp1)*length(subgrp2)
  # }
  cur_grps <- c(subgrp1,subgrp2)
  
  #dimnames of output vector
  chi_deltas <- matrix(0,combos,choose(num_feat,2))
  temp1 <- unlist(sapply(17:1, function(i) { rep(marker_names[18-i],each=i)}))
  temp2 <- unlist(sapply(2:18, function(i) { rep(marker_names[i:num_feat],each=1)}))
  #rep(marker_names[2:num_feat],times=choose(18,2)/(num_feat-1))
  bio_pairs <- paste0(temp1,";",temp2)
  rm(temp1,temp2)
  # if (is.null(compar_size)){ #if (is.null(group2)&(is.null(compar_size))){
  #   grp_pairs <- combn(subgrp1,2, FUN = paste, collapse = ";")
  # }else{
  #   grp_pairs <- apply(expand.grid(subgrp1, subgrp2), 1, paste, collapse=";") #paste0(temp1,";",temp2)
  # }
  grp_pairs <- paste0(group_list[group1],";",group_list[group1])
  dimnames(chi_deltas) <- list(grp_pairs,bio_pairs)
  
  #store exchangeability factor
  S0 <- rep(s0_val,num_clusts)
  
  
  #fisher z transform
  fisherz <- function(x) {
    0.5*(log(1+x) - log(1-x))
  }
  
  #count which patient in each group we're on
  count <- rep(0,10)
  
  # find pooled correlation matrices at every node
  for (n in 1:num_clusts) {
    for (file_idx in 1:length(groups[,"Identifier"])){
      if (groups[file_idx,"Identifier"] %in% cur_grps){
        if (groups[file_idx,"Tube"] == tube){
          count[as.numeric(groups[file_idx,"Group Number"])] <- count[as.numeric(groups[file_idx,"Group Number"])] + 1
          loc <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,which,"/",groups[file_idx,"Identifier"])
          
          # create pooled correlation matrices
          cur_file <- paste0(loc,"/node",n,".Rda")
          load(cur_file) #tempMat (patient "file_idx" info at node "clust_idx")
          if (!is.matrix(tempMat)){ #for size 1 x num_feat frames where dim() won't work
            cur_frame <- t(as.matrix(tempMat))[,markerInd,drop=FALSE] # drop forces this to stay a matrix
          }else{
            cur_frame <- tempMat[,markerInd]
          }
          rm(tempMat)
        
          # build group by group data frame (this builds a frame with all of the data for each group at a given node)
          if(groups[file_idx,"Identifier"] %in% subgrp1){
            if (!exists("pool_temp1")){
              pool_temp1 <- cur_frame 
            }else{
              pool_temp1 <- rbind(pool_temp1,cur_frame)
            }
          }else if(groups[file_idx,"Identifier"] %in% subgrp2){
            if (!exists("pool_temp2")){
              pool_temp2 <- cur_frame
            }else{
              pool_temp2 <- rbind(pool_temp2,cur_frame)
            }
          }
        }
      }
    }
    for (subgrp_idx in c(1,2)){
      bollocks <- paste0("pooled",subgrp_idx,"_cm_", n, sep='')
      cur_pool <- get(paste0("pool_temp",subgrp_idx))
      assign(bollocks,correlate(cur_pool,num_feat))
    }
    rm(bollocks,pool_temp1,pool_temp2)
  }
  for (subgrp_idx in c(1,2)){
    #build bundle arrays
    stuff <- array(NA, dim=c(num_feat,num_feat,num_clusts))
    zstuff <- array(NA, dim=c(num_feat,num_feat,num_clusts))
    cur_subgrp <- get(paste0("subgrp",subgrp_idx))
    Ntemp <- get(paste0("counts",tube))[,cur_subgrp]
    N <- as.matrix(apply(Ntemp,1,hm))
    N[N<4] = 4
    #N <- rowSums(get(paste0("counts",tube))[,cur_subgrp])
    #N[N<4] = 4
    for (n in 1:num_clusts) {
      #store number of cells per node. make sure it's at least 4
      bollocks <- paste0("pooled",subgrp_idx,"_cm_", n, sep='')
      stuff[,,n] <- get(bollocks)
      # fisher z transform correlation coefficients
      ztemp <- get(bollocks)
      ztemp[row(ztemp)>col(ztemp)] <- fisherz(ztemp[row(ztemp)>col(ztemp)])
      ztemp[row(ztemp) <= col(ztemp)] <- NA # as rho -> 1, z-> 3/inf
      zstuff[,,n] <- ztemp
    }
    dimnames(stuff) <- list(marker_names,marker_names,1:num_clusts)
    assign(paste("group_",subgrp_idx, sep=""), list(stuff,N))
    dimnames(zstuff) <- list(marker_names,marker_names,1:num_clusts)
    assign(paste("zgroup_",subgrp_idx, sep=""), list(zstuff,N))
    rm(stuff, zstuff, ztemp, bollocks) 
  }
  rm(list = ls(pattern = 'pooled*'))
  
  first_grp = 1
  sec_grp = 2
  
  z1 <- get(paste0("zgroup_",first_grp))[[1]]
  z2 <- get(paste0("zgroup_",sec_grp))[[1]]
  N1 <- get(paste0("zgroup_",first_grp))[[2]]
  N2 <- get(paste0("zgroup_",sec_grp))[[2]]
  
  cur_row <- paste0(group_list[group1],";",group_list[group1])
  
  # compute chi squared values
  temp <- sweep((z1 - z2)^2,MARGIN = 3,(sqrt(1/(N1-3)+1/(N2-3))+S0)^2,"/")
  chi_sq <- apply(temp, 1:2, sum)
  chi_deltas[cur_row,]<-chi_sq[lower.tri(chi_sq)]
  
  return(chi_deltas)
}

chi_delta_groups_vs_normals <- function(tube,which,s0_val=0.03268146){
  
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
    loc <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,which,"/grp_cms/group",g)
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

chi_delta_pooled_subsample2 <- function(group1,group2=NULL,compar_size=NULL,equal_compar=FALSE,tube,s0_val=0.03268146){
  
  marker_names <- get(paste0("tube",tube,"_panel"))[featInd-2]
  
  group_list <- c('CBF-AML',
                  'APL',
                  'NK-AML FLT3-ITD',
                  'NK-AML FLT3wt',
                  'Adverse-risk',
                  'NK-AML FLT3-TKD',
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
  if (is.null(group2)){
    if (is.null(compar_size)){
      stop("Need a comparison size for now")
      #subgrp1 <- unique(groups[groups[,"Group Number"] == group1,"Identifier"])
      #subgrp2 <- subgrp1
      #combos <- choose(num_subtypes[[toString(group1)]],2)
    }else{
      if ((equal_compar = TRUE)&&(2*compar_size <= num_subtypes[[toString(group1)]])){
        subgrp1 <- sample(unique(groups[groups[,"Group Number"] == group1,"Identifier"]),compar_size)
        subgrp2 <- sample(setdiff(unique(groups[groups[,"Group Number"] == group1,"Identifier"]),subgrp1),compar_size)
        combos <- 1 #length(subgrp1)*length(subgrp2)
      }else{
        subgrp1 <- sample(unique(groups[groups[,"Group Number"] == group1,"Identifier"]),compar_size)
        subgrp2 <- setdiff(unique(groups[groups[,"Group Number"] == group1,"Identifier"]),subgrp1)
        combos <- 1 #length(subgrp1)*length(subgrp2)
      }
    }
  }else{
    stop("Only one group for now")
  }
  #   subgrp1 <- unique(groups[groups[,"Group Number"] == group1,"Identifier"])
  #   subgrp2 <- unique(groups[groups[,"Group Number"] == group2,"Identifier"])
  #   combos <- length(subgrp1)*length(subgrp2)
  # }
  cur_grps <- c(subgrp1,subgrp2)
  
  #dimnames of output vector
  chi_deltas <- matrix(0,combos,choose(num_feat,2))
  temp1 <- unlist(sapply(17:1, function(i) { rep(marker_names[18-i],each=i)}))
  temp2 <- unlist(sapply(2:18, function(i) { rep(marker_names[i:num_feat],each=1)}))
  #rep(marker_names[2:num_feat],times=choose(18,2)/(num_feat-1))
  bio_pairs <- paste0(temp1,";",temp2)
  rm(temp1,temp2)
  # if (is.null(compar_size)){ #if (is.null(group2)&(is.null(compar_size))){
  #   grp_pairs <- combn(subgrp1,2, FUN = paste, collapse = ";")
  # }else{
  #   grp_pairs <- apply(expand.grid(subgrp1, subgrp2), 1, paste, collapse=";") #paste0(temp1,";",temp2)
  # }
  grp_pairs <- paste0(group_list[group1],";",group_list[group1])
  dimnames(chi_deltas) <- list(grp_pairs,bio_pairs)
  
  #store exchangeability factor
  S0 <- rep(s0_val,num_clusts)
  
  
  #fisher z transform
  fisherz <- function(x) {
    0.5*(log(1+x) - log(1-x))
  }
  
  #count which patient in each group we're on
  count <- rep(0,10)
  
  # find pooled correlation matrices at every node
  for (n in 1:num_clusts) {
    for (file_idx in 1:length(groups[,"Identifier"])){
      if (groups[file_idx,"Identifier"] %in% cur_grps){
        if (groups[file_idx,"Tube"] == tube){
          count[as.numeric(groups[file_idx,"Group Number"])] <- count[as.numeric(groups[file_idx,"Group Number"])] + 1
          loc <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,"/",groups[file_idx,"Identifier"])
          
          # create pooled correlation matrices
          cur_file <- paste0(loc,"/node",n,".Rda")
          load(cur_file) #tempMat (patient "file_idx" info at node "clust_idx")
          if (!is.matrix(tempMat)){ #for size 1 x num_feat frames where dim() won't work
            cur_frame <- t(as.matrix(tempMat))[,markerInd,drop=FALSE] # drop forces this to stay a matrix
          }else{
            cur_frame <- tempMat[,markerInd]
          }
          rm(tempMat)
          
          # build group by group data frame (this builds a frame with all of the data for each group at a given node)
          if(groups[file_idx,"Identifier"] %in% subgrp1){
            if (!exists("pool_temp1")){
              pool_temp1 <- cur_frame 
            }else{
              pool_temp1 <- rbind(pool_temp1,cur_frame)
            }
          }else if(groups[file_idx,"Identifier"] %in% subgrp2){
            if (!exists("pool_temp2")){
              pool_temp2 <- cur_frame
            }else{
              pool_temp2 <- rbind(pool_temp2,cur_frame)
            }
          }
        }
      }
    }
    for (subgrp_idx in c(1,2)){
      bollocks <- paste0("pooled",subgrp_idx,"_cm_", n, sep='')
      cur_pool <- get(paste0("pool_temp",subgrp_idx))
      assign(bollocks,correlate(cur_pool,num_feat))
    }
    rm(bollocks,pool_temp1,pool_temp2)
  }
  for (subgrp_idx in c(1,2)){
    #build bundle arrays
    stuff <- array(NA, dim=c(num_feat,num_feat,num_clusts))
    zstuff <- array(NA, dim=c(num_feat,num_feat,num_clusts))
    cur_subgrp <- get(paste0("subgrp",subgrp_idx))
    Ntemp <- get(paste0("counts",tube))[,cur_subgrp]
    N <- as.matrix(apply(Ntemp,1,hm))
    N[N<4] = 4
    #N <- rowSums(get(paste0("counts",tube))[,cur_subgrp])
    #N[N<4] = 4
    for (n in 1:num_clusts) {
      #store number of cells per node. make sure it's at least 4
      bollocks <- paste0("pooled",subgrp_idx,"_cm_", n, sep='')
      stuff[,,n] <- get(bollocks)
      # fisher z transform correlation coefficients
      ztemp <- get(bollocks)
      ztemp[row(ztemp)>col(ztemp)] <- fisherz(ztemp[row(ztemp)>col(ztemp)])
      ztemp[row(ztemp) <= col(ztemp)] <- NA # as rho -> 1, z-> 3/inf
      zstuff[,,n] <- ztemp
    }
    dimnames(stuff) <- list(marker_names,marker_names,1:num_clusts)
    assign(paste("group_",subgrp_idx, sep=""), list(stuff,N))
    dimnames(zstuff) <- list(marker_names,marker_names,1:num_clusts)
    assign(paste("zgroup_",subgrp_idx, sep=""), list(zstuff,N))
    rm(stuff, zstuff, ztemp, bollocks) 
  }
  rm(list = ls(pattern = 'pooled*'))
  
  first_grp = 1
  sec_grp = 2
  
  z1 <- get(paste0("zgroup_",first_grp))[[1]]
  z2 <- get(paste0("zgroup_",sec_grp))[[1]]
  N1 <- get(paste0("zgroup_",first_grp))[[2]]
  N2 <- get(paste0("zgroup_",sec_grp))[[2]]
  
  cur_row <- paste0(group_list[group1],";",group_list[group1])
  
  # compute chi squared values
  temp <- sweep((z1 - z2)^2,MARGIN = 3,(sqrt(1/(N1-3)+1/(N2-3))+S0)^2,"/")
  chi_sq <- apply(temp, 1:2, sum)
  chi_deltas[cur_row,]<-chi_sq[lower.tri(chi_sq)]
  
  return(chi_deltas)
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

# biaxial plots

group_biaxial_plot <- function(grp_idx,marker1,marker2,tube,sub_clusts = seq(1,num_clusts),which = all_adjusts[1],list = short_list,add_line = FALSE,annotate_plot = FALSE){
  count <- 0
  for (file_idx in 1:length(groups[,"Identifier"])){
    if (groups[file_idx,"Tube"] == tube){
      if (groups[file_idx,"Group Number"] == grp_idx){
        loc <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,which,"/",groups[file_idx,"Identifier"])
        for (n in sub_clusts){
          load(paste0(loc,"/node",n,".Rda"))
          if (!is.matrix(tempMat)){ #for size 1 x num_feat frames where dim() won't work
            curMat <- tempMat[c(marker1,marker2)]  # check whether add a comma
          }else{
            curMat <- tempMat[,c(marker1,marker2)]  # check whether add a comma
          }
          #curMat <- tempMat[,c(marker1,marker2)]  # check whether add a comma
          if (count == 0){
            cur_frame <- curMat
            count <- 1
            #cat('Hmmm')
          }else{
            cur_frame <- rbind(cur_frame,curMat)
          }
        }
      }
    }
  }
  smoothScatter(cur_frame[,marker1],cur_frame[,marker2], xlab = marker1, ylab = marker2,ylim = c(-7,7),xlim = c(-7,5),main = list[grp_idx] )
  if(add_line){
    abline(lm(cur_frame[,marker2] ~ cur_frame[,marker1]),col = "orange")
  }
  if(annotate_plot){
    text(-6,6, labels = "*",col = "red",cex = 2)
  }
}


indiv_biaxial_plot <- function(file_idx,marker1,marker2,tube,sub_clusts = seq(1,num_clusts),which = all_adjusts[1]){
  for (n in sub_clusts){
    loc <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,which,"/",groups[file_idx,"Identifier"])
    load(paste0(loc,"/node",n,".Rda"))
    curMat <- tempMat[c(marker1,marker2)]
    if (n == 1){
      cur_frame <- curMat
    }else{
      cur_frame <- rbind(cur_frame,curMat)
    }
  }
  smoothScatter(cur_frame[,marker1],cur_frame[,marker2],xlab = marker1, ylab = marker2, main = paste0(groups[file_idx,"Identifier"],": ",marker1," vs. ",marker2))
}

local_chi_delta_groups_vs_normals <- function(tube,which,sub_clusts = seq(1,num_clusts),s0_val=0.03268146){
  
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
  num_sub_clusts <- length(sub_clusts)

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
  S0 <- rep(s0_val,num_sub_clusts)
  
  
  #fisher z transform
  fisherz <- function(x) {
    0.5*(log(1+x) - log(1-x))
  }
  
  #count which patient in each group we're on
  count <- rep(0,10)
  
  #build bundle arrays
  for (g in 1:num_grps) {
    loc <- paste0(root,"/CyTOF_Analysis/data/",dataset,"/clusters",run_date,"/",tube,which,"/grp_cms/group",g)
    stuff <- array(NA, dim=c(num_feat,num_feat,num_sub_clusts))
    zstuff <- array(NA, dim=c(num_feat,num_feat,num_sub_clusts))
    
    #store number of cells per node. make sure it's at least 4
    N <- get(paste0("counts",tube))[sub_clusts,paste0("Group ",g)]
    N[N<4] = 4
    
    for (n in sub_clusts) {
      clust_idx <- which(n == sub_clusts)
      # create array of correlation coefficients
      bollocks <- paste("grp", g, "_cm_", n, sep='')
      fil <- file.path(loc, paste(bollocks, ".Rda", sep=''))
      load(fil)
      stuff[,,clust_idx] <- get(bollocks)
      # fisher z transform correlation coefficients
      ztemp <- get(bollocks)
      ztemp[row(ztemp)>col(ztemp)] <- fisherz(ztemp[row(ztemp)>col(ztemp)])
      ztemp[row(ztemp) <= col(ztemp)] <- NA # as rho -> 1, z-> 3/inf
      zstuff[,,clust_idx] <- ztemp
      rm(list=bollocks)
    }
    dimnames(stuff) <- list(marker_names,marker_names,1:num_sub_clusts)
    assign(paste("group_", g, sep=""), list(stuff,N))
    dimnames(zstuff) <- list(marker_names,marker_names,1:num_sub_clusts)
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
