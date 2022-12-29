## helper function to get cycle-centroid
getCentroid <- function(cycle, view) {
  pts <- cycleSupport(cycle, view)
  apply(pts, 2, mean)
}

## helper function to compute angles of data points around centroid
getAngles <- function(dset, centroid) {
  recentered <- sweep(dset, 2, centroid, "-") # centroid now at (0,0)
  atan2(recentered[,1], recentered[,2])
}

## Function to calculate mean and SD of section of graph starting from
## centroid
angleMean <- function(view, rips, cycle = NULL, angleWidth = 20, incr = 15) {
  if (is.null(cycle)) {
    cycle <- getCycle(rips, 1) # longest cycle
  }
  ctr <- getCentroid(cycle, view)
  theta <- getAngles(view, ctr)
  degrees <- theta*180/pi
  deg <- c(degrees, degrees + 360)
  od <- order(deg)
  deg <- deg[od]
  magic <- rbind(view, view)[od,]
  partition <- seq(0, 360 - incr, incr)
  GM <- t(sapply(partition, function(center) {
    lb <- center - 10
    ub <- center + 10
    set <- subset(magic, deg > lb & deg < ub)
    m.gene <- apply(set, 2, mean, na.rm = TRUE)
    m.gene
  }))
  rownames(GM) <- partition
  GM
}

if (FALSE) {
  angle.df <- GM[, pal <- c("Ki-67", "PCNA", "CycB", "pRb", "CycA", "CD99")]
  opar <- par(cex = 1.5, mai = c(0, 0, 0, 0))
  circos.clear()
  circos.par(track.height = 0.08, 
             start.degree = 90)
  ## For loop for each track/gene
  for(i in 1:length(angle.df[1,])) {
    ## If statement to add angle designations on outside of 
    ## first track only
    if(i == 1) {
      data <- as.data.frame(angle.df[,i])
      col <- colorRamp2(c(min(data), max(data)), c("white", col.ls[i]))
      circos.heatmap(data, rownames.side = "outside", 
                     cluster = FALSE, col = col, rownames.cex = 1)
    } else {
      data <- as.data.frame(angle.df[,i])
      col <- colorRamp2(c(min(data), max(data)), c("white", col.ls[i]))
      circos.heatmap(data, cluster = FALSE, col = col)
    }
  }
  ## Generate legend
  lgd <- Legend(at = colnames(angle.df), title = "Genes", type = "points",
                title_position = "topleft", legend_gp = gpar(col = col.ls))
  ## Draw legend
  draw(lgd, x = unit(0.97, "npc"), y = unit(0.05, "npc"), just = c("right", "bottom"))
  par(opar)
}
