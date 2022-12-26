## Copright (C) 2022 Kevin R. Coombes, RB McGee, and Jake Reed

## helper function to get cycle-centroid
getCentroid <- function(cycle) {
  x <- as.vector(c(cycle[,1,1], cycle[,2,1]))
  y <- as.vector(c(cycle[,1,2], cycle[,2,2]))
  c(x = mean(x), y = mean(y))
}
## helper function to compute angles of data points around centroid
getAngles <- function(dset, centroid) {
  recentered <- sweep(dset, 2, centroid, "-") # centroid now at (0,0)
  atan2(recentered[,1], recentered[,2])
}

## Function to calculate mean and SD of section of graph starting from
## centroid
angleMean <- function(dset, rips, lb_angle = 0, ub_angle = 30, incr = 15) {
  # get the longest 1D cycle.
  target_cycle <- getLongCycle(rips)
  ## Pull out x and y coordinates from cycle and create centroid
  centroid <- getCentroid(target_cycle)
  ## Find angles with respect to centroid for all points
  radians <- getAngles(dset, centroid)
  alpha_v <- 180 + radians*180/pi
  ## Assign angles to dset data frame
  dset$angle <- alpha_v
  ## For loop calculating mean and sd of different angles
  result.df <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(result.df) <- c("AngleSpread", "UpperBound", "Mean", "SD")
  lb <- lb_angle
  for(ub in seq(ub_angle, 360, incr)) {
    angle_spread <- ub - lb
    set <- subset(dset, angle > lb & angle < ub)
    m.gene <- mean(set[, 3])
    sd.gene <- sd(set[, 3])
    result.row <- data.frame(AngleSpread = angle_spread, UpperBound = ub, 
                             Mean = m.gene, SD = sd.gene)
    result.df <- rbind(result.df, result.row)
    lb <- lb + incr
  }
  return(result.df)
}
