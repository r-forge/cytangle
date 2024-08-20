past_max_scale <- function(rips, maxscale, dimension) { #takes in rips list, maxscale int and dimension int
  diag <- rips$diagram
  whichD <- which(diag[,"dimension"] == dimension)
  suppressWarnings(max_death <- max(diag[,3][whichD]))
  suppressWarnings(max_birth <- max(diag[,2][whichD]))
  
  if (max_death >= maxscale | max_birth >= maxscale | max_death == -Inf | max_birth == -Inf) {
    return(TRUE) #if latest death/birth == maxscale, the loop is likely getting cut off by the filtration ending too quickly 
  }
  else {
    return(FALSE)
  }
}

get_max_death <- function(rips, dimension) {
  diag <- rips$diagram
  whichD <- which(diag[,"dimension"] == dimension)
  suppressWarnings(max_death <- max(diag[,3][whichD]))
  if (max_death == -Inf) {
    max_death = NaN
  }
  return(max_death)
  
}

get_max_scale <- function(data) {
  max = 0
  for (columns in colnames(data)) {
    diff <- range(columns)
    if (diff > max) {
      max = diff
    }
  }
  return(max / 2)
  
}