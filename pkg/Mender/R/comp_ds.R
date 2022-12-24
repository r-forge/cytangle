# Downsample function
compDensity <- function(dset, sigma = 1) {
  EG <- dset
  weighted <- rep(0, nrow(EG))
  for (i in 1:nrow(dset)) {
    mu <- dset[i,]
    euc <- apply(sweep(EG, 2, mu, '-')^2, 1, sum)
    spread <- exp(-euc/(2*sigma^2))/(2*pi*sigma)^(3/2)
    weighted <- weighted + spread
  }
  weighted
}

downSample <- function(dset, targetNum = 100, sigma = 1) {
  dense <- compDensity(dset, sigma)
  sdense <- sort(dense)
  cdf <- rev(cumsum(1.0 / rev(sdense) ))
  targets <- (targetNum - 1 : length(sdense)) / cdf 
  boundary <- targets[which.min(targets - sdense > 0)]
  selector <- boundary/dense > runif(length(dense))
  keepers <- usable[selector,]
  keepers
}