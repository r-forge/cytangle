############################################
## Gumbel Distribution

dgumbel <- function(x, mu, beta) 1/beta*exp((mu-x)/beta)*exp(-exp((mu-x)/beta)) #PDF
pgumbel <- function(q, mu, beta) exp(-exp((mu-q)/beta)) #CDF
qgumbel <- function(p, mu, beta) mu-beta*log(-log(p)) #quantile

GumbelFit <- function(data, resn = 200){
  breaks <- seq(0, max(data), length = resn)
  hist <- hist(data,breaks =breaks, plot = F)
  X0 <- hist$mids
  pdf <- hist$density
  fit <- fitdist(data, "gumbel", start = list(mu = mean(data),beta = sd(data)))
  loc <- as.numeric(fit$estimate["mu"])
  scale <- as.numeric(fit$estimate["beta"])
  val <- list(X0 = X0,pdf = pdf, hist = hist, loc = loc,
              scale = scale, data = data)
  val
}

## Emirical Bayes approach to Gumbel
## 'object' should be obtained from the GumbelFit function
EBGumbel <- function(object, resn = 200) { 
  X0 <- object$X0
  theo <- dgumbel(X0,object$loc, object$scale)
  weird <- theo - object$pdf
  Y0 <- log(object$pdf) - log(theo)
  click <- !is.infinite(Y0)
  Y <- Y0[click]
  X <- X0[click]
  Z <- lm(Y ~ bs(X, df = 5))
  YP <- predict(Z, data.frame(X = X0))
  unravel <- exp(YP + log(theo))
  val <- list(gumbel = object, theoretical.pdf = theo, unravel = unravel)
  val
}

############################################
## Log Gumbel Distribution

plgumbel <- function(x) 1 - exp(-exp(x))
dlgumbel <- function(x) exp(x - exp(x))
lgumbelpdf <- dlgumbel # should be deprecated
lgumbelcdf <- plgumbel # should be deprecated

setClass('EBLGumbel',
         slots = c(X0 = "numeric",
                   pdf = "numeric",
                   gumbel = "numeric",
                   theoretical.pdf = "numeric",
                   unravel = "numeric"))

EBLGumbel <- function(data, resn = 200, method=c("loess", "density", "spline")) {
  # take in processed data
  method <- match.arg(method)
  breaks <- seq(min(data), max(data), length = resn)
  hist <- hist(data, breaks = breaks, plot = FALSE)
  X0 <- hist$mids[hist$mids > -3]
  pdf <- hist$density[hist$mids > -3]
  theo <- lgumbelpdf(X0)
  loessMethod <- function(X0, pdf) {
    LO <- loess(pdf ~ X0)
    predict(LO, X0)
  }
  densityMethod <- function(X0, data) {
    dd <- density(data)
    A <- approx(dd$x, dd$y, X0)
    A$y
  }
  splineMethod <- function(X0, pdf, theo) {
    Y0 <- pdf - theo
    Z <- lm(Y0 ~ bs(X0, df = 5))
    YP <- predict(Z, data.frame(X = X0))
    unravel <- YP + theo
    M <- min(theo)
    unravel[unravel < 0] <- M/2
    unravel
  }
  unravel = switch(method,
                   spline = splineMethod(X0, pdf, theo),
                   density = densityMethod(X0, data),
                   loess = loessMethod(X0, pdf))
  
  val <- new("EBLGumbel",X0 = X0, pdf = pdf,
             gumbel = data,theoretical.pdf = theo, unravel = unravel)
}

setMethod("hist", "EBLGumbel", function(x, xlab="", ylab="Prob(Interesting | X)", main="", ...) {
  top <- max(c(x@unravel * 2, x@`theoretical.pdf` * 2))
  xvals <- x@X0
  hist(x@gumbel, probability=TRUE, breaks=25, ylim=c(0, top),
       xlim=c(min(xvals), max(xvals)), xlab=xlab, main=main)
  lines(xvals, x@`theoretical.pdf`, col="green", lwd=2)
  lines(xvals, x@unravel, col="blue", lwd=2)
  legend("topright", c('Empirical', 'Theoretical'),
         col=c("blue", "green"), lwd=2)
  invisible(x)
})


fitPrior <- function (object, minp, maxp, resn = 1000)  #get minimum p0 to keep postriors nonnegative
{
  sap <- sapply(sip <- seq(minp, maxp, length = resn),
                function(x) min(postEst(x,object)))
  prior <- sip[max(which(sap >= 0))]
  prior
}

postEst <- function(p0, object) {
  1 - p0 * object@theoretical.pdf/object@unravel
}
old_cutoff <- function(target,prior,object) {
  pd <- postEst(prior, object)
  X0  <- object@X0
  cut <- X0[min(which(pd > target))]
  if (is.na(cut)) {
    return(max(X0) + 1)
  }
  else {
    return(cut)
  }

}
new_cutoff <- function(target, prior, object) {
  X0 <- object@X0
  pd <- postEst(prior,object)
  for (index in seq_along(pd)) {
    value <- pd[index]
    if (value >= target) {
      if (X0[index] > 0) {
        cut <- X0[index]
        return(cut)
      }
    }
  }
  return(max(X0 + 1))
}

credible_cutoff <- function(target = 0.99,prior,object) {
  X0 <- object@X0
  pd <- postEst(prior,object)
  top <- quantile(pd,target)
  
  for (index in seq_along(pd)) {
    value <- pd[index]
    if (value >= top){
      if (X0[index] > 1)
        return(X0[index])
    }
  }
}

calculate_scores <- function(true,pred) {
  true_pos = 0
  false_pos = 0
  false_neg = 0
  for (index in 1:length(true)) {
    real <- true[index]
    guess <- pred[index] 
    
    if (real == guess) {
      true_pos = true_pos + real
    }
    else if (guess > real) {
      true_pos = true_pos + real
      false_pos = false_pos + (guess - real)
    }
    else if (real > guess) {
      true_pos = true_pos + guess
      false_neg = false_neg + (real - guess)
    }
  }
  precision = true_pos / (true_pos + false_pos)
  recall = true_pos / (true_pos + false_neg)
  f1 = (2 * precision * recall) / (precision + recall)
  return(list(precision, recall, f1))
}

binary_calculate_scores <- function(true,pred) {
  true_pos = 0
  false_pos = 0
  false_neg = 0
  for (index in 1:length(true)) {
    real <- true[index]
    guess <- pred[index]
    
    if (real == guess) {
      true_pos = true_pos + 1
    }
    else if (guess > real) {
      false_pos = false_pos + 1
    }
    else if (real > guess) {
      false_neg = false_neg + 1
    }
  }
  precision = true_pos / (true_pos + false_pos)
  recall = true_pos / (true_pos + false_neg)
  f1 = (2 * precision * recall) / (precision + recall)
  return(list(precision, recall, f1))
}

calculate_counts <- function(count, sig_count,pred_num,list) {
  if(pred_num == sig_count) {
    list[[1]] <- list[[1]] + sig_count
    list[[3]] <- list[[3]] + (count - sig_count)
  }
  else if (pred_num > sig_count) {
    list[[1]] <- list[[1]] + (sig_count)
    list[[2]] <- list[[2]] + (pred_num - sig_count)
    list[[3]] <- list[[3]] + (count - pred_num)
  }
  else if (sig_count > pred_num) {
    list[[1]] <- list[[1]] + pred_num
    list[[3]] <- list[[3]] + (count - sig_count)
    list[[4]] <- list[[4]] + (sig_count - pred_num)
  }
  return(list)
}
