
library('Mender')
library('TDA')
source('eb_gumbel_fxns.R')
main <- readRDS('rips_ls/1_rips.rds')


dataset <- main[[J <- sample(100, 1)]]
plot(dataset$diagram)
proc <- Mender:::process(dataset$diagram)

dgm <- dataset$diagram
dgm1 <- dgm[dgm[,1] == 1,]
rats <- dgm1[,3]/dgm1[,2]    # ratio of death to birth
LLR <- log(log(rats))      # log(log(ratio))
beta <- sqrt(6)/pi*sd(LLR) # estimate beta from standard deviation
other <- mean(LLR) + beta*digamma(1)
alpha <- median(LLR) + beta*log(log(2))
zed <- (LLR - alpha)/beta
var(zed) - pi^2/6          # should be zero
mean(zed) + digamma(1)     # should be zero, if no outliers
median(zed) + log(log(2))  # should be zero

ded <- (LLR - other)/beta
var(ded) - pi^2/6          # should be zero
mean(ded) + digamma(1)     # should be zero, if no outliers
median(ded) + log(log(2))  # should be zero

var(proc) - pi^2/6          # should be zero
mean(proc) + digamma(1)     # should be zero, if no outliers
median(proc) + log(log(2))  # should be zero


vals <- EBLGumbel(proc, resn = 200)
hist(vals)
prior <- fitPrior(vals,0,1)
prior
eek <-  postEst(prior, vals)
plot(vals@X0,eek, pch = 16, main = J)
abline(h=c(0,1), col = "gray")

plot(vals@theoretical.pdf/vals@unravel)
abline(h=0, col = "gray")
plot(vals@unravel); abline(h=0, col="gray")

G.proc <- EBLGumbel(proc)
P.proc <- postEst(fitPrior(G.proc, 0, 1), G.proc)
G.zed <- EBLGumbel(zed)
P.zed <- postEst(fitPrior(G.zed, 0, 1), G.zed)
G.ded <- EBLGumbel(ded)
P.ded <- postEst(fitPrior(G.ded, 0, 1), G.ded)

opar <- par(mfrow=c(2,3))
hist(G.proc, main = paste(J, ", process default", sep=""))
hist(G.zed, main = "Median standardized")
hist(G.ded, main = "Mean standardized")
plot(G.proc@X0, P.proc, pch = 16, main = "process default")
abline(h = c(0, 1), col = "gray")
plot(G.zed@X0, P.zed, pch = 16, main = "median std")
abline(h = c(0, 1), col = "gray")
plot(G.ded@X0, P.ded, pch = 16, main = "mean std")
abline(h = c(0, 1), col = "gray")
par(opar)



G.low <- EBLGumbel(proc, method = "loess")
P.low <- postEst(fitPrior(G.low, 0, 1), G.low)
G.den <- EBLGumbel(proc, method = "density")
P.den <- postEst(fitPrior(G.den, 0, 1), G.den)
G.spl <- EBLGumbel(proc, metho = "spline")
P.spl <- postEst(fitPrior(G.spl, 0, 1), G.spl)

opar <- par(mfrow=c(2,3))
hist(G.low, main = paste(J, "; lowess", sep=""))
hist(G.den, main = "density")
hist(G.spl, main = "splines")
plot(G.low@X0, P.low, type = "l", lwd = 2, main = "lowess")
abline(h = c(0, 1), col = "gray")
plot(G.den@X0, P.den, type = "l", lwd = 2, main = "density")
abline(h = c(0, 1), col = "gray")
plot(G.spl@X0, P.spl, type = "l", lwd = 2, main = "spline")
abline(h = c(0, 1), col = "gray")
par(opar)
