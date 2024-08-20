if (!requireNamespace("gumbel")) {
  install.packages("gumbel")
}
library(gumbel)

if (!requireNamespace("dgumbel")) {
  install.packages("dgumbel")
}
library(dgumbel)

library(EBGLIDE)

EBGLIDE:::dgumbel
EBGLIDE:::dlgumbel
args(gumbel::dgumbel)
args(dgumbel::dgumbel)
x <- seq(-5, 5, length = 500)
rp <- EBGLIDE:::dgumbel(x, 0, 1)
g <- dgumbel::dgumbel(x, 0, 1)
summary(g - rp) # equal to 15 decimal places
lrp <- EBGLIDE:::dlgumbel(x)
par(lwd=2)
plot(x, lrp, type="l")
lines(x, rp, col = "orange")
legend("topleft", c("LeftSkewed", "RightSkewed"), lwd=2,
       col = c("black", "orange"))XS

EBGLIDE:::pgumbel
EBGLIDE:::plgumbel
args(gumbel::pgumbel)
args(dgumbel::pgumbel)
q <- seq(0.01, 0.99, 0.01)
rp <- EBGLIDE:::pgumbel(q, 0, 1)
g <- dgumbel::pgumbel(q, 0, 1)
summary(g - rp) # equal

EBGLIDE:::qgumbel
args(dgumbel::qgumbel)

EBGLIDE:::process
