if (!requireNamespace("gumbel")) {
  install.packages("gumbel")
}
library(gumbel)

if (!requireNamespace("dgumbel")) {
  install.packages("dgumbel")
}
library(dgumbel)

library(RPointCloud)

RPointCloud:::dgumbel
RPointCloud:::dlgumbel
args(gumbel::dgumbel)
args(dgumbel::dgumbel)
x <- seq(-4, 2, length = 100)
rp <- RPointCloud:::dgumbel(x, 0, 1)
g <- dgumbel::dgumbel(x, 0, 1)
summary(g - rp) # equal to 15 decimal places
lrp <- RPointCloud:::dlgumbel(x)

RPointCloud:::pgumbel
RPointCloud:::plgumbel
args(gumbel::pgumbel)
args(dgumbel::pgumbel)
q <- seq(0.01, 0.99, 0.01)
rp <- RPointCloud:::pgumbel(q, 0, 1)
g <- dgumbel::pgumbel(q, 0, 1)
summary(g - rp) # equal

RPointCloud:::qgumbel
args(dgumbel::qgumbel)

RPointCloud:::process
