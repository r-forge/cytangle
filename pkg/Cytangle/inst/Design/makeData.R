layoutTable <- as.matrix(read.table("layout.table"))
### transform into user coordindates
xform <- function(v) {
  X <- max(v)
  N <- min(v)
  a <- 2/(X-N)
  b <- (N+X)/(N-X)
  function(y) b + a*y
}
f1 <- xform(layoutTable[,1])
u1 <- f1(layoutTable[,1])
f2 <- xform(layoutTable[,2])
u2 <- f2(layoutTable[,2])
ucoords <- cbind(u1, u2)
rm(f1, u1, f2, u2, xform)


save(layoutTable, ucoords, file = "../../data/layoutTable.rda")
save(layoutTable, ucoords, file = "../../R/sysdata.rda")
