CMD <-function(p,q){
  # http://publik.tuwien.ac.at/files/pub-et_10119.pdf
  n_p <- norm(p,type = c("f"))
  n_q <- norm(q,type = c("f"))
  tr_pq <- sum(diag(p %*% q))
  dis <- 1 - tr_pq/(n_p*n_q)
  return(dis)
}

CholD <-function(p,q){
  P <- chol(p,pivot = TRUE) #t(chol(p))
  Q <- chol(q,pivot = TRUE) #t(chol(p))
  #covDif <- sqrt(t(chol(cov(p))-chol(cov(q)))*(chol(cov(p))-chol(cov(q)))) # http://stats.stackexchange.com/questions/71236/covariance-and-correlation-matrix-comparison
  #dis <- sqrt(sum(diag(cov(covDif))))
  radicand <- t(P-Q) %*% (P-Q) # 4/17/16: I don't think the pivoting will hurt this step
  dis <- sqrt(sum(diag(radicand))) # http://projecteuclid.org/euclid.aoas/1254773280 
  # http://projecteuclid.org/download/pdfview_1/euclid.aoas/1254773280
  return(dis)
}

DetD <-function(p,q){
  #http://stats.stackexchange.com/questions/14673/measures-of-similarity-or-distance-between-two-covariance-matrices
  # compute determinants
  dp <- det(p)
  dq <- det(q)
  d_pq <- det((p+q)/2)
  # positivity checks
  dp   <- be_positive(dp)
  dq   <- be_positive(dq)
  d_pq <- be_positive(d_pq)
  # norm
  if ( all(p == q) | ((dp == 0) & (dq == 0) & (d_pq == 0)) ){
    dis <- 0
  }else{
    dis <- log(d_pq/sqrt(dp*dq)) # From arithmetic geometric mean
  }
  return(dis)
}

CondD <-function(p,q){
  #http://stats.stackexchange.com/questions/14673/measures-of-similarity-or-distance-between-two-covariance-matrices
  p.eig <- eigen(p,symmetric = TRUE)
  p.eig$values <- be_positive(p.eig$values)
  #http://realizationsinbiostatistics.blogspot.com/2008/08/matrix-square-roots-in-r_18.html
  p.sqrt <- p.eig$vectors %*% diag(sqrt(p.eig$values)) %*% t(p.eig$vectors)
  conjp <- solve(p.sqrt)
  # going this route might be problematic if there's a zero eigenvalue 
  sigma <- conjp %*% q %*% conjp
  sigma_sym <- (sigma+t(sigma))/2 # To prevent numerical error from breaking symmetry
  sigma.eig <- eigen(sigma_sym,symmetric = TRUE)
  evalues <- sigma.eig$values
  pos_evalues <- evalues[evalues>0]
  lam_max <- max(pos_evalues)
  lam_min <- min(pos_evalues)
  dis <- log(lam_max) - log(lam_min)
  return(dis)
}