
s0_comp <- function (resid, sd, s0.perc = seq(0, 1, by = 0.025)){
  br = unique(quantile(sd, seq(0, 1, len = 101))) #quantiles
  nbr = length(br)
  a <- cut(sd, br, labels = F) # place standard deviations into the quantiles
  a[is.na(a)] <- 1 #what are the nodes causing the NAs. Are they big, small? maybe set them to nbr-1?
  cv.sd <- rep(0, length(s0.perc)) #storage for a CoeffVar for each candidate s0
  mu.sd <- cv.sd
  sd.sd <- cv.sd
  #print(s0.perc)
  for (j in 1:length(s0.perc)) { #for each candidate, compute CoeffVar
    w <- quantile(sd, s0.perc[j])
    w[j == 1] <- 0
    w_vec <- w*rep(1,length(sd))
    tt2 <- resid^2 * 1/(sd + w_vec)^2#tt * sd/(sd + w)
    tt2[tt2 == Inf] = NA
    sds <- rep(0, nbr - 1) #storage for a mad for each quantile
    for (i in 1:(nbr - 1)) { #for each quantile, compute mad of the modified statistic
      sds[i] <- mad(tt2[a == i], na.rm = TRUE)
    }
    hist(sds)
    print('mean')
    print(mu.sd[j]<-mean(sds,na.rm = TRUE))
    print('std. dev.')
    print(sd.sd[j]<-sqrt(var(sds,na.rm = TRUE)))
    cv.sd[j] <- sqrt(var(sds, na.rm = TRUE))/mean(sds, na.rm = TRUE)
  }
  candidates <- as.vector(quantile(sd[sd != 0], s0.perc))
  constr_cvs <- cv.sd[cv.sd+(1/0.97)*candidates<1]
  o = (1:length(s0.perc))[constr_cvs == min(constr_cvs)]
  #o = (1:length(s0.perc))[cv.sd == min(cv.sd)]
  s0.hat = quantile(sd[sd != 0], s0.perc[o])
  return(list(s0.perc = s0.perc, cv.sd = cv.sd, s0.hat = s0.hat))
}

candidates <- as.vector(quantile(sd[sd != 0], s0.perc))
plot(candidates,cv.sd)
Browse[2]> lines(candidates,1-(1/0.97)*candidates)

set.seed(3221987)
g <- function(s_alpha,delta,sigma){abs(mean((delta/(sigma+s_alpha))^2)-1)}
v <- function(s_alpha,delta,sigma){abs(var((delta/(sigma+s_alpha))^2)-2)}

num_rows <- sample(dim(all_rs)[1],1e4) # num rows
num_cols <- sample(dim(all_sigmas)[2],1e3)

rs_cur <- as.vector(all_rs[num_rows,num_cols])
sigs_cur <- as.vector(all_sigmas[num_rows,num_cols])

min_loc_via_mean <- optimize(g,c(0,2),delta=rs_cur,sigma=sigs_cur)
min2_loc_via_mean <- optimize(g,c(0,2),delta=rs_cur,sigma=sigs_cur,tol=1e-14)
min_glob_via_mean <- optimize(g,c(0,2),delta=as.vector(rs_cur),sigma=as.vector(sigs_cur),tol=1e-10)

min_via_var <- optimize(g,c(0,2),delta=rs_cur,sigma=sigs_cur)


alot_pts <- sapply(seq(0,2,0.01),g,delta=rs_cur,sigma=sigs_cur)
plot(seq(0,2,0.02),alot_pts,ylim = c(0.9998,1))
plot(seq(0,2,0.01),alot_pts)
