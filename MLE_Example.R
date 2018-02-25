# Simulate the data first:

fi = 10
#fi = 1

ranges_a <- vector()
means_a <- vector()
ranges_v <- vector()
means_v <- vector()

for (i in (1:20)){

  x <- seq(0, fi, length = 50)
  dists <- rdist(x)
  Sigma <- Matern(dists, range = 0.2, smoothness = 1.5)
  Sigma.c = chol(Sigma)
  y <- t(Sigma.c) %*% rnorm(50)
  plot(x,y)
  lines(x,y)

  # MLE calculation

  ell <- function(a){
    Sigma.c <- t(chol(Matern(dists,range=a[1],nu=a[2]))) # lower triangular
    out <- forwardsolve(Sigma.c,y)
    quad.form <- sum(out^2)
    det.part <- 2*sum(log(diag(Sigma.c)))
  
    # up to the normalizing constants
    det.part + quad.form
  }

  out <- optim(par=c(0.5,1),fn=ell,hessian=TRUE,method="L-BFGS-B",lower=c(0.01,0.01),
             upper=c(1,2))

  ## MLEs
  # truth is (0.2, 1.5)
  means_a[i] <- out$par[1]
  means_v[i] <- out$par[2]

  ## MLEs with standard errors
  I.inv <- solve(out$hessian)
  SEs <- sqrt(diag(I.inv))

  CIs <- cbind(out$par - 1.96*SEs,out$par + 1.96*SEs)

  ranges_a[i] <- CIs[1,2] - CIs[1,1]
  ranges_v[i] <- CIs[2,2] - CIs[2,1]
}

print(means_a)
print(ranges_a)
print(means_v)
print(ranges_v)


