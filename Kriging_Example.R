################################################################################################
## Kriging uncertainty: 1D (simple kriging)
################################################################################################

library(fields)

##
## Generate some 1-d observations Y(s) = Z(s), E Z(s) = 0
##

n <- 1000

grd <- seq(0,20,length.out=n)
dist.mat <- rdist(grd)

## Matern: a=1, nu=1.5
Sigma <- Matern(dist.mat,range=0.5,nu=1.5)
Sigma.c <- chol(Sigma)

z <- t(Sigma.c) %*% rnorm(n)

plot(z~grd,type="b")

## Observations
ns <- 15
#these <- seq(1,n,length = 15) # uniform obs (a)
these <- sort(sample(1:n,ns)) # random places (b)
obs.grd <- grd[these]
y <- z[these]

plot(y~obs.grd,pch=19,xlim=range(grd),ylim=range(z))
lines(z~grd)

##
## Simple kriging with predictive uncertainty
##

## Simple kriging
Sigma <- Sigma[these,these]
dim(Sigma)

dist.to.pred <- rdist(obs.grd,grd)
Sigma0 <- Matern(dist.to.pred,range=0.5,nu=1.5)
dim(Sigma0)

y.hat <- t(Sigma0) %*% solve(Sigma) %*% y

# Kriging variance and standard error
y.hat.var <- 1 - diag(t(Sigma0) %*% solve(Sigma) %*% Sigma0)
y.hat.var[y.hat.var < 0] <- 0
y.hat.se <- sqrt(y.hat.var)

y.hat.up <- y.hat + 1.96 * y.hat.se
y.hat.dn <- y.hat - 1.96 * y.hat.se

plot(y~obs.grd,pch=19,xlim=range(grd),ylim=range(y.hat.up,y.hat.dn),
     main = "Krieged Estimate of Isotropic Process",
     xlab = "t", ylab = "Z")
polygon(x=c(grd,rev(grd)),y=c(y.hat.up,rev(y.hat.dn)),col="grey")
points(y~obs.grd,pch=19)

lines(z~grd,col="blue")
lines(y.hat~grd)

legend(x=locator(1),legend=c("Real Process","Kriged Estimate","95% Interval"),lty=1,col=c("Blue","black","grey"))