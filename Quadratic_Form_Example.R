library(fields)

n <- c(100,500,1000,1500,2000,2500,3000)
lazy <- vector()
good <- vector()

for(i in (1:7)){
  y <- rnorm(n[i])
  x <- seq(0, 10, length = n[i])
  dists <- rdist(x)
  Sigma <- exp(-dists)
  lazy[i] = system.time(t(y) %*% solve(Sigma) %*% y)[3]
  good[i] = system.time(sum(forwardsolve(t(chol(Sigma)),y)^2))[3]
}

plot(n,lazy,col = "red", xlab = "n", ylab = "time (seconds)", main = "Time to Solve Quadratic Form")
lines(n,lazy,col = "red")

points(n,good,col = "blue")
lines(n,good,col = "blue")

legend(x = "topleft", legend=c("Lazy approach","Cholesky approach"),lty=1,col=c("red","blue"))
