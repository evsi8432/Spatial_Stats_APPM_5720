library(fields)

n <- c(100,500,1000,1500,2000,2500,3000)
lazy <- vector()
good <- vector()

for(i in (1:7)){
  x <- seq(0, 10, length = n[i])
  dists <- rdist(x)
  Sigma <- exp(-dists)
  lazy[i] = system.time(determinant(Sigma,log = TRUE)$mod)[3]
  good[i] = system.time(2*sum(log(diag(chol(Sigma)))))[3]
}

plot(n,lazy,col = "red", xlab = "n", ylab = "time (seconds)", main = "Time to Solve Determinant")
lines(n,lazy,col = "red")

points(n,good,col = "blue")
lines(n,good,col = "blue")

legend(x = "topleft", legend=c("Lazy approach","Cholesky approach"),lty=1,col=c("red","blue"))