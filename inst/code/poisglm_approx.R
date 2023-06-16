# TO DO: Explain here what this script is for, and how to use it.
library(pracma)

# Simulate data from a Poisson GLM model,
#
#   y ~ Pois(r)
#   r = exp(x*b)
#
set.seed(1)
n <- 80
x <- abs(rnorm(n))
b <- -1
r <- exp(x*b)
y <- rpois(n,r)

# Plot the log-likelihood surface.
compute_loglik <- function (b, x, y)
  sum(x*y*b - exp(x*b))
b <- seq(-5,0.75,length.out = 1000)
n <- length(b)
ll <- rep(0,n)
for (i in 1:n)
  ll[i] <- compute_loglik(b[i],x,y)
plot(b,ll - max(ll),type = "l",col = "black",lwd = 2,ylab = "loglik")
b.mle <- optimize(compute_loglik,c(-4,4),maximum = TRUE,x = x,y = y)
points(b.mle$maximum,0,pch = 4,col = "black")

# Plot the approximate log-likelihood surface.
#
#   ans <- polyApprox(exp,a = -4,b = 0,n = 2)
#   a1  <- ans$p[2]
#   a2  <- ans$p[1]
#
a1 <- 0.5594
a2 <- 0.08637
f <- function (x)
  a1*x + a2*x^2
compute_loglik_approx <- function (b, x, y) {
  i <- which(y > 0)
  j <- which(y == 0)
  return(compute_loglik(b,x[i],y[i]) - sum(f(b*x[j])))
}
n <- length(b)
ll1 <- rep(0,n)
for (i in 1:n)
  ll1[i] <- compute_loglik_approx(b[i],x,y)
lines(b,ll1 - max(ll1),type = "l",col = "red",lwd = 2)
b.approx <- optimize(compute_loglik_approx,c(-4,4),maximum = TRUE,x = x,y = y)
points(b.approx$maximum,0,pch = 4,col = "red")
