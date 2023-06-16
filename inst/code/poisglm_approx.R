# TO DO: Explain here what this script is for, and how to use it.
library(pracma)

# Simulate data from a Poisson GLM model, y ~ Pois(r), r = exp(x*b).
n <- 100
x <- rnorm(n)
b <- -1
r <- exp(x*b)
y <- rpois(n,r)

# Plot the log-likelihood surface.
compute_loglik <- function (x, y, b)
  sum(x*y*b - exp(x*b))
b <- seq(-4,4,length.out = 1000)
n <- length(b)
ll <- rep(0,n)
for (i in 1:n)
  ll[i] <- compute_loglik(x,y,b[i])
plot(b,ll,type = "l",col = "black",lwd = 2)
