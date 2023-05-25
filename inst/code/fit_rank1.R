library(Matrix)
library(rootSolve)

# Simulate a 200 x 400 counts matrix.
set.seed(1)
n <- 200
m <- 400
dat <- list(LL = rnorm(n,1,0.3),FF = rnorm(m,1,0.3))
X <- matrix(rpois(n*m,with(dat,exp(tcrossprod(LL,FF)))),n,m)

# Fit a rank-1 GLM-PCA model using the algorithm in the plash package.
fit <- fit_glmpca(X,K = 1,max_iter = 10,algorithm = "ccd",tol = 1e-6)

# Compare the estimates against the ground-truth.
par(mfrow = 1:2)
plot(dat$LL,-fit$LL,pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")
plot(dat$FF,-fit$FF,pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")

# Here's a simple alternative algorithm for fitting a rank-1 GLM-PCA
# model, x[i,j] ~ Pois(exp(l[i],f[j])).
n <- nrow(X)
m <- ncol(X)
l <- rep(1,n)
exact <- FALSE
numiter <- 4
t0 <- proc.time()
b <- 8
for (iter in 1:numiter) {
  l0 <- l
  if (exact) {
  
    # Solve for f given l exactly.  
    gf <- function (f)
      drop(l %*% (X - exp(tcrossprod(l,f))))
    f <- multiroot(gf,start = rep(0,m))$root

    # Solve for l given f exactly.
    gl <- function (l)
      drop((X - exp(tcrossprod(l,f))) %*% f)
    l <- multiroot(gl,start = rep(0,n))$root
  } else {

    # Solve for f given l then l given f approximately.
    u <- l/sum(l)  
    f <- drop(u %*% (X - b))
    u <- f/sum(f)
    l <- drop((X - b) %*% u)
  }
  cat(max(abs(c(l - l0))),"\n")
}
t1 <- proc.time()
print(t1 - t0)

# Rescale l and f.
d <- sqrt(abs(mean(l)/mean(f)))
f <- f*d
l <- l/d

# Compare the two solutions.
plot(fit$LL,-l,pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")
plot(fit$FF,-f,pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")

#
# a <- 0.1
# plot(x,exp(a*x),type = "l",lwd = 1.25)
# lines(x,1 + a*x,type = "l",col = "dodgerblue",lty = "dotted",lwd = 1.25)
#

# This implements the biwhitening procedure as described in
# Algorithm 1 of Landa et al (2021).
scale.rows <- function (A, b)
  A * b
scale.cols <- function (A, b)
  t(t(A) * b)
rows <- rep(1,n)
for (iter in 1:numiter) {
  cols <- 1/colMeans(scale.rows(X,rows))
  rows <- 1/rowMeans(scale.cols(X,cols))
}
plot(1/rows,l,pch = 20)
plot(1/cols,f,pch = 20)
print(cor(1/rows,l))
print(cor(1/cols,f))
