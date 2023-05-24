library(Matrix)
library(rootSolve)

# Simulate a 200 x 400 counts matrix.
set.seed(1)
dat <- generate_data_simple(200,400,K = 1)
X <- dat$Y

# Fit a rank-1 GLM-PCA model using the algorithm in the plash package.
fit <- fit_glmpca(X,K = 1,max_iter = 10,algorithm = "ccd",tol = 1e-6)

# Compare the estimates against the ground-truth.
par(mfrow = 1:2)
plot(dat$LL,fit$LL,pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")
plot(dat$FF,fit$FF,pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")

# Here's a simple alternative algorithm for fitting a rank-1 GLM-PCA
# model, x[i,j] ~ Pois(exp(l[i],f[j])).
n <- nrow(X)
m <- ncol(X)
l <- rep(1,n)
exact <- TRUE
numiter <- 4
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
    f <- drop(l %*% (X - 1))/sum(l)
    l <- drop((X - 1) %*% f)/sum(f)
  }
  cat(max(abs(c(l - l0))),"\n")
}

# Rescale l and f.
d <- sqrt(abs(mean(l)/mean(f)))
f <- f*d
l <- l/d

# Compare the two solutions.
plot(fit$LL,l,pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")
plot(fit$FF,f,pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")

#
# a <- 0.1
# plot(x,exp(a*x),type = "l",lwd = 1.25)
# lines(x,1 + a*x,type = "l",col = "dodgerblue",lty = "dotted",lwd = 1.25)
#
