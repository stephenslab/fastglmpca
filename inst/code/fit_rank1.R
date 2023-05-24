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
for (iter in 1:10) {
  l0 <- l

  # Solve for f given l.  
  gf <- function (f)
    drop(l %*% (X - exp(tcrossprod(l,f))))
  f <- multiroot(gf,rep(0,m))$root

  # Solve for l given f.
  gl <- function (l)
    drop((X - exp(tcrossprod(l,f))) %*% f)
  l <- multiroot(gl,rep(0,n))$root
  
  cat(max(abs(c(l - l0))),"\n")
}
plot(fit$LL,l,pch = 20)
plot(fit$FF,f,pch = 20)
