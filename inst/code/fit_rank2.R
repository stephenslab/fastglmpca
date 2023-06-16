library(Matrix)
library(rootSolve)

# Simulate a 100 x 200 counts matrix.
set.seed(1)
n <- 100
m <- 200
dat <- list(LL = matrix(rnorm(2*n,1,0.3),2,n),
            FF = matrix(rnorm(2*m,1,0.3),2,m))
X <- matrix(rpois(n*m,with(dat,exp(crossprod(LL,FF)))),n,m)

# Fit a rank-2 GLM-PCA model using the algorithm in the plash
# package. (To nudge GLM-PCA model toward the truth, here I initialize
# the factorization to the true values.)
fit0 <- init_glmpca(X,LL = dat$LL,FF = dat$FF)
fit <- fit_glmpca(X,fit0 = fit0,max_iter = 40,algorithm = "ccd",tol = 1e-6)

# Compare the estimates against the ground-truth.
par(mfrow = 1:2)
plot(dat$LL,fit$LL,pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")
plot(dat$FF,fit$FF,pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")

# Here's a simple alternative algorithm for fitting a rank-2 GLM-PCA
# model, x[i,j] ~ Pois(exp(L[i,1]*F[j,1] + L[i,2]*F[j,2])). (Note that
# this code should work for any K.)
K <- 2
n <- nrow(X)
m <- ncol(X)
L <- matrix(dat$LL[1,])
F <- matrix(dat$FF[1,])
exact <- FALSE
numiter <- 4
t0 <- proc.time()
for (k in 1:2) {
  cat("k =",k,"\n")
  if (k > 1) {
    L <- cbind(L,dat$LL[,k])
    F <- cbind(F,dat$FF[,k])
  }
  for (iter in 1:numiter) {
    L0 <- L
    if (exact) {

      # Solve for F[,k] given L, F[,-k] exactly.  
      gf <- function (f) {
        Fnew     <- F
        Fnew[,k] <- f
        return(drop(L[,k] %*% (X - exp(tcrossprod(L,Fnew)))))
      }
      F[,k] <- multiroot(gf,start = rep(0,m))$root
      
      # Solve for L[,k] given F, L[,-k] exactly.
      gl <- function (l) {
        Lnew     <- L
        Lnew[,k] <- l
        return(drop((X - exp(tcrossprod(Lnew,F))) %*% F[,k]))
      }
      L[,k] <- multiroot(gl,start = rep(0,n))$root
    } else {

      # Solve for F[,k] given L, F[,-k] approximately.
      # Then solve for L[,k] given F, L[,-k] approximately.
      if (k == 1) {
        B <- matrix(1,n,m)
      } else {
        ks <- seq(1,k-1)
        B <- tcrossprod(L[,ks,drop = FALSE],F[,ks,drop = FALSE])
      } 
      F[,k] <- drop((L[,k] %*% (X - B))/(L[,k]^2 %*% B))
      L[,k] <- drop(((X - B) %*% F[,k])/(B %*% F[,k]^2))
    }
    cat(max(abs(c(L - L0))),"\n")
  }
}
t1 <- proc.time()
print(t1 - t0)

# Compare the two solutions.
par(mfrow = c(2,2))
x <- fit$LL[1,]
y <- L[,1]
out <- lm(y ~ x)
plot(x,y,pch = 20)
abline(a = coef(out)[1],b = coef(out)[2],lty = "dashed",col = "magenta")

x <- fit$FF[1,]
y <- F[,1]
out <- lm(y ~ x)
plot(x,y,pch = 20)
abline(a = coef(out)[1],b = coef(out)[2],lty = "dashed",col = "magenta")

x <- fit$LL[2,]
y <- L[,2]
out <- lm(y ~ x)
plot(x,y,pch = 20)
abline(a = coef(out)[1],b = coef(out)[2],lty = "dashed",col = "magenta")

x <- fit$FF[2,]
y <- F[,2]
out <- lm(y ~ x)
plot(x,y,pch = 20)
abline(a = coef(out)[1],b = coef(out)[2],lty = "dashed",col = "magenta")
