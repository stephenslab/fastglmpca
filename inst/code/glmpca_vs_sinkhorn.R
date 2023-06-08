# Simulate a 100 x 200 counts matrix.
set.seed(1)
n <- 100
m <- 200
dat <- list(l = rnorm(n,1,1/3),f = rnorm(m,1,1/3))
X <- matrix(rpois(n*m,with(dat,exp(tcrossprod(l,f)))),n,m)
hist(X,n = 32)

# Run biwhitening procedure of Landa et al (2021).
numiter <- 20
rows <- rep(1,n)
for (iter in 1:numiter) {
  cols <- 1/colMeans(X * rows)
  rows <- 1/rowMeans(t(t(X) * cols))
}
rows <- 1/rows
cols <- 1/cols

# Approximately fit rank-1 GLM-PCA.
b <- 8
l <- rep(1,n)
for (iter in 1:numiter) {
  f <- drop((l/sum(l)) %*% (X - b))
  l <- drop((X - b) %*% (f/sum(f)))
}

# Compare the solutions.
par(mfrow = c(2,2))
plot(dat$l,l,pch = 20)
plot(dat$f,f,pch = 20)
plot(rows,l,pch = 20)
abline(lm(l ~ rows),col = "magenta",lty = "dashed")
plot(cols,f,pch = 20)
abline(lm(f ~ cols),col = "magenta",lty = "dashed")
