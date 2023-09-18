context("init_glmpca_pois")

test_that("Initial fits should satisfy orthogonality constraints",{
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y
  X <- matrix(rnorm(2*n),n,2)
  Z <- matrix(rnorm(m),m,1)

  # Initialize a GLM-PCA fit with K = 3.
  fit <- init_glmpca_pois(Y,X = X,Z = Z,K = 3)
  d   <- diag(fit$D)
  expect_equivalent(crossprod(fit$U),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit$V),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(d,sort(d,decreasing = TRUE))
  
  # Initialize a GLM-PCA fit with K = 1.
  fit <- init_glmpca_pois(Y,X = X,Z = Z,K = 1)
  expect_equivalent(crossprod(fit$U),1,scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit$V),1,scale = 1,tolerance = 1e-8)
})
