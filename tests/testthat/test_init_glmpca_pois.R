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
  expect_equivalent(crossprod(fit$U),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit$V),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(fit$d,sort(fit$d,decreasing = TRUE))
  
  # Initialize a GLM-PCA fit with K = 1.
  fit <- init_glmpca_pois(Y,X = X,Z = Z,K = 1)
  expect_equivalent(crossprod(fit$U),1,scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit$V),1,scale = 1,tolerance = 1e-8)
})

test_that("Initial fit works with no row intercept or column size factor",{
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y
  
  # Initialize a GLM-PCA model to the data.
  fit0 <- init_glmpca_pois(Y,K = 3, col_size_factor = FALSE,
                           row_intercept = FALSE)        
  
  expect_equal(fit0$X,numeric(0))
  expect_equal(fit0$B,numeric(0))
  expect_equal(fit0$Z,numeric(0))
  expect_equal(fit0$W,numeric(0))
})

test_that("Initial fit is the same with sparse and dense Y",{
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y
  Y_sp <- as(Y, "sparseMatrix")
  
  # Initialize a GLM-PCA model to the data with dense Y.
  set.seed(1)
  fit0 <- init_glmpca_pois(Y,K = 3)
  
  set.seed(1)
  # Initialize a GLM-PCA model to the data with sparse Y.
  fit0_sp <- init_glmpca_pois(Y_sp,K = 3)
  
  expect_equal(fit0,fit0_sp)
})
