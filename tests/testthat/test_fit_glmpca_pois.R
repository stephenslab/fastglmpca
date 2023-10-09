context("fit_glmpca_pois")

test_that("Some basic tests of fit_glmpca_pois",{

  # Simulate a 100 x 200 data set to factorize.
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y

  # Fit a GLM-PCA model to the data.
  fit0 <- init_glmpca_pois(Y,K = 3)
  suppressWarnings(capture.output(
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,max_iter = 20)))
  suppressWarnings(capture.output(
    fit <- fit_glmpca_pois(Y,fit0 = fit_quick,max_iter = 500,tol = 1e-8,
                           control = list(calc_deriv = TRUE,
                                          calc_max_diff = TRUE))))
  capture.output(print(summary(fit)))

  # All the updates should monotonically increase the likelihood.
  expect_nondecreasing(fit$progress$loglik)
  
  # The fixed row and column intercepts should not change.
  expect_equal(fit0$X,fit$X)
  expect_equal(fit0$B,fit$B)
  expect_equal(fit0$Z,fit$Z)
  expect_equal(fit0$W,fit$W)
  
  # Check that orthogonality constraints are satisfied.
  d <- diag(fit$D)
  expect_equivalent(crossprod(fit$U),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit$V),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(d,sort(d,decreasing = TRUE))
})

test_that("fit_glmpca_pois works with K = 1",{

  # Simulate a 100 x 200 data set to factorize.
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 1)$Y

  # Fit a GLM-PCA model to the data.
  fit0 <- init_glmpca_pois(Y,K = 1)
  suppressWarnings(capture.output(
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,max_iter = 20)))
  suppressWarnings(capture.output(
    fit <- fit_glmpca_pois(Y,fit0 = fit_quick,max_iter = 500,tol = 1e-8,
                           control = list(calc_deriv = TRUE,
                                          calc_max_diff = TRUE))))
  capture.output(print(summary(fit)))

  # All the updates should monotonically increase the likelihood.
  expect_nondecreasing(fit$progress$loglik)

  # The fixed row and column intercepts should not change.
  expect_equal(fit0$X,fit$X)
  expect_equal(fit0$B,fit$B)
  expect_equal(fit0$Z,fit$Z)
  expect_equal(fit0$W,fit$W)

  # Check that orthogonality constraints are satisfied.
  expect_equivalent(crossprod(fit$U),1,scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit$V),1,scale = 1,tolerance = 1e-8)
})

test_that("fit_glmpca_pois works with orthonormalize = FALSE",{

  # Simulate a 100 x 200 data set to factorize.
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y

  # Fit a GLM-PCA model to the data.
  fit0 <- init_glmpca_pois(Y,K = 3)
  suppressWarnings(capture.output(
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,max_iter = 20)))
  suppressWarnings(capture.output(
    fit <- fit_glmpca_pois(Y,fit0 = fit_quick,max_iter = 500,tol = 1e-8,
                           control = list(calc_deriv = TRUE,
                                          calc_max_diff = TRUE))))
  capture.output(print(summary(fit)))

  # All the updates should monotonically increase the likelihood.
  expect_nondecreasing(fit$progress$loglik)

  # The fixed row and column intercepts should not change.
  expect_equal(fit0$X,fit$X)
  expect_equal(fit0$B,fit$B)
  expect_equal(fit0$Z,fit$Z)
  expect_equal(fit0$W,fit$W)

  # Check that orthogonality constraints are satisfied.
  d <- diag(fit$D)
  expect_equivalent(crossprod(fit$U),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit$V),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(d,sort(d,decreasing = TRUE))
})

test_that("Fit works with no row intercept or column size factor", {
  
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y
  
  # Fit a GLM-PCA model to the data.
  fit0 <- init_glmpca_pois(
    Y,
    K = 3, 
    col_size_factor = FALSE, 
    row_intercept = FALSE
  )
  
  suppressWarnings(capture.output(
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,max_iter = 20)))
  
  expect_nondecreasing(fit_quick$progress$loglik)
  
  # The fixed row and column intercepts should not change.
  expect_equal(fit0$X,fit_quick$X)
  expect_equal(fit0$B,fit_quick$B)
  expect_equal(fit0$Z,fit_quick$Z)
  expect_equal(fit0$W,fit_quick$W)
  
  # Check that orthogonality constraints are satisfied.
  expect_equivalent(crossprod(fit_quick$U),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit_quick$V),diag(3),scale = 1,tolerance = 1e-8)
  
})

test_that("Final fit is the same with sparse and dense Y", {
  
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y
  
  Y_sp <- as(Y, "sparseMatrix")
  
  set.seed(1)
  # Initialize a GLM-PCA model to the data with dense Y
  fit0 <- init_glmpca_pois(
    Y,
    K = 3
  )
  
  suppressWarnings(capture.output(
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,max_iter = 20)))
  
  set.seed(1)
  # Initialize a GLM-PCA model to the data with sparse Y
  fit0_sp <- init_glmpca_pois(
    Y_sp,
    K = 3
  )  
  
  suppressWarnings(capture.output(
    fit_quick_sp <- fit_glmpca_pois(Y_sp,fit0 = fit0_sp,max_iter = 20)))
  
  expect_equal(fit_quick$X,fit_quick_sp$X)
  expect_equal(fit_quick$B,fit_quick_sp$B)
  expect_equal(fit_quick$Z,fit_quick_sp$Z)
  expect_equal(fit_quick$W,fit_quick_sp$W)
  expect_equal(fit_quick$U,fit_quick_sp$U)
  expect_equal(fit_quick$V,fit_quick_sp$V)
  expect_equal(fit_quick$D,fit_quick_sp$D)
  expect_equal(fit_quick$loglik, fit_quick_sp$loglik)
  
})

test_that("Test fit works with input covariates",{
  
  # Simulate a 100 x 200 data set to factorize.
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y
  X <- matrix(data = rnorm(n * 2), nrow = n, ncol = 2)
  Z <- matrix(data = rnorm(m * 2), nrow = m, ncol = 2)
  
  # Fit a GLM-PCA model to the data.
  fit0 <- init_glmpca_pois(Y,K = 3, X = X, Z = Z, fixed_b_cols = c(1))
  
  suppressWarnings(capture.output(
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,max_iter = 20)))
  
  expect_nondecreasing(fit_quick$progress$loglik)
  
  expect_equal(fit0$X,fit_quick$X)
  expect_equal(fit0$Z,fit_quick$Z)
  
  expect_equal(fit0$B[, fit0$fixed_b_cols],fit_quick$B[, fit0$fixed_b_cols])
  expect_equal(fit0$W[, fit0$fixed_w_cols],fit_quick$W[, fit0$fixed_w_cols])
  
  # Check that orthogonality constraints are satisfied.
  expect_equivalent(crossprod(fit_quick$U),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit_quick$V),diag(3),scale = 1,tolerance = 1e-8)
  
})
