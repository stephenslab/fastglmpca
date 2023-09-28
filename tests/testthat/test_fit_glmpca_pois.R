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
    fit <- fit_glmpca_pois(Y,fit0 = fit0,
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
    fit <- fit_glmpca_pois(Y,fit0 = fit0,
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
