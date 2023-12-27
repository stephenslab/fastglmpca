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
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,
                                 control = list(maxiter = 20,
                                                calc_deriv = TRUE,
                                                calc_max_diff = TRUE))))
  suppressWarnings(capture.output(
    fit <- fit_glmpca_pois(Y,fit0 = fit_quick,
                           control = list(maxiter = 500,tol = 1e-8,
                                          calc_deriv = TRUE,
                                          calc_max_diff = TRUE))))
  capture.output(print(summary(fit)))

  # All the updates should monotonically increase the likelihood.
  expect_nondecreasing(fit$progress$loglik)

  # The fixed row and column intercepts should not change.
  expect_equal(fit0$X,fit$X)
  expect_equal(fit0$B,fit$B)
  expect_equal(fit0$Z,fit$Z)

  # Check that orthogonality constraints are satisfied.
  expect_equivalent(crossprod(fit$U),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit$V),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(fit$d,sort(fit$d,decreasing = TRUE))

  # loglik, tail(progress$loglik,n = 1) and manual calculation of the
  # log-likelihood should all be the same.
  loglik <-
    lik_glmpca_pois_log(Y,
                        LL = with(fit,t(cbind(U %*% diag(sqrt(d)),X,W))),
                        FF = with(fit,t(cbind(V %*% diag(sqrt(d)),B,Z))),
                        const = sum(lfactorial(Y)))
  expect_equal(fit$loglik,tail(fit$progress$loglik,n = 1))
  expect_equal(fit$loglik,loglik)
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
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,
                                 control = list(maxiter = 20,
                                                calc_deriv = TRUE,
                                                calc_max_diff = TRUE))))
  suppressWarnings(capture.output(
    fit <- fit_glmpca_pois(Y,fit0 = fit_quick,
                           control = list(maxiter = 500,tol = 1e-8,
                                          calc_deriv = TRUE,
                                          calc_max_diff = TRUE))))
  capture.output(print(summary(fit)))

  # All the updates should monotonically increase the likelihood.
  expect_nondecreasing(fit$progress$loglik)

  # The fixed row and column intercepts should not change.
  expect_equal(fit0$X,fit$X)
  expect_equal(fit0$B,fit$B)
  expect_equal(fit0$Z,fit$Z)

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
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,
                                 control = list(maxiter = 20,
                                                calc_deriv = TRUE,
                                                calc_max_diff = TRUE))))
  suppressWarnings(capture.output(
    fit <- fit_glmpca_pois(Y,fit0 = fit_quick,
                           control = list(maxiter = 500,tol = 1e-8,
                                          calc_deriv = TRUE,
                                          calc_max_diff = TRUE))))
  capture.output(print(summary(fit)))

  # All the updates should monotonically increase the likelihood.
  expect_nondecreasing(fit$progress$loglik)

  # The fixed row and column intercepts should not change.
  expect_equal(fit0$X,fit$X)
  expect_equal(fit0$B,fit$B)
  expect_equal(fit0$Z,fit$Z)

  # Check that orthogonality constraints are satisfied.
  expect_equivalent(crossprod(fit$U),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit$V),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(fit$d,sort(fit$d,decreasing = TRUE))
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
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,
                                 control = list(maxiter = 20,
                                                calc_deriv = TRUE,
                                                calc_max_diff = TRUE))))

  expect_nondecreasing(fit_quick$progress$loglik)

  # The fixed row and column intercepts should not change.
  expect_equal(fit0$X,fit_quick$X)
  expect_equal(fit0$B,fit_quick$B)
  expect_equal(fit0$Z,fit_quick$Z)

  # Check that orthogonality constraints are satisfied.
  expect_equivalent(crossprod(fit_quick$U),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit_quick$V),diag(3),scale = 1,tolerance = 1e-8)
})

test_that("Final fit is the same with sparse and dense Y",{
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y
  Y_sp <- as(Y,"sparseMatrix")

  # Fit a GLM-PCA model to the data with dense Y.
  set.seed(1)
  fit0 <- init_glmpca_pois(Y,K = 3)
  suppressWarnings(capture.output(
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,
                                 control = list(maxiter = 20,
                                                calc_deriv = TRUE,
                                                calc_max_diff = TRUE))))

  # Fit a GLM-PCA model to the data with sparse Y.
  set.seed(1)
  fit0_sp <- init_glmpca_pois(Y_sp,K = 3)
  suppressWarnings(capture.output(
    fit_quick_sp <- fit_glmpca_pois(Y_sp,fit0 = fit0_sp,
                                    control = list(maxiter = 20,
                                                   calc_deriv = TRUE,
                                                   calc_max_diff = TRUE))))

  fit_quick$progress[,"time"]    <- 0
  fit_quick_sp$progress[,"time"] <- 0
  expect_equal(fit_quick,fit_quick_sp)
})

test_that("Final fit is the same with single thread or multiple threads",{
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y

  # Fit a GLM-PCA model to the data using 1 thread.
  set.seed(1)
  set_fastglmpca_threads(1)
  fit0 <- init_glmpca_pois(Y,K = 3)
  suppressWarnings(capture.output(
    fit1 <- fit_glmpca_pois(Y,fit0 = fit0,
                            control = list(maxiter = 20,
                                           calc_deriv = TRUE,
                                           calc_max_diff = TRUE))))

  # Fit a GLM-PCA model to the data using 2 threads.
  set.seed(1)
  set_fastglmpca_threads(2)
  fit0 <- init_glmpca_pois(Y,K = 3)
  suppressWarnings(capture.output(
    fit2 <- fit_glmpca_pois(Y,fit0 = fit0,
                            control = list(maxiter = 20,
                                           calc_deriv = TRUE,
                                           calc_max_diff = TRUE))))
  fit1$progress[,"time"] <- 0
  fit2$progress[,"time"] <- 0
  expect_equal(fit1,fit2)
})

test_that("Final fit is (roughly) the same with or without daarem",{
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y

  # Fit a GLM-PCA model to the data without daarem.
  fit0 <- init_glmpca_pois(Y,K = 3)
  suppressWarnings(capture.output(
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,
                                 control = list(use_daarem = TRUE,
                                                maxiter = 40,tol = 0.001,
                                                orthonormalize = FALSE,
                                                calc_deriv = TRUE,
                                                calc_max_diff = TRUE))))
  suppressWarnings(capture.output(
    fit1 <- fit_glmpca_pois(Y,fit0 = fit_quick,
                            control = list(use_daarem = FALSE,
                                           orthonormalize = TRUE,
                                           maxiter = 100,tol = 1e-8,
                                           calc_deriv = TRUE,
                                           calc_max_diff = TRUE))))

  # Fit a GLM-PCA model to the data with daarem.
  suppressWarnings(capture.output(
    fit2 <- fit_glmpca_pois(Y,fit0 = fit_quick,
                            control = list(use_daarem = TRUE,
                                           orthonormalize = FALSE,
                                           maxiter = 200,tol = 1e-8,
                                           calc_deriv = TRUE,
                                           calc_max_diff = TRUE))))
  fit1["progress"] <- NULL
  fit2["progress"] <- NULL
  expect_equal(fit1,fit2,scale = 1,tolerance = 0.001)

  # The fixed row and column intercepts should not change.
  expect_equal(fit0$X,fit1$X)
  expect_equal(fit0$B,fit1$B)
  expect_equal(fit0$Z,fit1$Z)
  expect_equal(fit0$X,fit2$X)
  expect_equal(fit0$B,fit2$B)
  expect_equal(fit0$Z,fit2$Z)
})

test_that("Test fit works with input covariates",{

  # Simulate a 100 x 200 data set to factorize.
  set.seed(1)
  n <- 100
  m <- 200
  Y <- generate_glmpca_data_pois(n,m,K = 3)$Y
  X <- matrix(data = rnorm(2*n),n,2)
  Z <- matrix(data = rnorm(2*m),m,2)

  # Fit a GLM-PCA model to the data and check the outputs.
  fit0 <- init_glmpca_pois(Y,K = 3, X = X, Z = Z, fixed_b_cols = c(1))
  suppressWarnings(capture.output(
    fit_quick <- fit_glmpca_pois(Y,fit0 = fit0,
                                 control = list(maxiter = 20,
                                                calc_deriv = TRUE,
                                                calc_max_diff = TRUE))))

  expect_nondecreasing(fit_quick$progress$loglik)
  expect_equal(fit0$X,fit_quick$X)
  expect_equal(fit0$Z,fit_quick$Z)
  expect_equal(fit0$B[,fit0$fixed_b_cols],fit_quick$B[,fit0$fixed_b_cols])
  expect_equal(fit0$W[,fit0$fixed_w_cols],fit_quick$W[,fit0$fixed_w_cols])

  # Check that the orthogonality constraints are satisfied.
  expect_equivalent(crossprod(fit_quick$U),diag(3),scale = 1,tolerance = 1e-8)
  expect_equivalent(crossprod(fit_quick$V),diag(3),scale = 1,tolerance = 1e-8)
})
