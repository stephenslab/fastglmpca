set.seed(1)
data <- plash::generate_glmpca_data(n = 2500, p = 1250, K = 3)

data$Y <- data$Y[rowSums(data$Y) > 0,]

sp_Y <- as(data$Y, "sparseMatrix")

tictoc::tic()
glmpca_fit <- glmpca::glmpca(
  Y = data$Y, L = 3, optimizer = "fisher", ctl = list(verbose = TRUE)
)
tictoc::toc()

fit0 <- plash::init_glmpca(
  Y = sp_Y, K = 3, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
)

tictoc::tic()
fast_glmpca_fit <- plash::fit_glmpca(
  Y = sp_Y, fit0 = fit0, tol = 1e-4, algorithm = "ccd"
)
tictoc::toc()

glmpca_fitted_vals <- predict(glmpca_fit)
glmpca_lik <- sum(dpois(drop(data$Y), drop(glmpca_fitted_vals), log = TRUE))
fast_glmpca_fitted_vals <- fitted(fast_glmpca_fit)
fast_glmpca_lik <- sum(dpois(drop(data$Y), drop(fast_glmpca_fitted_vals), log = TRUE))

cor(as.vector(glmpca_fitted_vals), as.vector(fast_glmpca_fitted_vals))
