# ultimately I should add more things to the output and more printing
# but for now I'm just trying to make sure that everything works correctly
# then I can get fancier
# and also add to the documentation
set.seed(1)
data <- plash::generate_glmpca_data(n = 2500, p = 1250, K = 5, link = "log")
#data <- plash:::generate_data_simple(n = 1000, p = 500, K = 1, link = "log1p")

# data <- fastTopics::simulate_poisson_gene_data(n = 2500, m = 1250, k = 5)
# data$Y <- data$X

data$Y <- data$Y[rowSums(data$Y) > 0,]
data$Y <- data$Y[, colSums(data$Y) > 0]

sp_Y <- as(data$Y, "sparseMatrix")

set.seed(1)
fit0 <- plash::init_glmpca(
  Y = sp_Y, K = 5, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
)

tictoc::tic()
fast_glmpca_fit_init <- plash:::fit_glmpca(
  Y = sp_Y, fit0 = fit0, tol = 1e-4, algorithm = "ccd", link = "log",
  control = list(line_search = TRUE, num_iter = 5), max_iter = 5
)
tictoc::toc()

tictoc::tic()
fast_glmpca_fit_daarem <- plash:::fit_glmpca(
  Y = data$Y, fit0 = fast_glmpca_fit_init, tol = 1e-4, algorithm = "ccd", link = "log",
  control = list(line_search = TRUE, num_iter = 5), max_iter = 10, use_daarem = TRUE
)
tictoc::toc()

# first, get glmpca fit
# tictoc::tic()
glmpca_fit <- glmpca::glmpca(
  Y = sp_Y, L = 5, optimizer = "avagrad", ctl = list(verbose = TRUE, maxIter = 100), minibatch = "stochastic"
)
# tictoc::toc()
# 
# fit0 <- plash::init_glmpca(
#   Y = data$Y, K = 1, fit_col_size_factor = FALSE, fit_row_intercept = FALSE
# )

#fit01 <- plash::init_glmpca(
#  Y = data$Y, K = 3, fit_col_size_factor = FALSE, fit_row_intercept = FALSE
#)

set.seed(6)
fit0 <- plash::init_glmpca(
  Y = sp_Y, K = 5, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
)

#fast <- plash:::fast_fit(
#  Y = data$Y, fit0 = fit0, n_iter = 10
#)

#set.seed(6)
#fit0 <- plash::init_glmpca(
#  Y = data$Y, K = 10
#)

tictoc::tic()
fast_glmpca_fit_log1p <- plash:::fit_glmpca(
  Y = data$Y, fit0 = fit0, tol = 1e-4, algorithm = "ccd", link = "log",
  control = list(line_search = TRUE, num_iter = 5), warmup = TRUE, warmup_steps = 5, max_iter = 10
)
tictoc::toc()

set.seed(6)

fit0 <- plash::init_glmpca(
  Y = data$Y, K = 5, fit_col_size_factor = FALSE, fit_row_intercept = FALSE
)

tictoc::tic()
fast_glmpca_fit_log <- plash:::fit_glmpca(
  Y = sp_Y, fit0 = fit0, tol = 1e-4, algorithm = "ccd", link = "log",
  control = list(line_search = TRUE, num_iter = 5)
)
tictoc::toc()

# tictoc::tic()
# fast_glmpca_fit_irls <- plash:::fit_glmpca(
#   Y = data$Y, fit0 = fit01, tol = 1e-4, algorithm = "irls",
#   control = list(num_iter = 3)
# )
# tictoc::toc()

glmpca_fitted_vals <- predict(glmpca_fit)
glmpca_lik <- sum(dpois(drop(data$Y), drop(glmpca_fitted_vals), log = TRUE))
fast_glmpca_fitted_vals <- exp(t(fast_glmpca_fit$LL) %*% fast_glmpca_fit$FF)
fast_glmpca_lik <- sum(dpois(drop(data$Y), drop(fast_glmpca_fitted_vals), log = TRUE))

# The code I have now is not perfect
# I am interested in developing an IRLS algorithm
# I am also interested in printing more verbose output if neccessary

# I will worry about the IRLS algorithm first
# It would be nice if at least I could have a comparison between these few algorithms

