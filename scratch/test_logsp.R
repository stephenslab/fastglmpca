# I want to use this script in order to test the logsp link function
# first I just want to make sure that it actually converges
# then, I can make sure that it will converge to the same solution when s = 1

set.seed(1)
data <- plash::generate_glmpca_data(n = 1000, p = 500, K = 1, link = "log1p")

data$Y <- data$Y[rowSums(data$Y) > 0,]
data$Y <- data$Y[, colSums(data$Y) > 0]

set.seed(6)
fit0_1 <- plash::init_glmpca(
  Y = data$Y, K = 1, fit_col_size_factor = FALSE, fit_row_intercept = FALSE
)

tictoc::tic()
fit_1 <- plash:::fit_glmpca(
  Y = data$Y, fit0 = fit0_1, tol = 1e-4, algorithm = "ccd", link = "log1p",
  control = list(line_search = TRUE, num_iter = 5)
)
tictoc::toc()

k <- 1
s <- 1
n <- nrow(data$Y)
p <- ncol(data$Y)
set.seed(6)
init_LL <- matrix(
  data = runif(n = n * (k + 1), max = .1), ncol = n, nrow = k + 1
)

init_LL[1, ] <- 1

init_FF <- matrix(
  data = runif(n =p * (k+1), max = .1), ncol = p, nrow = k + 1
)

init_FF[1, ] <- log(s)

fit0_2 <- plash::init_glmpca(
  Y = data$Y, LL = init_LL, FF = init_FF, fixed_loadings = c(1), fixed_factors = c(1)
)

fit_2 <- plash::fit_glmpca(
  Y = data$Y, fit0 = fit0_2, tol = 1e-4, algorithm = "ccd", link = "logsp",
  control = list(line_search = TRUE, num_iter = 5)
)

fit_1_fitted <- exp(crossprod(fit_1$LL, fit_1$FF)) - 1
fit_2_fitted <- exp(crossprod(fit_1$LL, fit_1$FF)) - s

################# TESTING EQUIVALENCE CONJECTURE #################

# first, I will test to see if my conjecture about NMF is correct
# this should be pretty easy using fastTopics

#first, generate data from an NMF model
set.seed(1)
data <- fastTopics::simulate_poisson_gene_data(n = 2000, m = 1000, k = 2)
data$Y <- data$X

data$Y <- data$Y[rowSums(data$Y) > 0,]
data$Y <- data$Y[, colSums(data$Y) > 0]

#nmf_fit <- fastTopics:::fit_pnmf_rank1(data$Y)
nmf_fit <- fastTopics::fit_poisson_nmf(X = as(data$Y, "sparseMatrix"), k = 2)

# set.seed(1)
# data <- plash::generate_glmpca_data(n = 4000, p = 2000, K = 1, link = "log1p")
# 
# data$Y <- data$Y[rowSums(data$Y) > 0,]
# data$Y <- data$Y[, colSums(data$Y) > 0]

k <- 2
s <- 1e4
n <- nrow(data$Y)
p <- ncol(data$Y)
set.seed(6)
init_LL <- matrix(
  data = runif(n = n * (k + 1), max = 1), ncol = n, nrow = k + 1
)

init_LL[1, ] <- 1

init_FF <- matrix(
  data = runif(n =p * (k+1), max = 1), ncol = p, nrow = k + 1
)

init_FF[1, ] <- log(s)

fit0_2 <- plash::init_glmpca(
  Y = data$Y, LL = init_LL, FF = init_FF, fixed_loadings = c(1), fixed_factors = c(1)
)

fit_2 <- plash::fit_glmpca(
  Y = data$Y, fit0 = fit0_2, tol = 1e-4, algorithm = "ccd", link = "logsp", s = s,
  control = list(line_search = TRUE, num_iter = 2)
)

# normal_fit0 <- plash::init_glmpca(
#   Y = data$Y, K = 1, fit_col_size_factor = FALSE, fit_row_intercept = FALSE
# )
# 
# normal_fit <- plash::fit_glmpca(
#   Y = data$Y, fit0 = normal_fit0,
# )

expected_H <- tcrossprod(nmf_fit$L, nmf_fit$F) / s
actual_H <- crossprod(fit_2$LL, fit_2$FF)

# I might be able to test the theorem from the overleaf just through 
get_nmf_lik <- function(Y, H) {
  
  sum(
    dpois(x = as.vector(Y), lambda = as.vector(H), log = TRUE)
  )
  
}

nmf_lik <- get_nmf_lik(data$Y, tcrossprod(nmf_fit$L, nmf_fit$F))

get_logsp_lik <- function(Y, H, s) {
  
  Lambda <- s * (exp(H) - 1)
  sum(
    dpois(x = as.vector(Y), lambda = as.vector(Lambda), log = TRUE)
  )
  
}

logsp_lik <- get_logsp_lik(data$Y, tcrossprod(nmf_fit$L, nmf_fit$F) / (1e6), 1e6)

get_glmpca_lik <- function(Y, H) {
  
  Lambda <- exp(H)
  sum(
    dpois(x = as.vector(Y), lambda = as.vector(Lambda), log = TRUE)
  )
  
}

data2 <- plash:::generate_glmpca_data(n = 3000, p = 1500, K = 1)

glmpca_fit <- plash::fit_glmpca(Y = data2$Y, K = 1)

glmpca_lik <- get_glmpca_lik(data$Y, crossprod(glmpca_fit$LL, glmpca_fit$FF))

logsp_lik <- get_logsp_lik(data$Y, crossprod(glmpca_fit$LL, glmpca_fit$FF) - log(1e-6), 1e-6)




