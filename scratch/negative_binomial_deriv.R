# Here, I want to make sure that I've taken the negative binomial derivates correctly
# once I've done this then I can move on to thinking about the actual implementation
# dealing with the disperion parameter could be slightly tricky
library(numDeriv)
neg_bin_obj <- function(X, y, b, alpha) {
  
  eta <- as.vector(X %*% b)
  lik_vec <- (y + (1 / alpha)) * log(1 + alpha * exp(eta)) - y * eta
  lik <- sum(lik_vec)
  return(lik)
  
}

neg_bin_deriv <- function(X, y, b, alpha, k) {
  
  eta <- as.vector(X %*% b)
  exp_eta <- exp(eta)
  
  deriv_vec <- ((y + (1 / alpha)) / (1 + alpha * exp_eta)) * (alpha * exp_eta * X[, k]) - y * X[, k]
  deriv <- sum(deriv_vec)
  return(deriv)
  
}

neg_bin_second_deriv <- function(X, y, b, alpha, k) {
  
  eta <- as.vector(X %*% b)
  exp_eta <- exp(eta)
  
  second_deriv_vec <- ((y + (1 / alpha)) / (1 + alpha * exp_eta)) * (alpha * exp_eta * (X[, k] ^ 2)) -
    ((y + (1 / alpha)) / (1 + alpha * exp_eta) ^ 2) * ((alpha * exp_eta * X[, k]) ^ 2)
  second_deriv <- sum(second_deriv_vec)
  return(second_deriv)
  
}

num_deriv <- function(X, y, b, alpha, k) {
  
  grad(
    func = neg_bin_obj, x = b, X = X, y = y, alpha = alpha
  )[k]
  
}

num_hess <- function(X, y, b, alpha, k) {
  
  hessian(
    func = neg_bin_obj, x = b, X = X, y = y, alpha = alpha
  )[k, k]
  
}


n <- 100
p <- 5

X_dat <- runif(n = n * p)
X <- matrix(data = X_dat, nrow = n, ncol = p)
b <- runif(p)

mu <- exp(X %*% b)
y <- rnbinom(n = n, mu = mu, size = 10)

#num_deriv(X, y, b, (1/10), 3)
#neg_bin_deriv(X, y, b, (1/10), 3)

num_hess(X, y, b, (1/10), 2)
neg_bin_second_deriv(X, y, b, (1/10), 2)



