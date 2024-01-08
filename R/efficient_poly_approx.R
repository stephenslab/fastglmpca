get_b_complete <- function(b, b_fixed, b_fixed_idx) {

  b_complete <- numeric(length(b) + length(b_fixed))
  b_complete[b_fixed_idx] <- b_fixed
  b_idx <- 1:length(b_complete)
  b_variable_idx <- b_idx[!(b_idx %in% b_fixed_idx)]
  b_complete[b_variable_idx] <- b
  return(b_complete)

}

get_approx_loglik <- function(
    b,
    b_fixed,
    b_fixed_idx,
    y,
    X,
    X_cross_y_a1,
    X_0_cross_X,
    a2,
    y_nonzero_idx
) {

  b_complete <- get_b_complete(b, b_fixed, b_fixed_idx)

  loglik <- crossprod(b, X_cross_y_a1) -
    a2 * (crossprod(b_complete, X_0_cross_X) %*% b_complete) -
    sum(exp(X[y_nonzero_idx, , drop = FALSE] %*% b_complete))

  return(loglik)

}

get_deriv_loglik <- function(
    b,
    b_fixed,
    b_fixed_idx,
    y,
    X,
    X_cross_y_a1,
    X_0_cross_X,
    a2,
    y_nonzero_idx
) {

  b_complete <- get_b_complete(b, b_fixed, b_fixed_idx)

  grad <- X_cross_y_a1 -
    2 * a2 * (X_0_cross_X[-b_fixed_idx, , drop = FALSE] %*% b_complete) -
    crossprod(
      X[y_nonzero_idx, -b_fixed_idx, drop = FALSE],
      exp(X[y_nonzero_idx, , drop = FALSE] %*% b_complete)
    )

  return(grad)

}

get_opt_b_approx <- function(
  b_init,
  b_fixed_idx,
  y,
  X,
  a1,
  a2
) {

  y_zero_idx <- which(y == 0)
  y_nonzero_idx <- which(y != 0)

  n <- length(y)
  a1_vec <- rep(a1, n)
  a1_vec[y_nonzero_idx] <- 0

  if (length(y_zero_idx) == 0) {

    X_0_cross_X <- matrix(
      data = 0, nrow = length(b_init), ncol = length(b_init)
    )

  } else {

    X_0_cross_X <- crossprod(X[y_zero_idx, , drop = FALSE])

  }

  opt_out <- optim(
    par = b_init[-b_fixed_idx],
    fn = get_approx_loglik,
    gr = get_deriv_loglik,
    method = "BFGS",
    control = list(fnscale=-1), # maximize
    b_fixed = b_init[b_fixed_idx],
    b_fixed_idx = b_fixed_idx,
    y = y,
    X = X,
    X_cross_y_a1 = crossprod(X[, -b_fixed_idx], y - a1_vec),
    X_0_cross_X = X_0_cross_X,
    a2 = a2,
    y_nonzero_idx = y_nonzero_idx
  )

  b_opt <- get_b_complete(opt_out$par, b_init[b_fixed_idx], b_fixed_idx)
  return(b_opt)

}

# set.seed(1)
# n <- 1000
# x1 <- abs(rnorm(n))
# x2 <- abs(rnorm(n))
# X <- matrix(
#   data = c(rep(1, n), x1, x2),
#   ncol = 3
# )
# b <- c(-.1, -0.5, -0.3)
# r <- exp(X %*% b)
# y <- rpois(n,r)
#
# ans <- pracma::polyApprox(exp,a = -3.5,b = 0.25,n = 2)
# a1  <- ans$p[2]
# a2  <- ans$p[1]
#
# optim_val <- get_opt_b_approx(
#   b_init = c(-.1, 0, 0),
#   b_fixed_idx = c(1),
#   y = y,
#   X = X,
#   a1 = a1,
#   a2 = a2
# )

#g <- glm.fit(x = X, y = y, family = poisson())


# get_approx_loglik <- function(
#     y, X, b, b_fixed, b_fixed_idx, a1, a2, y_zero_idx, y_nonzero_idx
#   ) {
#
#   n <- length(y)
#   b_total <- c(b_fixed, b)
#   a1_vec <- rep(a1, n)
#   a1_vec[y_nonzero_idx] <- 0
#
#   loglik <- t(b) %*% crossprod(X[, -b_fixed_idx], y - a1_vec) -
#     a2 * t(b_total) %*% crossprod(X[y_zero_idx, ]) %*% b_total -
#     sum(exp(X[y_nonzero_idx, ] %*% b_total))
#
#   return(loglik)
#
# }


