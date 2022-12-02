# Here, code to solve poisson regression with unknown offset
pois_reg_unknown_offset_obj_grad <- function(params, X, y) {
  
  n <- length(y)
  p <- length(params)
  cc <- params[1]
  b <- params[2:p]
  
  exp_eta <- exp(X %*% b)
  obj <- mean(-y * log(exp_eta - cc) + exp_eta) - cc
  
  # now, the gradient 
  grad_b <- (drop(t(((-y / (exp_eta - cc)) + 1) * exp_eta) %*% X)) / n
  grad_c <- mean(y / (exp_eta - cc)) - 1
  grad <- c(grad_c, grad_b)
  
  return(
    list(
      "objective" = obj,
      "gradient" = grad
    )
  )
  
}

pois_reg_unknown_offset_const_jac <- function(params, X, y) {
  
  n <- length(y)
  p <- length(params)
  cc <- params[1]
  b <- params[2:p]
  
  exp_eta <- exp(X %*% b)
  # constraint on the regression
  const <- .Machine$double.eps + cc - exp_eta
  # constraint on c
  const <- c(-cc, const)
  
  jac <- -X * matrix(data = rep(exp_eta, (p - 1)), nrow = n, ncol = p - 1)
  jac <- rbind(rep(0, p - 1), jac)
  jac <- cbind(c(-1, rep(1, n)), jac)
  
  return(
    list(
      "constraints" = const,
      "jacobian" = jac
    )
  )
  
}


solve_pois_reg_unknown_offset <- function(X, y) {
  
  opts <- list("algorithm"="NLOPT_LD_MMA",
               "xtol_rel"=1.0e-8)
  
  p <- ncol(X)
  init <- c(0, rnorm(p))
  
  sol <- nloptr::nloptr(
    x0 = init,
    eval_f = pois_reg_unknown_offset_obj_grad,
    eval_g_ineq = pois_reg_unknown_offset_const_jac,
    opts = opts,
    X = X, 
    y = y
  )
  
  return(
    list(
      cc = sol$solution[1],
      b = sol$solution[2:(p+1)]
    )
  )
  
}

generate_pois_reg_data <- function(n, p) {
  
  # first, will generate without the offset
  X <- matrix(
    data = runif(n * p), nrow = n, ncol = p
  )
  
  b <- rnorm(p)
  
  exp_eta <- exp(X %*% b)
  
  y <- rpois(
    n = n, lambda = exp_eta
  )
  
  return(
    list(
      y = y,
      X = X,
      b = b
    )
  )
  
}

n <- 1000
p <- 5

data <- generate_pois_reg_data(n, p)

sol <- solve_pois_reg_unknown_offset(data$X, data$y)
