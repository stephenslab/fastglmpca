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
  const <- 1e-5 + cc - exp_eta
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


solve_pois_reg_unknown_offset <- function(X, y, init) {
  
  opts <- list("algorithm"="NLOPT_LD_SLSQP",
               "xtol_rel"=1.0e-6)
  
  sol <- nloptr::nloptr(
    x0 = init,
    eval_f = pois_reg_unknown_offset_obj_grad,
    eval_g_ineq = pois_reg_unknown_offset_const_jac,
    opts = opts,
    X = X, 
    y = y
  )
  
  # first element is c, the rest is b
  return(sol$solution)
  
}
