# Here, code to solve poisson regression with unknown offset
pois_reg_unknown_offset_obj_grad_fixed_b <- function(params, fixed_b, X, y, size) {
  
  n <- length(y)
  p <- length(params)
  fixed_p <- length(fixed_b)
  cc <- params[1]
  b <- params[2:p]
  cat_b <- c(fixed_b, b)
  p_cat_b <- length(cat_b)
  
  exp_eta <- exp(X %*% cat_b)
  obj <- mean(-y * log(exp_eta - size * cc) + exp_eta - size * cc) 
  
  # now, the gradient 
  grad_cat_b <- (drop(t(((-y / (exp_eta - size * cc)) + 1) * exp_eta) %*% X)) / n
  grad_b <- grad_cat_b[(fixed_p + 1):p_cat_b]
  grad_c <- mean(y / (exp_eta - size * cc) - size)
  grad <- c(grad_c, grad_b)
  
  return(
    list(
      "objective" = obj,
      "gradient" = grad
    )
  )
  
}

pois_reg_unknown_offset_const_jac_fixed_b <- function(params, fixed_b, X, y, size) {
  
  n <- length(y)
  p <- length(params)
  fixed_p <- length(fixed_b)
  cc <- params[1]
  b <- params[2:p]
  cat_b <- c(fixed_b, b)

  exp_eta <- exp(X %*% cat_b)
  # constraint on the regression
  const <- 1e-5 + size * cc - exp_eta
  # constraint on c
  const <- c(-cc, const)
  
  jac <- -X[, (fixed_p + 1):ncol(X)] * matrix(data = rep(exp_eta, (p - 1)), nrow = n, ncol = p - 1)
  jac <- rbind(rep(0, p - 1), jac)
  jac <- cbind(c(-1, size), jac)
  
  return(
    list(
      "constraints" = const,
      "jacobian" = jac
    )
  )
  
}


solve_pois_reg_unknown_offset_fixed_b <- function(X, y, size, fixed_b, init) {
  
  opts <- list("algorithm"="NLOPT_LD_SLSQP",
               "xtol_rel"=1.0e-6)
  
  sol <- nloptr::nloptr(
    x0 = init,
    eval_f = pois_reg_unknown_offset_obj_grad_fixed_b,
    eval_g_ineq = pois_reg_unknown_offset_const_jac_fixed_b,
    opts = opts,
    X = X, 
    y = y,
    fixed_b = fixed_b,
    size = size
  )
  
  # first element is c, the rest is b
  sol_out <- c(sol$solution[1], fixed_b, sol$solution[2:length(sol$solution)])
  return(sol_out)
  
}

# Now, I want to attempt to use constrOptim to solve the poisson regression
# with known offset
# Here, code to solve poisson regression with unknown offset
pois_reg_known_offset_obj_grad_fixed_b <- function(params, fixed_b, X, y, size, cc) {
  
  n <- length(y)
  fixed_p <- length(fixed_b)
  b <- params
  cat_b <- c(fixed_b, b)
  p_cat_b <- length(cat_b)
  
  exp_eta <- exp(X %*% cat_b)
  obj <- sum(-y * log(exp_eta - size * cc) + exp_eta - size * cc) 
  
  # now, the gradient 
  grad_cat_b <- (drop(t(((-y / (exp_eta - size * cc)) + 1) * exp_eta) %*% X))
  grad_b <- grad_cat_b[(fixed_p + 1):p_cat_b]
  grad <- grad_b
  
  return(
    list(
      "objective" = obj,
      "gradient" = grad
    )
  )
  
}

pois_reg_known_offset_const_jac_fixed_b <- function(params, fixed_b, X, y, size, cc) {
  
  n <- length(y)
  p <- length(params)
  fixed_p <- length(fixed_b)
  b <- params
  cat_b <- c(fixed_b, b)
  
  if (fixed_p > 0) {
    
    ci_fixed <- tcrossprod(X[ , 1:fixed_p], t(fixed_b))
    
  } else {
    
    ci_fixed <- 0
    
  }
  
  if (ncol(X) - (fixed_p + 1) == 0) {
    
    ui <- matrix(data = X[, (fixed_p + 1):ncol(X)], ncol = 1)
    
  } else {
    
    ui <- X[, (fixed_p + 1):ncol(X)]
    
  }
  
  const <- -(ui %*% b - (log(1e-5 + cc) - ci_fixed))
  jac <- -ui
  
  return(
    list(
      "constraints" = const,
      "jacobian" = jac
    )
  )
  
}


solve_pois_reg_known_offset_fixed_b_v2 <- function(X, y, size, fixed_b, init, cc) {
  
  opts <- list("algorithm"="NLOPT_LD_SLSQP",
               "xtol_rel"=1.0e-6)
  
  sol <- nloptr::nloptr(
    x0 = init,
    eval_f = pois_reg_known_offset_obj_grad_fixed_b,
    eval_g_ineq = pois_reg_known_offset_const_jac_fixed_b,
    opts = opts,
    X = X, 
    y = y,
    fixed_b = fixed_b,
    size = size, 
    cc = cc
  )
  
  # first element is c, the rest is b
  sol_out <- c(fixed_b, sol$solution)
  return(sol_out)
  
}
