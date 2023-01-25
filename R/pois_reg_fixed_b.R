# Here, write code to run glmpca with c = 0
# The objective and derivative should be much less complex
pois_reg_objective_fixed_b <- function(b, fixed_b, X_T, y) {
  
  cat_b <- c(fixed_b, b)
  eta <- crossprod(X_T, cat_b)
  obj <- sum(-y * eta + exp(eta))
  return(obj)
  
}

pois_reg_gradient_fixed_b <- function(b, fixed_b, X_T, y) {
  
  fixed_p <- length(fixed_b)
  cat_b <- c(fixed_b, b)
  p <- length(cat_b)
  exp_eta <- exp(crossprod(X_T, cat_b))
  grad_cat <- tcrossprod(t(exp_eta - y), X_T)
  grad <- grad_cat[(fixed_p + 1):p]
  return(grad)
  
}

solve_pois_reg_fixed_b <- function(X_T, X, y, fixed_b = NULL, b_init = NULL) {
    
  sol <- optim(
    par = b_init,
    fn = pois_reg_objective_fixed_b,
    gr = pois_reg_gradient_fixed_b,
    X_T = X_T,
    y = y,
    fixed_b = fixed_b,
    method = "BFGS",
    control = list(trace = 0) # verbose for now
  )

  if (sol$convergence != 0) {
    
    warning("Optimization did not converge")
    
  }
  
  return(c(fixed_b, sol$par))
  
}
