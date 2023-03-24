pois_reg_logsp_objective_b <- function(b, X_T, y, s) {
  
  exp_eta <- exp(crossprod(X_T, b))
  exp_eta_minus_one <- pmax(exp_eta - 1, .Machine$double.eps)
  obj <- mean(-y * log(exp_eta_minus_one) + s * exp_eta)
  if (is.na(obj)){
    
    print(0)
    
  }
  return(obj)
  
}

pois_reg_logsp_gradient_b <- function(b, X_T, y, s) {
  
  exp_eta <- exp(crossprod(X_T, b))
  exp_eta_minus_one <- pmax(exp_eta - 1, .Machine$double.eps)
  grad <- as.vector(tcrossprod(t((s - (y / (exp_eta_minus_one))) * exp_eta), X_T))
  return(grad / length(y))
  
}

solve_pois_reg_logsp_b <- function(X_T, y, s, b_init = NULL) {
  
  if (is.null(b_init)) {
    
    p <- nrow(X_T)
    b_init <- rep(0, p)
    
  }
    
  sol <- optim(
    par = b_init,
    fn = pois_reg_logsp_objective_b,
    gr = pois_reg_logsp_gradient_b,
    X_T = X_T,
    y = y,
    s = s,
    method = "L-BFGS-B",
    lower = rep(.Machine$double.eps, length(b_init)),
    control = list(trace = 0, maxit = 5) # verbose for now
  )

  return(sol$par)
  
}

