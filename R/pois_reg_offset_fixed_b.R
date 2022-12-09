# functions for poisson regression with fixed parameters and modified link function
pois_reg_offset_objective_fixed_b <- function(b, fixed_b, X_T, y, c) {
  
  cat_b <- c(fixed_b, b)
  exp_eta <- exp(crossprod(X_T, cat_b))
  # fix here to deal with numerical issues at initialization
  exp_eta <- pmax(exp_eta, .Machine$double.eps)
  obj <- mean(-y * log(exp_eta - c) + exp_eta)
  return(obj)
  
}

pois_reg_offset_gradient_fixed_b <- function(b, fixed_b, X_T, y, c) {
  
  fixed_p <- length(fixed_b)
  cat_b <- c(fixed_b, b)
  p <- length(cat_b)
  exp_eta <- exp(crossprod(X_T, cat_b))
  grad_cat <- drop(tcrossprod(t(((-y / (exp_eta - c)) + 1) * exp_eta), X_T))
  grad <- grad_cat[(fixed_p + 1):p]
  return(grad / length(y))
  
}

solve_pois_reg_offset_fixed_b <- function(X_T, X, y, c, fixed_b = NULL, b_init = NULL) {
  
  if (all(c <= 0)) {
    
    sol <- optim(
      par = b_init,
      fn = pois_reg_offset_objective_fixed_b,
      gr = pois_reg_offset_gradient_fixed_b,
      X_T = X_T,
      y = y,
      c = c,
      fixed_b = fixed_b,
      method = "BFGS",
      control = list(trace = 0) # verbose for now
    )
    
  } else {
    
    is_constrained <- (c > 0)
    
    fixed_p <- length(fixed_b)
    
    if (fixed_p > 0) {
      
      ci_fixed <- tcrossprod(X[is_constrained, 1:fixed_p], t(fixed_b))
      
    } else {
      
      ci_fixed <- 0
      
    }
    
    if (ncol(X) - (fixed_p + 1) == 0) {
      
      ui <- matrix(data = X[is_constrained, (fixed_p + 1):ncol(X)], ncol = 1)
      
    } else {
      
      ui <- X[is_constrained, (fixed_p + 1):ncol(X)]
      
    }
    
    sol <- constrOptim(
      theta = b_init,
      f = pois_reg_offset_objective_fixed_b,
      grad = pois_reg_offset_gradient_fixed_b,
      ui = ui,
      ci = log(.Machine$double.eps + c[is_constrained]) - ci_fixed,
      X_T = X_T,
      y = y,
      c = c,
      fixed_b = fixed_b,
      method = "BFGS",
      control = list(trace = 0) # verbose for now
    )
    
    # tryCatch(
    #   {
    #     
    #     sol <- constrOptim(
    #       theta = b_init,
    #       f = pois_reg_offset_objective_b,
    #       grad = pois_reg_offset_gradient_b,
    #       ui = X[is_constrained, ],
    #       ci = log(.Machine$double.eps + c[is_constrained]),
    #       X_T = X_T,
    #       y = y,
    #       c = c,
    #       method = "BFGS",
    #       control = list(trace = 0) # verbose for now
    #     )
    #     
    #   }, 
    #   error = function(e) {
    #     browser()
    #   }
    # )
    
  }
  
  if (sol$convergence != 0) {
    
    warning("Optimization did not converge")
    
  }
  
  return(c(fixed_b, sol$par))
  
}

