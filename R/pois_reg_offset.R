pois_reg_offset_objective_b <- function(b, X_T, y, c) {

  exp_eta <- exp(crossprod(X_T, b))
  obj <- mean(-y * log(exp_eta - c) + exp_eta)
  return(obj)

}

pois_reg_offset_gradient_b <- function(b, X_T, y, c) {

  exp_eta <- exp(crossprod(X_T, b))
  grad <- drop(tcrossprod(t(((-y / (exp_eta - c)) + 1) * exp_eta), X_T))
  return(grad / length(y))

}

solve_pois_reg_offset_b <- function(X_T, X, y, c, b_init = NULL) {

  if (is.null(b_init)) {

    p <- nrow(X_T)
    b_init <- rep(0, p)

  }

  if (all(c <= 0)) {

    sol <- optim(
      par = b_init,
      fn = pois_reg_offset_objective_b,
      gr = pois_reg_offset_gradient_b,
      X_T = X_T,
      y = y,
      c = c,
      method = "BFGS",
      control = list(trace = 0) # verbose for now
    )

  } else {

    is_constrained <- (c > 0)
    
    sol <- constrOptim(
      theta = b_init,
      f = pois_reg_offset_objective_b,
      grad = pois_reg_offset_gradient_b,
      ui = X[is_constrained, ],
      ci = log(1e-10 + c[is_constrained]),
      X_T = X_T,
      y = y,
      c = c,
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

  return(sol$par)

}

pois_reg_offset_objective_c <- function(c, exp_eta, y, size) {

  if (any(exp_eta - size * c <= 0)) {
    
    print(0)
    
  }
  obj <- mean(-y * log(exp_eta - size * c) - size * c)
  return(obj)

}

pois_reg_offset_gradient_c <- function(c, exp_eta, y, size) {

  grad <- mean(y * size / (exp_eta - c) - size)
  return(grad)

}

solve_pois_reg_offset_c <- function(X_T, y, b, size, c_init = NULL) {

  if (is.null(c_init)) {

    c_init <- 0

  }

  exp_eta <- exp(crossprod(X_T, b))
  sol <- optim(
    par = c_init,
    fn = pois_reg_offset_objective_c,
    gr = pois_reg_offset_gradient_c,
    method = "L-BFGS-B",
    upper =  min(exp_eta / size) - 5e-10,
    lower = 0,
    exp_eta = exp_eta,
    y = y,
    size = size
  )
  
  if (any((exp_eta - sol$par) == 0)) {
    
    browser()
    
  }
  
  return(sol$par)

}

