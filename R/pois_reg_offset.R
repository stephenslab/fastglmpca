pois_reg_offset_objective_b <- function(b, X, y, c) {

  exp_eta <- exp(X %*% b)
  obj <- mean(-y * log(exp_eta - c) + exp_eta)
  return(obj)

}

pois_reg_offset_gradient_b <- function(b, X, y, c) {

  exp_eta <- exp(X %*% b)
  grad <- drop(t(((-y / (exp_eta - c)) + 1) * exp_eta) %*% X)
  return(grad / length(y))

}

solve_pois_reg_offset_b <- function(X, y, c, b_init = NULL) {

  if (is.null(b_init)) {

    p <- ncol(X)
    b_init <- rep(0, p)

  }

  if (all(c <= 0)) {

    sol <- optim(
      par = b_init,
      fn = pois_reg_offset_objective_b,
      gr = pois_reg_offset_gradient_b,
      X = X,
      y = y,
      c = c,
      method = "BFGS",
      control = list(trace = 0) # verbose for now
    )

  } else {

    is_constrained <- (c > 0)

    # TODO: ensure feasible initialization
    sol <- constrOptim(
      theta = b_init,
      f = pois_reg_offset_objective_b,
      grad = pois_reg_offset_gradient_b,
      ui = X[is_constrained, ],
      ci = log(.Machine$double.eps + c[is_constrained]),
      X = X,
      y = y,
      c = c,
      method = "BFGS",
      control = list(trace = 0) # verbose for now
    )

  }

  if (sol$convergence != 0) {

    warning("Optimization did not converge")

  }

  return(sol$par)

}

pois_reg_offset_objective_c <- function(c, exp_eta, y) {

  obj <- mean(-y * log(exp_eta - c) - c)
  return(obj)

}

pois_reg_offset_gradient_c <- function(c, exp_eta, y) {

  grad <- mean(y / (exp_eta - c)) - 1
  return(grad)

}

solve_pois_reg_offset_c <- function(X, y, b, c_init = NULL) {

  if (is.null(c_init)) {

    c_init <- 0

  }

  exp_eta <- exp(X %*% b)
  sol <- optim(
    par = c_init,
    fn = pois_reg_offset_objective_c,
    gr = pois_reg_offset_gradient_c,
    method = "L-BFGS-B",
    upper =  min(exp_eta) - .Machine$double.neg.eps,
    exp_eta = exp_eta,
    y = y
  )
  return(sol$par)

}

