pois_reg_offset_objective_b <- function(b, X, y, c) {

  exp_eta <- exp(X %*% b)
  obj <- mean(-y * log(exp_eta - c) + exp_eta)
  return(obj)

}

pois_reg_offset_gradient_b <- function(b, X, y, c) {

  n <- length(y)
  exp_eta <- exp(X %*% b)
  grad <- drop(t((-y / (exp_eta - c)) * exp_eta) %*% X + t(exp_eta) %*% X)
  return(grad / n)

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
      control = list(trace = 6) # verbose for now
    )

  } else {

    is_constrained <- (y == 0) & (c > 0)

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
      control = list(trace = 6) # verbose for now
    )

  }

  if (sol$convergence != 0) {

    stop("Optimization did not converge")

  }

  return(sol$par)

}

pois_reg_offset_objective_c <- function(c, exp_eta, y) {

  obj <- mean(-y * log(exp_eta - c) - c)
  return(obj)

}

solve_pois_reg_offset_c <- function(X, y, b) {

  exp_eta <- exp(X %*% b)
  sol <- optimize(
    f = pois_reg_offset_objective_c,
    interval = c(-max(y), min(exp_eta) - .Machine$double.neg.eps),
    exp_eta = exp_eta
  )$minimum
  return(sol)

}

