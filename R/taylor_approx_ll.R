approx_reg_ccd <- function(
  X,
  y,
  w,
  a1,
  a2,
  update_indices,
  n_updates,
  alpha,
  beta
) {
  
  zero_idx <- which(y == 0)
  non_zero_idx <- which(y != 0)
  num_zero <- length(zero_idx)
  y_nz <- y[non_zero_idx]
  X_0 <- X[zero_idx, ]
  X_nz <- X[non_zero_idx, ]
  Z_0 <- crossprod(X_0)
  m_0 <- crossprod(X_0[, update_indices], rep(-a1, num_zero))
  m_nz <- crossprod(X_nz[, update_indices], y_nz)
  m <- m_0 + m_nz
  sqrd_factor <- 2 * a2
  
  p <- length(w)
  
  eta <- X_nz %*% w
  exp_eta <- exp(eta)
  
  current_nonlinear_lik <- a2 * tcrossprod(w %*% Z_0, w) + sum(exp_eta)
  
  for (updates in 1:n_updates) {
    
    for (idx in 1:length(update_indices)) {
      
      j = update_indices[idx]

      current_lik <- current_nonlinear_lik - w[j] * m[idx]
      deriv_prod <- exp_eta * X_nz[, j]
      first_deriv <- sqrd_factor * sum(w * Z_0[,j]) + sum(deriv_prod) - m[idx]
      second_deriv <- sqrd_factor * Z_0[j, j] + sum(deriv_prod * X_nz[, j])
      
      newton_dir = first_deriv / second_deriv
      
      t = 1.0
      step_accepted = FALSE
      newton_dec = alpha * first_deriv * newton_dir
      w_j_og = w[j]
      eta_og = eta
      exp_eta_og = exp_eta
      
      while(!step_accepted) {
        
        w[j] = w_j_og - t * newton_dir
        eta_proposed = eta + (w[j] - w_j_og) * X_nz[, j]
        exp_eta_proposed = exp(eta_proposed)
        f_proposed = a2 * tcrossprod(w %*% Z_0, w) + sum(exp_eta_proposed) - w[j] * m[idx]
        
        if (f_proposed <= current_lik - t * newton_dec) {
          
          step_accepted = TRUE
          current_nonlinear_lik = f_proposed + w[j] * m[idx]
          eta = eta_proposed
          exp_eta = exp_eta_proposed
          
        } else {
          
          t = beta * t
          if (abs((t * newton_dir) / w_j_og) < 1e-16) {
            
            w[j] = w_j_og
            step_accepted = TRUE
            eta = eta_og
            exp_eta = exp_eta_og
            
          }
          
        }
        
      }
      
    }
    
  }
  
  return(w)
  
}

