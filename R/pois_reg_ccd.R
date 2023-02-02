# here, I want to solve a poisson regression with fixed coordinates
# using the ccd algorithm. I'll need a number of helper functions along the way

# j is the coordinate number of b to get the gradient wrt
pois_reg_first_deriv_fixed_b <- function(exp_eta, x, y) {
  
  deriv <- sum((exp_eta - y) * x)
  return(deriv)
  
}

pois_reg_second_deriv_fixed_b <- function(exp_eta, x) {
  
  deriv <- sum((x ^ 2) * exp_eta)
  return(deriv)
  
}

solve_pois_reg_fixed_b_ccd <- function(
    X_T, 
    y, 
    b_init, 
    max_updates = 5, 
    tol = 1e-6,
    fixed_b = NULL,
    line_search = TRUE
  ) {
  
  num_updates <- 0
  convergence <- FALSE
  p <- length(b_init)
  fixed_p <- length(fixed_b)
  b <- c(fixed_b, b_init)
  
  start_update_lik <- pois_reg_objective(b, X_T, y)
  
  while(num_updates < max_updates && !convergence) {
    
    for (j in c((fixed_p + 1):p)) {
      
      eta <- crossprod(X_T, b)
      exp_eta <- exp(eta)
      current_lik <- pois_reg_objective_fast(exp_eta, eta, y)
      
      first_deriv <- pois_reg_first_deriv_fixed_b(
        exp_eta = exp_eta, x = X_T[j, ], y = y
      )
      
      second_deriv <- pois_reg_second_deriv_fixed_b(
        exp_eta = exp_eta, x = X_T[j, ]
      )
      
      newton_direction <- (first_deriv / second_deriv)
      
      if (line_search) {
        
        alpha <- .25
        beta <- .5
        
        t <- 1
        step_accepted <- FALSE
        
        dir_deriv <- sum(-first_deriv * newton_direction)
        
        b_j_og <- b[j]
        
        while (!step_accepted) {
          
          print(glue::glue("Attempting line search with t = {t}"))
          
          b[j] <- b[j] - t * newton_direction
          f_proposed <- pois_reg_objective(b, X_T, y)
          
          if (f_proposed <= current_lik + alpha * t * dir_deriv) {
            
            step_accepted <- TRUE
            
          } else {
            
            t <- beta * t
            b[j] <- b_j_og
            if (t < 1e-10) {
              
              print(0)
              
            }
            
          }
          
        }
        
      }
      
    }
    
    # now, assess convergence 
    new_lik <- f_proposed
    
    rel_improvement <- abs((new_lik - start_update_lik) / start_update_lik)
    if (rel_improvement < tol) {
      
      convergence <- TRUE
      
    }
    
    num_updates <- num_updates + 1
    
  }
  
  return(b)
  
}

