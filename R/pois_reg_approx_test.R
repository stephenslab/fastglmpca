# I think it would make more sense to try to write this in R first
# and see if I can get it to work

approx_pois_reg_expected_ll <- function(
  X,
  y,
  n,
  linear_term_const,
  b,
  num_iter
) {
  
  alpha = .01
  beta = .25
  
  p <- length(b)
  
  for (update_num in 1:num_iter) {
    
    for (j in 1:p) {
      
      if (update_num == 62 && j == 5) {
        
        print(0)
        
      }
      
      sum_b_sqrd <- sum(b ^ 2)
      current_lik = n * exp(.5 * sum_b_sqrd) - b[j] * linear_term_const[j]
      first_deriv = n * exp(.5 * sum_b_sqrd) * b[j] - linear_term_const[j]
      second_deriv = n * exp(.5 * sum_b_sqrd) * (1 + b[j] * b[j])
      
      newton_dir = first_deriv / second_deriv
      
      # start line search
      step_accepted = FALSE
      newton_dec = alpha * first_deriv * newton_dir
      b_j_og = b[j]
      t = 1
      
      while(!step_accepted) {
        
        b[j] = b_j_og - t * newton_dir
        sum_b_sqrd <- sum(b ^ 2)
        f_proposed = n * exp(.5 * sum_b_sqrd) - b[j] * linear_term_const[j]
        
        if (f_proposed <= current_lik - t * newton_dec) {
          
          step_accepted = TRUE
          
        } else {
          
          t = beta * t
          
        }
        
      }
      
    }
    
  }
  
  return(b)
  
}
