#include <RcppArmadillo.h>
#include <omp.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec solve_pois_reg_cpp (
    const arma::mat X, 
    const arma::vec y, 
    arma::vec b, 
    unsigned int update_start_idx,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
  ) {
  
  double current_lik; // used to store log likelihood of each iteration
  double first_deriv;
  double second_deriv;
  double newton_dir;
  arma::vec eta;
  arma::vec exp_eta;
  double t;
  bool step_accepted;
  double b_j_og;
  double f_proposed;
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
    for (int j = update_start_idx; j < b.size(); j++) {
      
      eta = X * b;
      exp_eta = exp(eta);
      current_lik = sum(exp(eta) - (y % eta));
      
      // Now, take derivatives
      first_deriv = sum((exp_eta - y) % X.col(j));
      second_deriv = sum(exp_eta % pow(X.col(j), 2));
      
      newton_dir = first_deriv / second_deriv;
      
      if (line_search) {
        
        t = 1.0;
        step_accepted = false;
        b_j_og = b(j);
        
        while(!step_accepted) {
          
          b(j) = b(j) - t * newton_dir;
          eta = X * b;
          f_proposed = sum(exp(eta) - (y % eta));
          
          if (f_proposed <= current_lik - alpha * t * first_deriv * newton_dir) {
            
            step_accepted = true;
            
          } else {
            
            t = beta * t;
            b(j) = b_j_og;
            
          }
          
        }
        
      } else {
        
        // take a full Newton step
        b(j) = b(j) - newton_dir;
        
      }
      
    }
    
  }
  
  return(b);
  
}

// Now, want to write a function that calls the above function in parallel
// I think that at first I can write a function just for updating loadings
// And then one for updating factors
// And once I figure those both out, I can try to combine them for code efficiency
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_loadings (
    const arma::mat& F_T,
    arma::mat& L,
    const arma::mat& Y,
    unsigned int update_start_idx,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for shared(F_T, Y, L, update_start_idx, num_iter, line_search, alpha, beta) 
  for (int i = 0; i < Y.n_rows; i++) {
    
    L.col(i) = solve_pois_reg_cpp (
      F_T, 
      conv_to< colvec >::from(Y.row(i)),
      L.col(i), 
      update_start_idx,
      num_iter,
      line_search,
      alpha,
      beta
    );
    
  }
  
  return(L);
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_factors (
    const arma::mat& L_T,
    arma::mat& FF,
    const arma::mat& Y,
    unsigned int update_start_idx,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for shared(L_T, Y, FF, update_start_idx, num_iter, line_search, alpha, beta) 
  for (int j = 0; j < Y.n_cols; j++) {
    
    FF.col(j) = solve_pois_reg_cpp (
      L_T, 
      Y.col(j),
      FF.col(j), 
      update_start_idx,
      num_iter,
      line_search,
      alpha,
      beta
    );
    
  }
  
  return(FF);
  
}
