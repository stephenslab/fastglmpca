#include <RcppArmadillo.h>
#include <omp.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;

inline arma::vec solve_pois_reg_logsp_cpp (
    const arma::mat X, 
    const arma::vec y, 
    arma::vec b, 
    const double s,
    const std::vector<int> update_indices,
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
  int j;
  
  int num_indices = update_indices.size();
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
    for (int idx = 0; idx < num_indices; idx++) {
      
      j = update_indices[idx];
      eta = X * b;
      exp_eta = exp(eta);
      current_lik = sum(exp_eta - (y % log(exp_eta - s)));
      
      // Now, take derivatives
      first_deriv = sum((1 - (y / (exp_eta - s))) % exp_eta % X.col(j));
      second_deriv = sum((1 - (y / (exp_eta - s)) + (y / pow(exp_eta - s, 2)) % exp_eta) % exp_eta % pow(X.col(j), 2));
      
      newton_dir = first_deriv / second_deriv;
      
      if (line_search) {
        
        // start line search such that update remains non-negative
        if (newton_dir > 0) {
          
          // If search direction is negative and value is already at the boundary,
          // keep going in the loop
          if (b(j) <= 1e-12) {
            
            continue;
            
          }
          
          t = std::min(1.0, (b(j) - 1e-12) / newton_dir);
          
        } else {
          
          t = 1.0;
          
        }
        
        step_accepted = false;
        b_j_og = b(j);
        
        while(!step_accepted) {
          
          b(j) = b(j) - t * newton_dir;
          eta = X * b;
          exp_eta = exp(eta);
          f_proposed = sum(exp_eta - (y % log(exp_eta - 1)));
          
          if (f_proposed <= current_lik - alpha * t * first_deriv * newton_dir) {
            
            step_accepted = true;
            
          } else {
            
            t = beta * t;
            b(j) = b_j_og;
            
          }
          
        }
        
      } else {
        
        // take a full Newton step, thresholding at 0
        b(j) = std::max(b(j) - newton_dir, 1e-12);
        
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
arma::mat update_loadings_logsp (
    const arma::mat& F_T,
    arma::mat& L,
    const arma::mat& Y_T,
    const double s,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for shared(F_T, Y_T, L, num_iter) 
  for (int i = 0; i < Y_T.n_cols; i++) {
    
    L.col(i) = solve_pois_reg_logsp_cpp (
      F_T, 
      Y_T.col(i),
      L.col(i), 
      s,
      update_indices,
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
arma::mat update_factors_logsp (
    const arma::mat& L_T,
    arma::mat& FF,
    const arma::mat& Y,
    const double s,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for shared(L_T, Y, FF, num_iter) 
  for (int j = 0; j < Y.n_cols; j++) {
    
    FF.col(j) = solve_pois_reg_logsp_cpp (
      L_T, 
      Y.col(j),
      FF.col(j), 
      s,
      update_indices,
      num_iter,
      line_search,
      alpha,
      beta
    );
    
  }
  
  return(FF);
  
}
