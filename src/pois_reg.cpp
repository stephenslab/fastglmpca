#include <RcppArmadillo.h>
#include <omp.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;

inline arma::vec solve_pois_reg_cpp (
    const arma::mat X, 
    const arma::mat X_sqrd,
    const arma::vec y, 
    const arma::vec deriv_const,
    arma::vec b, 
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta,
    const double ccd_iter_tol
  ) {
  
  double current_lik; // used to store log likelihood of each iteration
  double first_deriv;
  double second_deriv;
  double newton_dir;
  double newton_dec;
  arma::vec eta = X * b;
  arma::vec eta_proposed;
  arma::vec exp_eta_proposed;
  arma::vec exp_eta = exp(eta);
  double t;
  bool step_accepted;
  double f_proposed;
  int j;
  double b_j_og;
  double start_iter_lik;
  
  int num_indices = update_indices.size();
  current_lik = sum(exp_eta - (y % eta));
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
    start_iter_lik = current_lik;
    
    for (int idx = 0; idx < num_indices; idx++) {
      
      j = update_indices[idx];
      
      // Now, take derivatives
      first_deriv = sum(exp_eta % X.col(j)) + deriv_const(j);
      second_deriv = sum(exp_eta % X_sqrd.col(j));
      
      newton_dir = first_deriv / second_deriv;
      
      if (line_search) {
        
        t = 1.0;
        step_accepted = false;
        newton_dec = alpha * first_deriv * newton_dir;
        b_j_og = b(j);

        while(!step_accepted) {
          
          b(j) = b_j_og - t * newton_dir;
          eta_proposed = eta + (b(j) - b_j_og) * X.col(j);
          exp_eta_proposed = exp(eta_proposed);
          f_proposed = sum(exp_eta_proposed - (y % eta_proposed));
          
          if (f_proposed <= current_lik - t * newton_dec) {
            
            step_accepted = true;
            current_lik = f_proposed;
            eta = eta_proposed;
            exp_eta = exp_eta_proposed;
            
          } else {
            
            t = beta * t;
            
          }
          
        }
        
      } else {
        
        // take a full Newton step
        b(j) = b(j) - newton_dir;
        
      }
      
    }
    
    // if the negative log likelihood didn't increase substantially, break
    if (start_iter_lik - current_lik < ccd_iter_tol) {
      
      break;
      
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
    const arma::mat& Y_T,
    const arma::mat& deriv_const_mat,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta,
    const double ccd_iter_tol
) {
  
  const arma::mat F_T_sqrd = arma::pow(F_T, 2);
  
  #pragma omp parallel for shared(F_T, Y_T, L, num_iter) 
  for (int i = 0; i < Y_T.n_cols; i++) {
    
    L.col(i) = solve_pois_reg_cpp (
      F_T, 
      F_T_sqrd,
      Y_T.col(i),
      deriv_const_mat.col(i),
      L.col(i), 
      update_indices,
      num_iter,
      line_search,
      alpha,
      beta,
      ccd_iter_tol
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
    const arma::mat& deriv_const_mat,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta,
    const double ccd_iter_tol
) {
  
  const arma::mat L_T_sqrd = arma::pow(L_T, 2);
  
  #pragma omp parallel for shared(L_T, Y, FF, num_iter) 
  for (int j = 0; j < Y.n_cols; j++) {
    
    FF.col(j) = solve_pois_reg_cpp (
      L_T, 
      L_T_sqrd,
      Y.col(j),
      deriv_const_mat.col(j),
      FF.col(j), 
      update_indices,
      num_iter,
      line_search,
      alpha,
      beta,
      ccd_iter_tol
    );
    
  }
  
  return(FF);
  
}
