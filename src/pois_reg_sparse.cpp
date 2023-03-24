#include <RcppArmadillo.h>
#include <omp.h>
#include <Rcpp.h>
#include <Rinterface.h>
using namespace arma;
using namespace Rcpp;

inline arma::vec solve_pois_reg_cpp_sp (
    const arma::mat X, 
    const arma::mat X_sqrd,
    const arma::sp_vec y, 
    const arma::sp_vec deriv_const,
    arma::vec b, 
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  double first_deriv;
  double second_deriv;
  double newton_dir;
  double newton_const;
  arma::vec eta = X * b;
  arma::vec eta_proposed;
  arma::vec exp_eta = exp(eta);
  double *current_lik;
  current_lik = (double *) malloc(sizeof(double));
  arma::vec old_b_sub;
  double t;
  bool step_accepted;
  double b_j_og;
  double *f_proposed;
  f_proposed = (double *) malloc(sizeof(double));
  int j;
  
  *current_lik = sum(exp_eta) - arma::dot(y, eta);
  
  int num_indices = update_indices.size();
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
    for (int idx = 0; idx < num_indices; idx++) {
      
      j = update_indices[idx];
      
      // Take derivatives
      first_deriv = sum(exp_eta % X.col(j)) + deriv_const(j);
      second_deriv = sum(exp_eta % X_sqrd.col(j));
      
      newton_dir = first_deriv / second_deriv;
      
      if (line_search) {
        
        t = 1.0;
        newton_const = alpha * first_deriv * newton_dir;
        step_accepted = false;
        b_j_og = b(j);
        old_b_sub = b_j_og * X.col(j);
        
        while(!step_accepted) {
          
          Rprintf("%f\n", t);
          b(j) = b_j_og - t * newton_dir;
          //eta = X * b;
          eta_proposed = eta - old_b_sub + b(j) * X.col(j);
          *f_proposed = sum(exp(eta_proposed)) - arma::dot(y, eta_proposed);
          
          if (*f_proposed <= *current_lik - t * newton_const) {
            
            step_accepted = true;
            
          } else {
            
            t = beta * t;
            
          }
          
        }
        
        // here, want to update params that have been computed
        *current_lik = *f_proposed;
        eta = eta_proposed;
        exp_eta = exp(eta);
        
      } else {
        
        // take a full Newton step
        b(j) = b(j) - newton_dir;
        eta = X * b;
        exp_eta = exp(eta);
        
      }
      
    }
    
  }
  
  free(current_lik);
  free(f_proposed);
  return(b);
  
}

// Now, want to write a function that calls the above function in parallel
// I think that at first I can write a function just for updating loadings
// And then one for updating factors
// And once I figure those both out, I can try to combine them for code efficiency
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_loadings_sp (
    const arma::mat& F_T,
    const arma::mat& F_T_sqrd,
    arma::mat& L,
    const arma::sp_mat& Y_T,
    const arma::sp_mat& deriv_const_mat,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for shared(F_T, F_T_sqrd, Y_T, L, num_iter) 
  for (int i = 0; i < Y_T.n_cols; i++) {
    
    L.col(i) = solve_pois_reg_cpp_sp (
      F_T, 
      F_T_sqrd,
      Y_T.col(i),
      deriv_const_mat.col(i),
      L.col(i), 
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
arma::mat update_factors_sp (
    const arma::mat& L_T,
    const arma::mat& L_T_sqrd,
    arma::mat& FF,
    const arma::sp_mat& Y,
    const arma::sp_mat& deriv_const_mat,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for shared(L_T, L_T_sqrd, Y, FF, num_iter) 
  for (int j = 0; j < Y.n_cols; j++) {
    
    FF.col(j) = solve_pois_reg_cpp_sp (
      L_T, 
      L_T_sqrd,
      Y.col(j),
      deriv_const_mat.col(j),
      FF.col(j), 
      update_indices,
      num_iter,
      line_search,
      alpha,
      beta
    );
    
  }
  
  return(FF);
  
}
