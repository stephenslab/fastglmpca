#include <RcppArmadillo.h>
#include <omp.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;

inline arma::vec solve_bin_reg_cpp (
    const arma::mat X, 
    const arma::mat X_sqrd,
    const arma::vec y,
    const arma::vec size,
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
  current_lik = sum((size % log(1 + exp_eta)) - (y % eta));
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
    start_iter_lik = current_lik;
    
    for (int idx = 0; idx < num_indices; idx++) {
      
      j = update_indices[idx];
      
      // Now, take derivatives
      first_deriv = sum((size % X.col(j) % (exp_eta / (1 + exp_eta))) - (y % X.col(j)));
      second_deriv = sum(size % X_sqrd.col(j) % (exp_eta / pow((1 + exp_eta), 2)));
      
      newton_dir = first_deriv / second_deriv;
      
      if (line_search) {
        
        t = 1.0;
        step_accepted = false;
        newton_dec = alpha * first_deriv * newton_dir;
        b_j_og = b(j);
        
        while(!step_accepted) {
          
          b(j) = b_j_og - t * newton_dir;
          eta_proposed = eta + (b(j) - b_j_og) * X.col(j);
          //Rprintf("max eta_proposed = %f\n", eta_proposed.max());
          //Rprintf("min eta_proposed = %f\n", eta_proposed.min());
          exp_eta_proposed = exp(eta_proposed);
          f_proposed = sum((size % log(1 + exp_eta_proposed)) - (y % eta_proposed));
          //Rprintf("max exp_eta_proposed = %f\n", exp_eta_proposed.max());
          //Rprintf("min exp_eta_proposed = %f\n", exp_eta_proposed.min());
          
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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
arma::mat update_loadings_bin (
    const arma::mat& F_T,
    arma::mat& L,
    const arma::mat& Y_T,
    const arma::mat& N_T,
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
    
    L.col(i) = solve_bin_reg_cpp (
      F_T, 
      F_T_sqrd,
      Y_T.col(i),
      N_T.col(i),
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
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
arma::mat update_factors_bin (
    const arma::mat& L_T,
    arma::mat& FF,
    const arma::mat& Y,
    const arma::mat& N,
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
    
    FF.col(j) = solve_bin_reg_cpp (
      L_T, 
      L_T_sqrd,
      Y.col(j),
      N.col(j),
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


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
arma::mat update_loadings_missing_bin (
    const arma::mat& F_T,
    arma::mat& L,
    const arma::mat& Y_T,
    const arma::mat& N_T,
    Rcpp::List nonmissing_index_list,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta,
    const double ccd_iter_tol
) {
  
  const arma::mat F_T_sqrd = arma::pow(F_T, 2);
  
  //#pragma omp parallel for shared(F_T, Y_T, L, num_iter) 
  for (int i = 0; i < Y_T.n_cols; i++) {
    
    arma::uvec nonmissing_idx = Rcpp::as<arma::uvec>(nonmissing_index_list[i]);
    arma::vec y = Y_T.col(i);
    arma::vec n = N_T.col(i);
    
    L.col(i) = solve_bin_reg_cpp (
      F_T.rows(nonmissing_idx), 
      F_T_sqrd.rows(nonmissing_idx),
      y.elem(nonmissing_idx),
      n.elem(nonmissing_idx),
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
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
arma::mat update_factors_missing_bin (
    const arma::mat& L_T,
    arma::mat& FF,
    const arma::mat& Y,
    const arma::mat& N,
    Rcpp::List nonmissing_index_list,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta,
    const double ccd_iter_tol
) {
  
  const arma::mat L_T_sqrd = arma::pow(L_T, 2);
  
  //#pragma omp parallel for shared(L_T, Y, FF, num_iter) 
  for (int j = 0; j < Y.n_cols; j++) {
    
    arma::uvec nonmissing_idx = Rcpp::as<arma::uvec>(nonmissing_index_list[j]);
    arma::vec y = Y.col(j);
    arma::vec n = N.col(j);
    
    FF.col(j) = solve_bin_reg_cpp (
      L_T.rows(nonmissing_idx), 
      L_T_sqrd.rows(nonmissing_idx),
      y.elem(nonmissing_idx),
      n.elem(nonmissing_idx),
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

