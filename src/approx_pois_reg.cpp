#include <RcppArmadillo.h>
#include <omp.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;

inline arma::vec approx_pois_reg_expected_ll_cpp (
    const int n,
    const arma::vec linear_term_const,
    arma::vec b, 
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
  double newton_dec;
  double t;
  bool step_accepted;
  double f_proposed;
  int j;
  double b_j_og;
  double sum_b_sqrd = arma::accu(arma::square(b));
  double exp_term;
  double sum_b_sqrd_og;
  
  int num_indices = update_indices.size();
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
    for (int idx = 0; idx < num_indices; idx++) {
      
      j = update_indices[idx];
      
      exp_term = n * exp(.5 * sum_b_sqrd);
      
      current_lik = exp_term - b(j) * linear_term_const(j);
      first_deriv = exp_term * b(j) - linear_term_const(j);
      second_deriv = exp_term * (1 + b(j) * b(j));
      
      newton_dir = first_deriv / second_deriv;
      
      if (line_search) {
        
        t = 1.0;
        step_accepted = false;
        newton_dec = alpha * first_deriv * newton_dir;
        b_j_og = b(j);
        sum_b_sqrd_og = sum_b_sqrd - (b_j_og * b_j_og);
        
        while(!step_accepted) {
          
          b(j) = b_j_og - t * newton_dir;
          sum_b_sqrd = sum_b_sqrd_og + (b(j) * b(j));
          f_proposed = n * exp(.5 * sum_b_sqrd) - b(j) * linear_term_const(j);
          
          if (f_proposed <= current_lik - t * newton_dec) {
            
            step_accepted = true;
            
          } else {
            
            t = beta * t;
            
            if (abs((t * newton_dir) / b_j_og) < 1e-16) {
              
              b(j) = b_j_og;
              step_accepted = true;
              sum_b_sqrd = sum_b_sqrd_og + (b_j_og * b_j_og);
              
            }
            
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

inline arma::vec approx_pois_reg_expected_ll_nz_mean_cpp (
    const int n,
    const arma::vec linear_term_const,
    arma::vec b, 
    const arma::vec x_mean,
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
  double newton_dec;
  double t;
  bool step_accepted;
  double f_proposed;
  int j;
  double b_j_og;
  double sum_b_sqrd = arma::accu(arma::square(b));
  double b_x_ip = arma::dot(b, x_mean);
  double exp_term;
  double sum_b_sqrd_og;
  double b_x_ip_og;
  double sum_b_x;
  
  int num_indices = update_indices.size();
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
    for (int idx = 0; idx < num_indices; idx++) {
      
      j = update_indices[idx];
      
      exp_term = n * exp(.5 * sum_b_sqrd + b_x_ip);
      sum_b_x = b(j) + x_mean(j);
      
      current_lik = exp_term - b(j) * linear_term_const(j);
      first_deriv = exp_term * sum_b_x - linear_term_const(j);
      second_deriv = exp_term * (1 + sum_b_x * sum_b_x);
      
      newton_dir = first_deriv / second_deriv;
      
      if (line_search) {
        
        t = 1.0;
        step_accepted = false;
        newton_dec = alpha * first_deriv * newton_dir;
        b_j_og = b(j);
        sum_b_sqrd_og = sum_b_sqrd - (b_j_og * b_j_og);
        b_x_ip_og = b_x_ip - b(j) * x_mean(j);
        
        while(!step_accepted) {
          
          b(j) = b_j_og - t * newton_dir;
          sum_b_sqrd = sum_b_sqrd_og + (b(j) * b(j));
          b_x_ip = b_x_ip_og + b(j) * x_mean(j);
          f_proposed = n * exp(.5 * sum_b_sqrd + b_x_ip) - b(j) * linear_term_const(j);
          
          if (f_proposed <= current_lik - t * newton_dec) {
            
            step_accepted = true;
            
          } else {
            
            t = beta * t;
            
            if (abs((t * newton_dir) / b_j_og) < 1e-16) {
              
              b(j) = b_j_og;
              step_accepted = true;
              sum_b_sqrd = sum_b_sqrd_og + (b_j_og * b_j_og);
              b_x_ip = b_x_ip_og + b_j_og * x_mean(j);
              
            }
            
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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
void update_loadings_approx (
    unsigned int n,
    unsigned int p,
    arma::mat& L,
    const arma::mat& linear_term_mat,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for shared(L, n, p, num_iter) 
  for (int i = 0; i < n; i++) {
    
    L.col(i) = approx_pois_reg_expected_ll_cpp (
      p,
      linear_term_mat.col(i),
      L.col(i), 
      update_indices,
      num_iter,
      line_search,
      alpha,
      beta
    );
    
  }
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
void update_factors_approx (
    unsigned int n,
    unsigned int p,
    arma::mat& FF,
    const arma::mat& linear_term_mat,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for shared(FF, n, p, num_iter) 
  for (int j = 0; j < p; j++) {
    
    FF.col(j) = approx_pois_reg_expected_ll_cpp (
      n,
      linear_term_mat.col(j),
      FF.col(j), 
      update_indices,
      num_iter,
      line_search,
      alpha,
      beta
    );
    
  }
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
void update_loadings_approx_nz_mean (
    unsigned int n,
    unsigned int p,
    arma::mat& L,
    const arma::vec& FF_row_means,
    const arma::mat& linear_term_mat,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for shared(L, n, p, num_iter) 
  for (int i = 0; i < n; i++) {
    
    L.col(i) = approx_pois_reg_expected_ll_nz_mean_cpp (
      p,
      linear_term_mat.col(i),
      L.col(i), 
      FF_row_means,
      update_indices,
      num_iter,
      line_search,
      alpha,
      beta
    );
    
  }
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
void update_factors_approx_nz_mean (
    unsigned int n,
    unsigned int p,
    arma::mat& FF,
    const arma::vec& LL_row_means,
    const arma::mat& linear_term_mat,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for shared(FF, n, p, num_iter, LL_row_means) 
  for (int j = 0; j < p; j++) {
    
    FF.col(j) = approx_pois_reg_expected_ll_nz_mean_cpp (
      n,
      linear_term_mat.col(j),
      FF.col(j), 
      LL_row_means,
      update_indices,
      num_iter,
      line_search,
      alpha,
      beta
    );
    
  }
  
}
