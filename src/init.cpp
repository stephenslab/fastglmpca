#include <RcppArmadillo.h>
#include <omp.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;

inline arma::vec approx_pois_reg_expected_ll_cpp (
    const int n,
    const arma::vec linear_term_const,
    arma::vec b, 
    double fixed_factor,
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
      
      current_lik = fixed_factor * exp_term - b(j) * linear_term_const(j);
      first_deriv = fixed_factor * exp_term * b(j) - linear_term_const(j);
      second_deriv = fixed_factor * exp_term * (1 + b(j) * b(j));
      
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
          f_proposed = n * fixed_factor * exp(.5 * sum_b_sqrd) - b(j) * linear_term_const(j);
          
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


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
void update_factors_approx (
    unsigned int n,
    unsigned int p,
    arma::mat& FF,
    const arma::vec& fixed_factor_vec,
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
      fixed_factor_vec(j),
      update_indices,
      num_iter,
      line_search,
      alpha,
      beta
    );
    
  }
  
}