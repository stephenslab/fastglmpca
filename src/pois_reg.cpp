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
    const double beta
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
  
  int num_indices = update_indices.size();
  current_lik = sum(exp_eta - (y % eta));
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
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
    
  }
  
  return(b);
  
}

inline arma::vec solve_pois_reg_faster_cpp (
    const arma::mat X, 
    const arma::vec m, 
    arma::vec b, 
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  double current_nonlinear_lik; // used to store log likelihood of each iteration
  double current_lik;
  double first_deriv;
  double second_deriv;
  double newton_dir;
  double newton_dec;
  arma::vec eta = X * b;
  arma::vec eta_proposed;
  arma::vec exp_deriv_term;
  arma::vec exp_eta_proposed;
  arma::vec exp_eta = exp(eta);
  double t;
  bool step_accepted;
  double f_proposed;
  int j;
  double b_j_og;
  
  int num_indices = update_indices.size();
  current_nonlinear_lik = sum(exp_eta);
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
    double lik_improvement = 0.0;
    
    for (int idx = 0; idx < num_indices; idx++) {
      
      j = update_indices[idx];
      
      // Now, take derivatives
      current_lik = current_nonlinear_lik - b(j) * m(idx);
      exp_deriv_term = exp_eta % X.col(j);
      first_deriv = sum(exp_deriv_term) - m(idx);
      second_deriv = sum(exp_deriv_term % X.col(j));
      
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
          f_proposed = sum(exp_eta_proposed) - b(j) * m(idx);
          
          if (f_proposed <= current_lik - t * newton_dec) {
            
            step_accepted = true;
            lik_improvement = lik_improvement + (current_lik - f_proposed);
            current_nonlinear_lik = f_proposed + b(j) * m(idx);
            eta = eta_proposed;
            exp_eta = exp_eta_proposed;
            
          } else {
            
            t = beta * t;
            
            if (abs((t * newton_dir) / b_j_og) < 1e-16) {
              
              b(j) = b_j_og;
              step_accepted = true;
              
            }
            
          }
          
        }
        
      } else {
        
        // take a full Newton step
        b(j) = b(j) - newton_dir;
        
      }
      
    }
    
    if (lik_improvement < 1e-8) {
      
      break;
      
    }
    
  }
  
  return(b);
  
}

inline arma::vec solve_pois_reg_cpp_calc_lik (
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
    double *end_lik
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
  
  int num_indices = update_indices.size();
  current_lik = sum(exp_eta - (y % eta));
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
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
    
  }
  
  *end_lik = f_proposed;
  return(b);
  
}

inline arma::vec solve_pois_reg_faster_calc_lik_cpp (
    const arma::mat X, 
    const arma::vec m, 
    arma::vec b, 
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta,
    double *end_lik
) {
  
  double current_nonlinear_lik; // used to store log likelihood of each iteration
  double current_lik;
  double first_deriv;
  double second_deriv;
  double newton_dir;
  double newton_dec;
  arma::vec eta = X * b;
  arma::vec eta_proposed;
  arma::vec exp_deriv_term;
  arma::vec exp_eta_proposed;
  arma::vec exp_eta = exp(eta);
  double t;
  bool step_accepted;
  double f_proposed;
  int j;
  double b_j_og;
  
  int num_indices = update_indices.size();
  current_nonlinear_lik = sum(exp_eta);
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
    double lik_improvement = 0.0;
    
    for (int idx = 0; idx < num_indices; idx++) {
      
      j = update_indices[idx];
      
      // Now, take derivatives
      current_lik = current_nonlinear_lik - b(j) * m(idx);
      exp_deriv_term = exp_eta % X.col(j);
      first_deriv = sum(exp_deriv_term) - m(idx);
      second_deriv = sum(exp_deriv_term % X.col(j));
      
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
          f_proposed = sum(exp_eta_proposed) - b(j) * m(idx);
          
          if (f_proposed <= current_lik - t * newton_dec) {
            
            step_accepted = true;
            lik_improvement = lik_improvement + (current_lik - f_proposed);
            current_nonlinear_lik = f_proposed + b(j) * m(idx);
            eta = eta_proposed;
            exp_eta = exp_eta_proposed;
            
          } else {
            
            t = beta * t;
            
            if (abs((t * newton_dir) / b_j_og) < 1e-16) {
              
              b(j) = b_j_og;
              step_accepted = true;
              
            }
            
          }
          
        }
        
      } else {
        
        // take a full Newton step
        b(j) = b(j) - newton_dir;
        
      }
      
    }
    
    // if likelihood didn't improve much after cycling through all params
    // move on
    if (lik_improvement < 1e-8) {
      
      break;
      
    }
    
  }
  
  double end_linear_lik = 0.0;
  
  for (int idx = 0; idx < num_indices; idx++) {
    
    j = update_indices[idx];
    end_linear_lik = end_linear_lik + b(j) * m(idx);
    
  }
  
  *end_lik = current_nonlinear_lik - end_linear_lik;
  return(b);
  
}

// Now, want to write a function that calls the above function in parallel
// I think that at first I can write a function just for updating loadings
// And then one for updating factors
// And once I figure those both out, I can try to combine them for code efficiency
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
void update_loadings (
    const arma::mat& F_T,
    arma::mat& L,
    const arma::mat& Y_T,
    const arma::mat& deriv_const_mat,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
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
      beta
    );
    
  }
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
void update_loadings_faster (
    const arma::mat& F_T,
    arma::mat& L,
    const arma::mat& M,
    const std::vector<int> update_indices,
    const int n,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for shared(F_T, M, L, num_iter) 
  for (int i = 0; i < n; i++) {
    
    L.col(i) = solve_pois_reg_faster_cpp (
      F_T, 
      M.col(i),
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
double update_factors (
    const arma::mat& L_T,
    arma::mat& FF,
    const arma::mat& Y,
    const arma::mat& deriv_const_mat,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  const arma::mat L_T_sqrd = arma::pow(L_T, 2);
  double end_lik_sum = 0.0;
  
  #pragma omp parallel for shared(L_T, Y, FF, num_iter) reduction(+:end_lik_sum)
  for (int j = 0; j < Y.n_cols; j++) {
    
    double* end_lik = (double*) malloc(sizeof(double));
    *end_lik = 0.0;
    
    FF.col(j) = solve_pois_reg_cpp_calc_lik (
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
      end_lik
    );
    
    // subtract end likelihood because inline function minimizes
    // the negative log-likelihood
    end_lik_sum -= *end_lik;
    free(end_lik);
    
  }
  
  return(end_lik_sum);
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
double update_factors_faster (
    const arma::mat& L_T,
    arma::mat& FF,
    const arma::mat& M,
    const std::vector<int> update_indices,
    const int p,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  double end_lik_sum = 0.0;
  
  #pragma omp parallel for shared(L_T, M, FF, num_iter) reduction(+:end_lik_sum)
  for (int j = 0; j < p; j++) {
    
    double* end_lik = (double*) malloc(sizeof(double));
    *end_lik = 0.0;
    
    FF.col(j) = solve_pois_reg_faster_calc_lik_cpp (
      L_T, 
      M.col(j),
      FF.col(j), 
      update_indices,
      num_iter,
      line_search,
      alpha,
      beta,
      end_lik
    );
    
    // subtract end likelihood because inline function minimizes
    // the negative log-likelihood
    end_lik_sum -= *end_lik;
    free(end_lik);
    
  }
  
  return(end_lik_sum);
  
}

