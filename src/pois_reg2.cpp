#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace arma;
using namespace Rcpp;
using namespace RcppParallel;

arma::vec solve_pois_reg_faster_cpp (
    const arma::mat X, 
    const arma::vec m, 
    arma::vec b, 
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const bool line_search,
    const double alpha,
    const double beta) {
  
  // Used to store log likelihood of each iteration.
  double current_nonlinear_lik; 

  double current_lik;
  double first_deriv;
  double second_deriv;
  double newton_dir;
  double newton_dec;
  vec eta = X * b;
  vec eta_proposed;
  vec exp_deriv_term;
  vec exp_eta_proposed;
  vec exp_eta = exp(eta);
  double t;
  double f_proposed;
  unsigned int i, j;
  double b_j_og;
  
  int num_indices = update_indices.size();
  current_nonlinear_lik = sum(exp_eta);
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    double lik_improvement = 0;
    for (i = 0; i < num_indices; i++) {
      j = update_indices[i];
      
      // Now, take derivatives.
      current_lik    = current_nonlinear_lik - b[j] * m[i];
      exp_deriv_term = exp_eta % X.col(j);
      first_deriv    = sum(exp_deriv_term) - m[i];
      second_deriv   = dot(exp_deriv_term,X.col(j));
      newton_dir     = first_deriv / second_deriv;
      
      if (line_search) {
        t             = 1;
        newton_dec    = alpha * first_deriv * newton_dir;
        b_j_og        = b[j];
        while (true) {
          b[j]             = b_j_og - t * newton_dir;
          eta_proposed     = eta + (b[j] - b_j_og) * X.col(j);
          exp_eta_proposed = exp(eta_proposed);
          f_proposed       = sum(exp_eta_proposed) - b[j] * m[i];
          if (f_proposed <= current_lik - t*newton_dec) {
            lik_improvement = lik_improvement + (current_lik - f_proposed);
            current_nonlinear_lik = f_proposed + b[j] * m[i];
            eta                   = eta_proposed;
            exp_eta               = exp_eta_proposed;
	    break;
          } else {
            t *= beta;

	    // Terminate backtracking line search because we have
	    // arrived at the smallest allowable step size.
	    if (t < 1e-4) {
	      t = 0;
	      break;
	    }
          }
        }
      } else
        
        // Take a full Newton step.
        b[j] -= newton_dir;
    }
    
    if (lik_improvement < 1e-8)
      break;
  }

  return(b);
}

struct FactorsUpdater : public Worker {
  const mat& L_T;
  const mat& M;
  mat& FF;
  const std::vector<int> update_indices;
  const unsigned int num_iter;
  const bool line_search;
  const double alpha;
  const double beta;
  
  FactorsUpdater (const mat& L_T, const mat& M, mat& FF,
                  const std::vector<int> update_indices, 
		  unsigned int num_iter, bool line_search, 
		  double alpha, double beta) 
    : L_T(L_T), M(M), FF(FF), update_indices(update_indices), 
      num_iter(num_iter), line_search(line_search), alpha(alpha), 
      beta(beta) { }
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++)
      FF.col(j) = solve_pois_reg_faster_cpp(L_T,M.col(j),FF.col(j),
					    update_indices,num_iter,
					    line_search,alpha,beta);
  }
};

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
void update_factors_faster_parallel(const arma::mat& L_T, arma::mat& FF, 
				    const arma::mat& M, 
				    const std::vector<int> update_indices, 
				    unsigned int num_iter, bool line_search, 
				    double alpha, double beta) {
  FactorsUpdater updater(L_T,M,FF,update_indices,num_iter,line_search,alpha,
			 beta);
  parallelFor(0,FF.n_cols,updater);
}
