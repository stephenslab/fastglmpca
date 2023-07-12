#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

inline arma::vec approx_reg_ccd_cpp (
    const arma::mat X,
    const arma::vec y,
    const arma::uvec non_zero_idx,
    const arma::uvec zero_idx,
    arma::vec w,
    const double a1,
    const double a2,
    const arma::uvec update_indices,
    const int n_updates,
    const double alpha,
    const double beta
) {
  
  arma::mat X_0 = arma::conv_to<arma::mat>::from(X.rows(zero_idx));
  arma::mat X_nz = arma::conv_to<arma::mat>::from(X.rows(non_zero_idx));
  arma::mat Z_0 = X_0.t() * X_0;
  arma::vec a1_bar(zero_idx.size()); //
  a1_bar.fill(-a1);
  arma::mat m_0 = X_0.cols(update_indices).t() * a1_bar;
  arma::mat m_nz = X_nz.cols(update_indices).t() * y;
  arma::mat m = m_0 + m_nz;
  double sqrd_factor = 2 * a2;
  
  arma::vec eta = X_nz * w;
  arma::vec exp_eta = exp(eta);
  
  double current_nonlinear_lik = a2 * as_scalar(w.t() * Z_0 * w) + sum(exp_eta);
  
  for (int updates = 1; updates <= n_updates; updates++) {
    for (arma::uword idx = 0; idx < update_indices.n_elem; idx++) {
      arma::uword j = update_indices(idx);
      
      double current_lik = current_nonlinear_lik - w(j) * m(idx);
      arma::vec deriv_prod = exp_eta % X_nz.col(j);
      double first_deriv = sqrd_factor * as_scalar(w.t() * Z_0.col(j)) + sum(deriv_prod) - m(idx);
      double second_deriv = sqrd_factor * Z_0(j, j) + sum(deriv_prod % X_nz.col(j));
      
      double newton_dir = first_deriv / second_deriv;
      
      double t = 1.0;
      bool step_accepted = false;
      double newton_dec = alpha * first_deriv * newton_dir;
      double w_j_og = w(j);
      arma::vec eta_og = eta;
      arma::vec exp_eta_og = exp_eta;
      
      while (!step_accepted) {

        w(j) = w_j_og - t * newton_dir;
        arma::vec eta_proposed = eta + (w(j) - w_j_og) * X_nz.col(j);
        arma::vec exp_eta_proposed = exp(eta_proposed);
        double f_proposed = a2 * as_scalar(w.t() * Z_0 * w) + sum(exp_eta_proposed) - w(j) * m(idx);
        
        if (f_proposed <= current_lik - t * newton_dec) {
          step_accepted = true;
          current_nonlinear_lik = f_proposed + w(j) * m(idx);
          eta = eta_proposed;
          exp_eta = exp_eta_proposed;
        } else {
          t = beta * t;
          if (std::abs((t * newton_dir) / w_j_og) < 1e-16) {
            w(j) = w_j_og;
            step_accepted = true;
            eta = eta_og;
            exp_eta = exp_eta_og;
          }
        }
      }
    }
  }
  
  return(w);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
void update_loadings_approx_reg (
    const arma::mat& F_T,
    arma::mat& L,
    const Rcpp::List non_zero_Y_idx_by_row,
    const Rcpp::List non_zero_Y_by_row,
    arma::uvec full_n_indices,
    const double a1,
    const double a2,
    const int n,
    const arma::uvec& update_indices,
    const int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  for (int i = 0; i < n; i++) {
    
    arma::uvec rows_to_shed = non_zero_Y_idx_by_row[i];
    arma::uvec full_indices = full_n_indices;
    
    for(uword l = rows_to_shed.n_elem; l>0; l--) {
      
      full_indices.shed_row(l);
      
    }
    
    // Note, I may have to do some casting here
    L.col(i) = approx_reg_ccd_cpp (
      F_T,
      non_zero_Y_by_row[i],
      non_zero_Y_idx_by_row[i],
      full_indices,
      L.col(i),
      a1,
      a2,
      update_indices,
      num_iter,
      alpha,
      beta
    );
    
  }
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
void update_factors_approx_reg (
    const arma::mat& L_T,
    arma::mat& FF,
    const Rcpp::List non_zero_Y_idx_by_col,
    const Rcpp::List non_zero_Y_by_col,
    arma::uvec full_p_indices,
    const double a1,
    const double a2,
    const int p,
    const arma::uvec& update_indices,
    const int num_iter,
    const bool line_search,
    const double alpha,
    const double beta
) {
  
  #pragma omp parallel for private(full_p_indices)
  for (int j = 0; j < p; j++) {
    
    arma::uvec rows_to_shed = non_zero_Y_idx_by_col[j];
    arma::uvec full_indices = full_p_indices; 

    for(uword l = rows_to_shed.n_elem; l>0; l--) {
      
      full_indices.shed_row(l);
      
    }
    
    // Note, I may have to do some casting here
    FF.col(j) = approx_reg_ccd_cpp (
      L_T,
      non_zero_Y_by_col[j],
      non_zero_Y_idx_by_col[j],
      full_indices,
      FF.col(j),
      a1,
      a2,
      update_indices,
      num_iter,
      alpha,
      beta
    );
    
  }
  
}

