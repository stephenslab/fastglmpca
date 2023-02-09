#include <RcppArmadillo.h>
#include <omp.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;

inline arma::vec solve_pois_reg_irls (
    const arma::mat X, 
    const arma::vec y, 
    arma::vec b, 
    unsigned int num_iter,
    const arma::vec update_vec,
    const arma::vec constant_vec
) {
  
  // Note, this code may run faster with malloc
  // Because I know the exact size of each of these elements
  arma::vec eta;
  arma::vec mu;
  arma::vec inv_mu;
  arma::vec z;
  arma::mat W_sqrt;
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
    eta = X * b;
    mu = exp(eta);
    inv_mu = 1 / mu;
    z = eta + (y - mu) % inv_mu;
    
    W_sqrt = arma::diagmat(sqrt(mu));
    b = arma::solve(W_sqrt * X, W_sqrt * z) % update_vec + b % constant_vec; 
        
  }
  
  return(b);
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_loadings_irls (
    const arma::mat& F_T,
    arma::mat& L,
    const arma::mat& Y_T,
    const arma::vec& update_vec,
    const arma::vec& constant_vec,
    unsigned int num_iter
) {

  #pragma omp parallel for shared(F_T, Y_T, L, num_iter)
  for (int i = 0; i < Y_T.n_cols; i++) {
    
    L.col(i) = solve_pois_reg_irls (
      F_T,
      Y_T.col(i),
      L.col(i),
      num_iter,
      update_vec,
      constant_vec
    );

  }

  return(L);

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_factors_irls (
    const arma::mat& L_T,
    arma::mat& FF,
    const arma::mat& Y,
    const arma::vec& update_vec,
    const arma::vec& constant_vec,
    unsigned int num_iter
) {

  #pragma omp parallel for shared(L_T, Y, FF, num_iter)
  for (int j = 0; j < Y.n_cols; j++) {

    FF.col(j) = solve_pois_reg_irls (
      L_T,
      Y.col(j),
      FF.col(j),
      num_iter,
      update_vec,
      constant_vec
    );

  }

  return(FF);

}

