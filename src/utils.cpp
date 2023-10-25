#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace arma;
using namespace Rcpp;
using namespace RcppParallel;

// Y is n x m
// L is K x n
// F is K x m
//
// This function returns
//
//  H = crossprod(LL,FF)
//  sum(Y*H - exp(H)) - loglik_const
//
// but the computation is done in a memory-efficient way.
//
// [[Rcpp::export]]
double lik_glmpca_pois_log_sp (const arma::sp_mat& Y,
			       const arma::mat& L, 
			       const arma::mat& F, 
			       double loglik_const) {
  unsigned int n = Y.n_rows;
  unsigned int m = Y.n_cols;
  unsigned int j;
  vec y(n);
  vec lf(n);
  vec logliks(m,fill::zeros);

  // Repeat for each column of Y.
  for (j = 0; j < m; j++) {
    y  = Y.col(j);
    lf = L.t() * F.col(j);
    logliks.at(j) = sum(y % lf) - sum(exp(lf));
  }

  return sum(logliks) - loglik_const;
}

// L is K x n
// F is K x m
//
// This function returns 
//
//   L %*% exp(t(L) %*% F)
// 
// but the computation is done in a memory-efficient way.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat deriv_prod (const arma::mat& L, const arma::mat& F) {
  const unsigned int n = L.n_cols;
  const unsigned int m = F.n_cols;
  const unsigned int K = L.n_rows;
  unsigned int i, k, j;
  mat out(K,m,fill::zeros);
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      for (k = 0; k < K; k++)
        out.at(k,j) += L.at(k,i) * exp(dot(L.col(i),F.col(j)));
  return(out);
}
