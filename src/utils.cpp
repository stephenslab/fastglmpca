#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace arma;
using namespace Rcpp;
using namespace RcppParallel;

// L is K x n
// F is K x p
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double big_exp_crossprod (const arma::mat& L, const arma::mat& F,
			  const int n, const int m) {
  double sum = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      sum = sum + exp(dot(L.col(i),F.col(j)));
  return(sum);
}

// L is K x n
// F is K x p
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double big_elementwise_mult_crossprod(
    const arma::mat& L,
    const arma::mat& F,
    const arma::vec& nonzero_y,
    const std::vector<int> nonzero_y_i_idx,
    const std::vector<int> nonzero_y_j_idx,
    const int num_nonzero_y) {
  double sum = 0;
  for (int r = 0; r < num_nonzero_y; r++)
    sum = sum + nonzero_y[r] * dot(L.col(nonzero_y_i_idx[r]),
				   F.col(nonzero_y_j_idx[r]));
  return(sum);
}

// L is K x n
// F is K x m
// compute L %*% exp(t(L) %*% F)
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat deriv_prod (const arma::mat& L, const arma::mat& F) {
  const unsigned int n = L.n_cols;
  const unsigned int m = F.n_cols;
  const unsigned int K = L.n_rows;
  unsigned int i, k, j;
  mat out;
  out.zeros(K,m);
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      for (k = 0; k < K; k++)
        out[k,j] += L[k,i] * exp(dot(L.col(i),F.col(j)));
  return(out);
}
