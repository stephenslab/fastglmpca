#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// L is K x n
// F is K x p
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double big_exp_crossprod(
  const arma::mat& L,
  const arma::mat& F,
  const int n,
  const int p
) {

  double sum = 0.0;

  // First, get this code to work without parallelism
  // Then, I can see if that would help
  for (int i = 0; i < n; i++)
    for (int j = 0; j < p; j++)
      sum = sum + exp(dot(L.col(i), F.col(j)));
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
    sum = sum + nonzero_y[r] * dot(L.col(nonzero_y_i_idx[r]), F.col(nonzero_y_j_idx[r]));
  return(sum);
}

// L is K x n
// F is K x p
// compute L %*% exp(t(L) %*% F)
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat deriv_product(const arma::mat& L, const arma::mat& F) {

  const int n = L.n_cols;
  const int K = L.n_rows;
  const int p = F.n_cols;
  mat prod;
  prod.zeros(K, p);
  double exp_arg;

  for (int k = 0; k < n; k++) {
    for (int i = 0; i < K; i++) {
      for (int j = 0; j < p; j++) {
        exp_arg = 0;
        for (int m = 0; m < K; m++) {
          exp_arg += L(m, k) * F(m, j);
        }

        prod(i, j) += L(i, k) * exp(exp_arg);
      }
    }
  }

  return(prod);

}
