#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace arma;
using namespace Rcpp;
using namespace RcppParallel;

double compute_loglik_glmpca_pois (const sp_mat& Y, const mat& L, 
				   const mat& F, unsigned int j) {
  unsigned int n = Y.n_rows;
  vec y(n);
  vec lf(n);
  y  = Y.col(j);
  lf = L.t() * F.col(j);
  return dot(y,lf) - sum(exp(lf));
}

struct lik_glmpca_pois_log_sp_worker : public RcppParallel::Worker {
  const sp_mat& Y;
  const mat& L;
  const mat& F;
  vec& logliks;

  // This is used to create a lik_glmpca_pois_log_sp_worker object.
  lik_glmpca_pois_log_sp_worker (const sp_mat& Y, const mat& L, 
				 const mat& F, vec& logliks) :
    Y(Y), L(L), F(F), logliks(logliks) { };

  // This performs the log-likelihood computations for a given range
  // of columns of Y.
  void operator() (std::size_t begin, std::size_t end) {
    for (unsigned int j = begin; j < end; j++)
      logliks(j) = compute_loglik_glmpca_pois(Y,L,F,j);
  }
};

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
  unsigned int m = Y.n_cols;
  vec logliks(m);
  lik_glmpca_pois_log_sp_worker worker(Y,L,F,logliks);
  parallelFor(0,m,worker);
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
