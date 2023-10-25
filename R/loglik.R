# Return the log-likelihood for the GLM-PCA model when Y is a dense
# matrix; i.e., when is.matrix(Y) is TRUE.
lik_glmpca_pois_log <- function(Y, LL, FF, const) {
  H <- crossprod(LL,FF)
  return(sum(Y*H - exp(H)) - const)
}
