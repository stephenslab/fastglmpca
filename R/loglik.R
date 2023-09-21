# Return the log-likelihood for the GLM-PCA model when Y is a dense
# matrix; i.e., when is.matrix(Y) is TRUE.
lik_glmpca_pois_log <- function(Y, LL, FF, const) {
  H <- crossprod(LL,FF)
  return(sum(Y*H - exp(H)) - const)
}

# Should return the same result as
#
#   lik_glmpca_pois_log(as.matrix(Y),LL,FF,const)
#
lik_glmpca_pois_log_sp <- function (Y, LL, FF, const) {
  out <- Matrix::summary(Y)
  return(
    big_elementwise_mult_crossprod(LL,FF,out$x,out$i-1,out$j-1,nrow(out))
  - big_exp_crossprod(LL,FF,nrow(Y),ncol(Y))
  - const)
}

