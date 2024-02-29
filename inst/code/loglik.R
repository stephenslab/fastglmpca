# Compute the log-likelihood for
#
#   y_ij ~ Pois(h_ij)
#
glmpca_pois_loglik <- function (Y, H, const = sum(lfactorial(as.matrix(Y))))
  sum(Y*H - exp(H)) - const

# Compute the log-rates from a fastglmpca output.
logrates_fastglmpca <- function (fit)
  with(fit,
       tcrossprod(U %*% diag(d),V) +
       tcrossprod(X,B) +
       tcrossprod(W,Z))

# Compute the log-rates from a glmpca output, assuming coefZ is NULL.
logrates_glmpca <- function (fit) {
  n <- nrow(fit$loadings)
  return(with(fit,
              tcrossprod(as.matrix(coefX),as.matrix(X)) +
              tcrossprod(rep(1,n),offsets) +
              tcrossprod(as.matrix(loadings),as.matrix(factors))))
}

# Compute thek log-rates from a scGBM output.
logrates_scgbm <- function (fit) {
  n <- nrow(fit$U)
  m <- nrow(fit$V)
  return(with(fit,
              tcrossprod(U %*% diag(D),V) +
              tcrossprod(alpha,rep(1,m)) +
              tcrossprod(rep(1,n),beta)))
}
