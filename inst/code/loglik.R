# TO DO: Explain here what this function does, and how to use it.
loglik_fastglmpca <- function (fit, Y, e = 1e-8) {
  Y <- as.matrix(Y)
  a <- sum(lfactorial(Y))
  H <- with(fit,
            tcrossprod(U %*% diag(d),V) +
            tcrossprod(X,B) +
            tcrossprod(W,Z))
  return(sum(Y*H - exp(H)) - a)
}

# TO DO: Explain here what this function does, and how to use it.
loglik_glmpca <- function (fit, Y, e = 1e-8) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  a <- sum(lfactorial(Y))
  H <- with(fit,
            tcrossprod(as.matrix(coefX),as.matrix(X)) +
            tcrossprod(rep(1,n),offsets) +
            tcrossprod(as.matrix(loadings),as.matrix(factors)))
  return(sum(Y*H - exp(H)) - a)
}

# TO DO: Explain here what this function does, and how to use it.
loglik_scgbm <- function (fit, Y, e = 1e-8) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  m <- ncol(Y)
  a <- sum(lfactorial(Y))
  H <- with(fit,
            tcrossprod(U %*% diag(D),V) +
            tcrossprod(alpha,rep(1,m)) +
            tcrossprod(rep(1,n),beta))
  return(sum(Y*H - exp(H)) - a)
}
