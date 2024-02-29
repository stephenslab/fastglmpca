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
