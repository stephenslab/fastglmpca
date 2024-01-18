fastglmpca_internal <- new.env(parent = emptyenv())
fastglmpca_internal$main_loop_iter <- 0
fastglmpca_internal$progress <- NA

#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
create_row_intercept <- function (Y)
  return(list(z = rep(1,ncol(Y)),
              w = log(rowSums(Y)/sum(colMeans(Y)))))

#' @importFrom Matrix colMeans
create_col_size_factor <- function (Y)
  return(list(x = rep(1,nrow(Y)),
              b = log(colMeans(Y))))

# Orthonormalizes a GLM-PCA fit object.
orthonormalize_fit <- function (fit) {
  out <- orthonormalize(fit$U,fit$V)
  fit$U <- out$U
  fit$V <- out$V
  fit$d <- out$d
  return(fit)
}

# Given factorization Y = tcrossprod(U,V), return 
# Y = tcrossprod(U1*diag(d),V1) satisfying the following
# constraints:
#
#   crossprod(U1) = I
#   crossprod(V1) = I
#   d = sort(d,decreasing = TRUE)
#   d > 0
#
orthonormalize <- function (U, V) {
  K <- ncol(U)
  if (K == 1) {
    U  <- drop(U)
    V  <- drop(V)
    du <- sqrt(sum(U^2))
    dv <- sqrt(sum(V^2))
    d  <- du*dv
    U  <- U/du
    V  <- V/dv
    return(list(U = matrix(U),V = matrix(V),d = d))
  } else {
    qr1 <- qr(U)
    qr2 <- qr(V)
    out <- svd(tcrossprod(qr.R(qr1),qr.R(qr2)))
    U   <- qr.Q(qr1) %*% out$u
    V   <- qr.Q(qr2) %*% out$v
    return(list(U = U,V = V,d = out$d))
  }
}

