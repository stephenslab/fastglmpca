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
  fit$D <- out$D
  return(fit)
}

# Given factorization Y = tcrossprod(U,V), return 
# Y = tcrossprod(U1*D,V1) satisfying the following
# constraints:
#
#   crossprod(U1) = I
#   crossprod(V1) = I
#   D = diag(d)
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
    return(list(U = matrix(U),V = matrix(V),D = matrix(d)))
  } else {
    qr1 <- qr(U)
    qr2 <- qr(V)
    out <- svd(tcrossprod(qr.R(qr1),qr.R(qr2)))
    U   <- qr.Q(qr1) %*% out$u
    V   <- qr.Q(qr2) %*% out$v
    return(list(U = U,V = V,D = diag(out$d)))
  }
}

add_dimnames_to_fit <- function (fit, Y) {
  rownames(fit$U) <- rownames(Y)
  rownames(fit$V) <- colnames(Y)
  colnames(fit$U) <- paste("k",1:K,sep = "_")
  colnames(fit$V) <- paste("k",1:K,sep = "_")
  rownames(fit$D) <- paste("k",1:K,sep = "_")
  colnames(fit$D) <- paste("k",1:K,sep = "_")
  if (length(fit$X) > 0) {
    rownames(fit$X) <- rownames(Y)
    rownames(fit$B) <- colnames(Y)
    colnames(fit$B) <- colnames(fit$X)
  }
  if (length(fit$Z) > 0) {
    rownames(fit$Z) <- colnames(Y)
    rownames(fit$W) <- rownames(Y)
    colnames(fit$W) <- colnames(fit$Z)
  }
  return(fit)
}
