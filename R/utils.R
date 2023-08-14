# remove element from list by name
safe_remove_elem <- function(l, name) {
  
  if (name %in% names(l)) {
    
    l <- l[names(l) != name]
    
  }
  
  return(l)
  
}


orthonormalize <- function(U, V) {
  
  K <- ncol(U)
  
  if (K == 1) {
    
    d <- sqrt(abs(mean(U[,1])/mean(V[,1])))
    U[,1] <- U[,1] / d
    V[,1] <- V[,1] * d

    return(
      list(
        U = U,
        V = V,
        D = matrix(data = d)
      )
    )
    
  } else {
    
    qr1 <- qr(U)
    qr2 <- qr(V)
    
    svd1 <- svd(tcrossprod(qr.R(qr1), qr.R(qr2)))
    
    U <- qr.Q(qr1) %*% svd1$u
    V <- qr.Q(qr2) %*% svd1$v
    
    return(
      list(
        U = U,
        V = V,
        D = diag(svd1$d)
      )
    )
    
  }
  
}

orthonormalize_fit_qr <- function(fit) {
  
  orthonormed <- orthonormalize(fit$U, fit$V)
  fit$U <- orthonormed$U
  fit$V <- orthonormed$V
  fit$D <- orthonormed$D
  
  return(fit)

}

postprocess_fit <- function(fit, n_x, n_z, K) {
  
  names(fit)[names(fit) == "LL"] <- "U"
  fit$U <- t(fit$U)
  
  names(fit)[names(fit) == "FF"] <- "V"
  fit$V <- t(fit$V)
  
  if (n_x > 0) {
    
    fit$X <- as.matrix(fit$U[,(K + 1):(K + n_x)])
    fit$B <- as.matrix(fit$V[,(K + 1):(K + n_x)])
    
    rownames(fit$X) <- rownames(fit$U)
    rownames(fit$B) <- rownames(fit$B)
    
    colnames(fit$X) <- colnames(fit$U[(K + 1):(K + n_x)])
    colnames(fit$B) <- colnames(fit$V[(K + 1):(K + n_x)])
    
  } else {
    
    fit$X <- numeric(0)
    fit$B <- numeric(0)
    
  }
  
  if (n_z > 0) {
    
    fit$Z <- as.matrix(fit$V[,(K + n_x + 1):(K + n_x + n_z)])
    fit$W <- as.matrix(fit$U[,(K + n_x + 1):(K + n_x + n_z)])
    
    rownames(fit$Z) <- rownames(fit$V)
    rownames(fit$W) <- rownames(fit$U)
    
    colnames(fit$Z) <- colnames(fit$V[(K + n_x + 1):(K + n_x + n_z)])
    colnames(fit$W) <- colnames(fit$U[(K + n_x + 1):(K + n_x + n_z)])
    
  } else {
    
    fit$Z <- numeric(0)
    fit$W <- numeric(0)
    
  }
  
  if (n_x + n_z > 0) {
    
    fit$U <- as.matrix(fit$U[,1:K])
    fit$V <- as.matrix(fit$V[,1:K])
    
  }
  
  if ("max_FF_deriv" %in% names(fit$progress)) {
    
    names(fit$progress)[names(fit$progress) == "max_FF_deriv"] <- "max_deriv_V"
    
  }
  
  if ("max_LL_deriv" %in% names(fit$progress)) {
    
    names(fit$progress)[names(fit$progress) == "max_LL_deriv"] <- "max_deriv_U"
    
  }
  
  if ("max_diff_FF" %in% names(fit$progress)) {
    
    names(fit$progress)[names(fit$progress) == "max_diff_FF"] <- "max_diff_V"
    
  }
  
  if ("max_diff_LL" %in% names(fit$progress)) {
    
    names(fit$progress)[names(fit$progress) == "max_diff_LL"] <- "max_diff_U"
    
  }
  
  class(fit) <- c("glmpca_pois_fit", "list")
  
  fit <- orthonormalize_fit_qr(fit)
  
  rownames(fit$D) <- colnames(fit$U)
  colnames(fit$D) <- colnames(fit$V)
  
  fit <- safe_remove_elem(fit, "fixed_loadings")
  fit <- safe_remove_elem(fit, "fixed_factors")
  
  return(fit)
  
}
