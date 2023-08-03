orthonormalize_fit_qr <- function(fit) {
  
  K <- ncol(fit$U)
  
  U_update_indices <- setdiff(1:K, fit$fixed_loadings) 
  V_update_indices <- setdiff(1:K, fit$fixed_factors) 
  
  joint_update_indices <- intersect(U_update_indices, V_update_indices)
  
  d <- rep(1, K)
  
  if (length(joint_update_indices) > 1) {
    
    qr1 <- qr(fit$U[, joint_update_indices])
    qr2 <- qr(fit$V[, joint_update_indices])
    
    svd1 <- svd(tcrossprod(qr.R(qr1), qr.R(qr2)))
    
    fit$U[, joint_update_indices] <- qr.Q(qr1) %*% svd1$u
    d[joint_update_indices] <- svd1$d
    fit$V[, joint_update_indices] <- qr.Q(qr2) %*% svd1$v
    
  } 
  
  fit$D <- diag(d)
  
  return(fit)
  
}

postprocess_fit <- function(fit) {
  
  names(fit)[names(fit) == "LL"] <- "U"
  fit$U <- t(fit$U)
  
  names(fit)[names(fit) == "FF"] <- "V"
  fit$V <- t(fit$V)
  
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
  
  return(fit)
  
}
