orthonormalize_fit <- function(fit) {
  
  K <- ncol(fit$U)
  
  U_update_indices <- setdiff(1:K, fit$fixed_loadings) 
  V_update_indices <- setdiff(1:K, fit$fixed_factors) 
  
  joint_update_indices <- intersect(U_update_indices, V_update_indices)
  
  if (length(joint_update_indices) > 1) {
    
    svd_1 <- svd(fit$U[, joint_update_indices])
    svd_2 <- svd(fit$V[, joint_update_indices])
    
    fit$U[, joint_update_indices] <- svd_1$u
    D_joint_update_indices <- diag(svd_1$d) %*% svd_1$v %*% t(svd_2$v) %*% t(diag(svd_2$d))
    d <- rep(1, K)
    d[joint_update_indices] <- diag(D_joint_update_indices)
    fit$D <- diag(d)
    fit$V[, joint_update_indices] <- svd_2$u
    
  }
  
  return(fit)
  
}

postprocess_fit <- function(fit) {
  
  names(fit)[names(fit) == "LL"] <- "U"
  fit$U <- t(fit$U)
  
  names(fit)[names(fit) == "FF"] <- "V"
  fit$V <- t(fit$V)
  
  # here, I should do some sort of normalization
  # can leave that for later
  
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
  
  fit <- orthonormalize_fit(fit)
  
  return(fit)
  
}
