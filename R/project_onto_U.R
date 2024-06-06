project_onto_U <- function(
    fit0, 
    new_Y, 
    new_Z
) {
  
  n <- nrow(new_Y)
  m <- ncol(new_Y)
  K <- ncol(fit0$U)
  
  new_V <- matrix(
    data = rnorm(m * K, sd = 1e-3),
    nrow = m,
    ncol = K
  )
  
  # add new size factors
  new_B <- matrix(log(Matrix::colMeans(new_Y)))
  colnames(new_B) <- "col_size_factor"
  n_b <- 1
  
  LL <- t(cbind(fit0$U %*% diag(sqrt(fit0$d)),fit0$X,fit0$W))
  FF <- t(cbind(new_V %*% diag(sqrt(fit0$d)),new_B,new_Z))
  
  update_indices_f <- 1:K
  
  FFnew <- matrix(FF,nrow(FF),m)
  # Now, project onto U...
  update_factors_faster_parallel(
    L_T = t(LL),
    FF = FFnew,
    M = as.matrix(LL[update_indices_f,,drop = FALSE] %*% new_Y),
    update_indices = update_indices_f - 1,
    num_iter = 1000,
    line_search = TRUE,
    alpha = .01,
    beta = .25
  )
  
  FF <- FFnew
  
  # now, I need to orthonormalize the entire fit
  fit <- list()
  
  fit$U <- t(LL[1:K,, drop = FALSE])
  fit$V <- t(FF[1:K,, drop = FALSE])
  
  # I need to add back the dimnames here
  fit$V <- rbind(fit$V, fit0$V %*% diag(sqrt(fit0$d)))
  
  fit <- orthonormalize_fit(fit)
  fit$W <- fit0$W
  fit$X <- fit0$X
  fit$Z <- rbind(new_Z, fit0$Z)
  fit$B <- rbind(new_B, fit0$B)
  
  rownames(fit$U) <- rownames(fit0$U)
  colnames(fit$U) <- colnames(fit0$U)
  
  rownames(fit$V) <- c(colnames(new_Y), rownames(fit0$V))
  colnames(fit$V) <- colnames(fit0$V)
  
  rownames(fit$W) <- rownames(fit0$W)
  colnames(fit$W) <- colnames(fit0$W)
  
  rownames(fit$Z) <- c(colnames(new_Y), rownames(fit0$Z))
  colnames(fit$Z) <- colnames(fit0$Z)
  
  rownames(fit$B) <- c(colnames(new_Y), rownames(fit0$B))
  colnames(fit$B) <- colnames(fit0$B)
  
  rownames(fit$X) <- rownames(fit0$X)
  colnames(fit$X) <- colnames(fit0$X)
  
  return(fit)
  
}


