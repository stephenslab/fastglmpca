# Here, I want to devise an approximate fitting algorithm
# This will be much the same as my previous code
# However, here I have to be very careful about the regressions
# and the shape of the X matrix

# The code below looks good.
# I need to work on two things
# (1) The likelihood calculation
# (2) I need to ensure that L and F are always 
# in the form that the algorithms need
# (3) Think more about fixed factors / loadings
# and how they might relate to this problem
# In particular orthogonalization is a bit difficult
# with respect to them
# There may also be a computational advantage when calculating
# the fixed linear term.

# I'm also slightly worried about the intercept term here and how
# that might throw things off. 
# I imagine that I can probably infer a special case with an intercept
# which may allow me to do a better job with the intercept.
# it may also make no difference because it has no variance, but I'm not sure

# I should probably first test things on simulated data and then
# go from there. Another option is ignoring the intercept in the approximate
# calculations and then performing likelihood estimation of the intercept
# at the end or intermittently

approx_glmpca_pois <- function(
    Y, 
    K, 
    fit0, 
    tol = 1e-4,
    min_iter = 1,
    max_iter = 100,
    verbose = TRUE,
    control = list()
) {
  
  verify.count.matrix(Y)
  
  if (!is.scalar(tol) || tol <= 0) {
    
    stop("Input argument \"tol\" must be a positive scalar")
    
  }
  
  if (min_iter > max_iter)
    stop("\"min_iter\" must be less than or equal to \"max_iter\"")
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  control <- modifyList(
    fit_glmpca_control_default(), 
    control, 
    keep.null = TRUE
  )
  
  calc_deriv <- control$calc_deriv
  calc_max_diff <- control$calc_max_diff
  
  if (missing(fit0)) {
    
    if (missing(K)) {
      
      stop("must provide either a \"fit0\" or a value of \"K\"")
      
    } else if (!is.scalar(K) || K < 1) {
      
      stop("\"K\" must be an integer greater than or equal to 1")
      
    }
    
    fit <- init_glmpca_pois(
      Y = Y,
      K = K
    )
    
  } else {
    
    if (!inherits(fit0,"glmpca_fit"))
      stop("Input argument \"fit0\" should be an object of class ",
           "\"glmpca_fit\", such as an output of init_glmpca_pois")
    
    verify.fit(fit0)
    fit <- fit0
    
  }
  
  LL_rownames <- rownames(fit$LL)
  FF_rownames <- rownames(fit$FF)
  
  K <- nrow(fit$LL)
  
  # get update indices, subtracting 1 for C++ compatibility
  LL_update_indices <- setdiff(1:K, fit$fixed_loadings) - 1
  FF_update_indices <- setdiff(1:K, fit$fixed_factors) - 1
  
  # calculate part of log likelihood that doesn't change
  if (!inherits(Y, "sparseMatrix")) {
    
    loglik_const <- sum(lfactorial(Y))
    loglik_func <- lik_glmpca_pois_log
    
  } else {
    
    loglik_const <- sum(mapSparse(Y, lfactorial))
    loglik_func <- lik_glmpca_pois_log_sp
    
  } 
  
  current_lik <- do.call(
    loglik_func,
    list(
      Y = Y, LL = fit$LL, FF = fit$FF, const = loglik_const
    )
  )
  
  converged <- FALSE
  
  t <- 1
  
  if (verbose) {
    cat(
      sprintf(
        "Fitting rank-%d GLM-PCA model to a %d x %d matrix.\n", K, n, p
      )
    )
  }
  
  fit$progress <- list()
  fit$progress[["iter"]] <- numeric(max_iter)
  fit$progress[["loglik"]] <- numeric(max_iter)
  fit$progress[["time"]] <- numeric(max_iter)
  
  if (calc_deriv) {
    
    fit$progress[["max_FF_deriv"]] <- numeric(max_iter)
    fit$progress[["max_LL_deriv"]] <- numeric(max_iter)
    
    LL_mask <- matrix(
      data = 1, nrow = nrow(fit$LL), ncol = ncol(fit$LL)
    )
    
    LL_mask[fit$fixed_loadings, ] <- 0
    if(!inherits(Y, "sparseMatrix")) {
      
      LL_mask <- t(LL_mask)
      
    }
    
    FF_mask <- matrix(
      data = 1, nrow = nrow(fit$FF), ncol = ncol(fit$FF)
    )
    
    FF_mask[fit$fixed_factors, ] <- 0
    if(!inherits(Y, "sparseMatrix")) {
      
      FF_mask <- t(FF_mask)
      
    }
    
  }
  
  if (calc_max_diff) {
    
    fit$progress[["max_diff_FF"]] <- numeric(max_iter)
    fit$progress[["max_diff_LL"]] <- numeric(max_iter)
    
  }
  
  fixed_rows <- union(fit$fixed_factors, fit$fixed_loadings)
  
  fit$progress$iter[1] <- 0
  fit$progress$loglik[1] <- current_lik
  fit$progress$time[1] <- 0
  
  if (calc_max_diff) {
    
    fit$progress$max_diff_FF[1] <- 0
    fit$progress$max_diff_LL[1] <- 0
    
  }
  
  if(inherits(Y, "sparseMatrix") && calc_deriv) {
    
    fit$progress$max_FF_deriv[1] <- max(abs((deriv_product(fit$LL, fit$FF) - fit$LL %*% Y) * FF_mask))
    fit$progress$max_LL_deriv[1] <- max(abs((deriv_product(fit$FF, fit$LL) - Matrix::tcrossprod(fit$FF, Y)) * LL_mask))
    
  } else if (calc_deriv) {
    
    fit$progress$max_FF_deriv[1] <- max(abs(crossprod(exp(crossprod(fit$LL, fit$FF)) - Y, t(fit$LL)) * FF_mask))
    fit$progress$max_LL_deriv[1] <- max(abs(crossprod(exp(crossprod(fit$FF, fit$LL)) - Y_T, t(fit$FF)) * LL_mask))
    
  }
  
  while (!converged && t <= max_iter) {
    
    start_time <- Sys.time()
    
    if (verbose) {
      
      cat(
        sprintf(
          "Iteration %d: Log-Likelihood = %+0.8e\n", t-1, current_lik
        )
      )
      
    }
  
    
    start_iter_LL <- fit$LL
    
    if (length(LL_update_indices) > 0) {
      
      # before the update of LL, get FF to be in the correct form
      svd_out <- svd(
        t(fit$FF), K 
      )
      
      fit$FF <- t(svd_out$u)
      fit$LL <- diag(svd_out$d) %*% t(svd_out$v) %*% fit$LL
      row_stds <- matrixStats::rowSds(fit$FF)
      
      fit$FF <- (1 / row_stds) * fit$FF
      fit$LL <- row_stds * fit$LL
      
      update_loadings_approx_nz_mean(
        n,
        p,
        fit$LL,
        rowMeans(fit$FF),
        as.matrix(Matrix::tcrossprod(fit$FF, Y)),
        LL_update_indices,
        num_iter = control$num_iter,
        line_search = control$line_search,
        alpha = control$alpha,
        beta = control$beta
      )
      
    }
    
    if (length(FF_update_indices) > 0) {
      
      # Transform LL so that it is in the correct form
      svd_out <- svd(
        t(fit$LL), K 
      )
      
      fit$LL <- t(svd_out$u)
      
      fit$FF <- diag(svd_out$d) %*% t(svd_out$v) %*% fit$FF
      
      row_stds <- matrixStats::rowSds(fit$LL)
      
      fit$FF <- row_stds * fit$FF
      fit$LL <- (1 / row_stds) * fit$LL
      
      update_factors_approx_nz_mean (
        n,
        p,
        fit$FF,
        rowMeans(fit$LL),
        as.matrix(fit$LL %*% Y),
        FF_update_indices,
        num_iter = control$num_iter,
        line_search = control$line_search,
        alpha = control$alpha,
        beta = control$beta
      )
      
    }
    
    # rescale loadings and factors for numerical stability
    #d <- sqrt(abs(rowMeans(fit$LL)/rowMeans(fit$FF)))
    #d[fixed_rows] <- 1
    #fit$FF <- fit$FF * d
    #fit$LL <- fit$LL / d
    
    new_lik <- do.call(
      loglik_func,
      list(
        Y = Y, LL = fit$LL, FF = fit$FF, const = loglik_const
      )
    )
    
    if (new_lik >= current_lik && t >= min_iter) {
      
      rel_improvement <- new_lik - current_lik
      if (rel_improvement < tol) {
        
        converged <- TRUE
        if (verbose) {
          
          cat(
            sprintf(
              "Iteration %d: Log-Likelihood = %+0.8e\n", t, new_lik
            )
          )
          
        }
        
      } 
      
    } 
    
    end_iter_time <- Sys.time()
    time_since_start <- as.numeric(difftime(end_iter_time, start_time, units = "secs"))
    current_lik <- new_lik
    fit$progress$time[t + 1] <- time_since_start
    fit$progress$iter[t + 1] <- t
    fit$progress$loglik[t + 1] <- current_lik
    
    if (calc_max_diff) {
      
      fit$progress$max_diff_LL[t + 1] <- max(abs(fit$LL - start_iter_LL))
      fit$progress$max_diff_FF[t + 1] <- max(abs(t(fit$FF) - FF_T))
      
    }
    
    if(inherits(Y, "sparseMatrix") && calc_deriv) {
      
      fit$progress$max_FF_deriv[t + 1] <- max(abs((deriv_product(fit$LL, fit$FF) - fit$LL %*% Y) * FF_mask))
      fit$progress$max_LL_deriv[t + 1] <- max(abs((deriv_product(fit$FF, fit$LL) - Matrix::tcrossprod(fit$FF, Y)) * LL_mask))
      
    } else if (calc_deriv) {
      
      fit$progress$max_FF_deriv[t + 1] <- max(abs(crossprod(exp(crossprod(fit$LL, fit$FF)) - Y, t(fit$LL)) * FF_mask))
      fit$progress$max_LL_deriv[t + 1] <- max(abs(crossprod(exp(crossprod(fit$FF, fit$LL)) - Y_T, t(fit$FF)) * LL_mask))
      
    }
    
    t <- t + 1
    
  }
  
  if (t > max_iter && !converged) {
    
    if (verbose) {
      
      cat(
        sprintf(
          "Iteration %d: Log-Likelihood = %+0.8e\n", t - 1, new_lik
        )
      )
      
    }
    
    warning(
      sprintf("Algorithm reached maximum iterations without convergence.")
    )
    
  }
  
  fit$progress$time <- fit$progress$time[1:t]
  fit$progress$iter <- fit$progress$iter[1:t]
  fit$progress$loglik <- fit$progress$loglik[1:t]
  
  if (calc_deriv) {
    
    fit$progress$max_FF_deriv <- fit$progress$max_FF_deriv[1:t]
    fit$progress$max_LL_deriv <- fit$progress$max_LL_deriv[1:t]
    
  }
  
  if (calc_max_diff) {
    
    fit$progress$max_diff_FF <- fit$progress$max_diff_FF[1:t]
    fit$progress$max_diff_LL <- fit$progress$max_diff_LL[1:t] 
    
  }
  
  rownames(fit$LL) <- LL_rownames
  rownames(fit$FF) <- FF_rownames
  
  return(fit)
  
}

