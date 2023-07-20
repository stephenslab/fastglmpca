lik_glmpca_pois_log <- function(Y, LL, FF, const) {
  
  H <- crossprod(LL, FF)
  lik <- sum(Y * H - exp(H)) - const
  return(lik)
  
}

lik_glmpca_pois_log_sp <- function(Y, LL, FF, const) {
  
  Y_summary <- summary(Y)
  lik <- big_elementwise_mult_crossprod(
    LL,
    FF,
    Y_summary$x,
    Y_summary$i - 1,
    Y_summary$j - 1,
    nrow(Y_summary)
  ) - big_exp_crossprod(
    LL,
    FF,
    nrow(Y),
    ncol(Y)
  ) - const
  return(lik)
  
}

#' @title Fit Poisson GLM-PCA Model to Count Data
#' 
#' @description Fit a Poisson GLM-PCA model to input matrix \code{Y}
#'   by maximum likelihood.
#'   
#' @details In generalized principal component analysis (GLM-PCA)
#' based on a Poisson likelihood (Townes et al, 2019), the counts
#' \eqn{y_{ij}} in the n x p matrix \eqn{Y} are modeled as
#' \deqn{y_{ij} \sim Poisson(\lambda_{ij}).} The logarithm of each
#' Poisson rate is defined as a linear combination of the parameters:
#' \deqn{\log \lambda_{ij} = \sum_{k=1}^K d_{k} u_{ik} v_{jk} = (UDV')_{ij}.} The model
#' parameters are stored as an n x K matrix \eqn{U} with entries
#' \eqn{u_{ik}}, a p x K matrix \eqn{V} with entries \eqn{v_{jk}}, and a diagonal matrix
#' D with kth diagonal entry \eqn{d_{k}}.
#' \eqn{K} is a tuning parameter specifying the rank of the matrices
#' \eqn{U} and \eqn{V}. \code{fit_glmpca_pois} computes maximum-likelihood
#' estimates (MLEs) of \eqn{U}, \eqn{D}, and \eqn{V}.
#' 
#' The algorithm works by repeatedly alternating between updating \eqn{U} with
#' \eqn{V} fixed and updating \eqn{V} with \eqn{U} fixed. Each update takes the 
#' form of a series of Poisson regressions solved using cyclic co-ordinate
#' descent (ccd). When the algorithm
#' terminates, we rotate the columns of \eqn{U} and \eqn{V} that are not fixed
#' so that they form an orthonormal set, and then we calculate \eqn{D} appropriately.
#' This rotation does not change the final log-likelihood.
#' 
#' The \code{control} argument is a list in which any of the following
#' named components will override the default optimization algorithm
#' settings (as they are defined by \code{fit_glmpca_control_default}):
#' 
#' \describe{
#'
#' \item{\code{num_iter}}{Number of ccd updates to be made to parameters
#'   at each iteration of the algorithm.}
#'
#' \item{\code{line_search}}{Boolean indicating if backtracking line search
#'   should be performed at each ccd iteration.}
#'
#' \item{\code{alpha}}{alpha value of line search between 0 and .5.}
#'
#' \item{\code{beta}}{beta value of line search between 0 and .5.}
#'   
#' \item{\code{calc_deriv}}{boolean indicated if maximum absolute derivatives of
#'  \eqn{U} and \eqn{V} should be calculated at each step of optimization. This may
#'  be useful for monitoring convergence though may have substantial computational
#'  cost for large matrices.}
#'   
#' \item{\code{calc_max_diff}}{boolean indicating if maximum absolute difference 
#' between successive updates of \eqn{U} and \eqn{V} should be calculated and
#' stored. This may be useful for monitoring convergence.}
#'   
#' }
#'
#' @param Y The n x p matrix of counts; all entries of Y should be
#'   non-negative. Y may be a sparse matrix from the \code{Matrix}
#'   package. 
#'   
#' @param K An integer 1 or greater giving the matrix rank. This
#'   argument will be ignored if the initial fit
#'   (\code{fit0}) is provided.
#'   
#' @param fit0 The initial model fit. It should be an object of class
#'   \dQuote{glmpca_fit}, such as an output from \code{init_glmpca_pois}, 
#'   or from a previous call to \code{fit_glmpca}.
#'   
#' @param tol Positive scalar determining relative tolerance for assessing convergence.
#'   Convergence is determined by comparing the log-likelihood at the previous
#'   iteration to the current iteration. 
#'   
#' @param min_iter Minimum number of updates to \eqn{L} and \eqn{F} to be run.
#' 
#' @param max_iter Maximum number of updates to \eqn{L} and \eqn{F} to be run.
#' 
#' @param verbose Boolean indicating if likelihood should be printed at each step.
#'   
#' @param control List of control parameters to modify behavior of \code{algorithm}.
#' 
#' @references
#' Townes, F. W., Hicks, S. C., Aryee, M. J. and Irizarry,
#' R. A. (2019). Feature selection and dimension reduction for
#' single-cell RNA-Seq based on a multinomial model. \emph{Genome Biology}
#' \bold{20}, 295. \url{https://doi.org/10.1186/s13059-019-1861-6}
#'
#' Collins, M., Dasgupta, S. and Schapire, R. E. (2002). A
#' generalization of principal components analysis to the exponential
#' family. In \emph{Advances in Neural Information Processing Systems} 14.
#'
#' @return An object capturing the final state of the model fit. It will contain
#' the final values of \eqn{U}, \eqn{D} and \eqn{V}, as well as a list \code{progress}
#' that records information about the algorithm progress at each iteration.
#'
#' @import Matrix
#' @importFrom utils modifyList
#' @importFrom MatrixExtra mapSparse
#' 
#' @export
#' 
#' @examples 
#' set.seed(1)
#' 
#' n <- 1000
#' p <- 500
#' K <- 1
#' 
#' data <- generate_glmpca_data_pois(n, p, K)
#' # initialize fit
#' fit0 <- init_glmpca_pois(
#'   Y = data$Y, 
#'   K = K,
#'   fit_col_size_factor = TRUE, 
#'   fit_row_intercept = TRUE
#' )
#' 
#' fit <- fit_glmpca_pois(Y = data$Y, fit0 = fit0)
#'
fit_glmpca_pois <- function(
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
    
    fit0 <- init_glmpca_pois(
      Y = Y,
      K = K
    )
    
  } else {
    
    if (!inherits(fit0,"glmpca_pois_fit"))
      stop("Input argument \"fit0\" should be an object of class ",
           "\"glmpca_fit\", such as an output of init_glmpca_pois")
    
    verify.fit(fit0)
    
  }
  
  fit <- list()
  fit$LL <- t(fit0$U %*% fit0$D)
  fit$FF <- t(fit0$V)
  fit$fixed_loadings <- fit0$fixed_loadings
  fit$fixed_factors <- fit0$fixed_factors
  
  # remove initial fit from local scope to preserve memory
  rm(fit0)
  
  LL_rownames <- rownames(fit$LL)
  FF_rownames <- rownames(fit$FF)
  
  K <- nrow(fit$LL)
    
  # get update indices, subtracting 1 for C++ compatibility
  LL_update_indices <- setdiff(1:K, fit$fixed_loadings) - 1
  FF_update_indices <- setdiff(1:K, fit$fixed_factors) - 1
  
  LL_update_indices_R <- LL_update_indices + 1
  FF_update_indices_R <- FF_update_indices + 1
  
  joint_update_indices_R <- intersect(LL_update_indices_R, FF_update_indices_R)

  Y_T <- Matrix::t(Y)
  
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
        "Fitting rank-%d GLM-PCA model to a %d x %d matrix\n with %d fixed factors and %d fixed loadings.\n",
        K, n, p, length(fit$fixed_factors), length(fit$fixed_loadings)
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
    
    if (calc_max_diff) {
      
      start_iter_LL <- fit$LL
      start_iter_FF <- fit$FF
      
    }
      
    if (length(LL_update_indices) > 0) {
      
      # orthonormalize rows of FF that are not fixed
      if (length(joint_update_indices_R) > 1) {
        
        svd_out <- svd(
          t(fit$FF[joint_update_indices_R, ])
        )
        
        fit$FF[joint_update_indices_R, ] <- t(svd_out$u)
        fit$LL[joint_update_indices_R, ] <- diag(svd_out$d) %*% t(svd_out$v) %*% fit$LL[joint_update_indices_R, ]
        
      }

      update_loadings_faster(
        F_T = t(fit$FF),
        L = fit$LL,
        M = as.matrix(MatrixExtra::tcrossprod(fit$FF[LL_update_indices_R, ], Y)),
        update_indices = LL_update_indices,
        n = n,
        num_iter = control$num_iter,
        line_search = control$line_search,
        alpha = control$alpha,
        beta = control$beta
      )

    }
    
    if (length(FF_update_indices) > 0) {
      
      if (length(joint_update_indices_R) > 1) {
        
        # orthonormalize rows of LL that are not fixed
        svd_out <- svd(
          t(fit$LL[joint_update_indices_R, ])
        )
        
        fit$LL[joint_update_indices_R, ] <- t(svd_out$u)
        
        fit$FF[joint_update_indices_R, ] <- diag(svd_out$d) %*% t(svd_out$v) %*% fit$FF[joint_update_indices_R, ]
        
        
      }

      if (length(fit$fixed_factors) > 0) {
        
        new_lik <- update_factors_faster(
          L_T = t(fit$LL),
          FF = fit$FF,
          M = as.matrix(MatrixExtra::tcrossprod(fit$LL[FF_update_indices_R, ], Y_T)),
          update_indices = FF_update_indices,
          p = p,
          num_iter = control$num_iter,
          line_search = control$line_search,
          alpha = control$alpha,
          beta = control$beta
        ) - loglik_const + sum(MatrixExtra::tcrossprod(fit$LL[fit$fixed_factors, ], Y_T) * fit$FF[fit$fixed_factors, ])
        
      } else {
        
        new_lik <- update_factors_faster(
          L_T = t(fit$LL),
          FF = fit$FF,
          M = as.matrix(MatrixExtra::tcrossprod(fit$LL, Y_T)),
          update_indices = FF_update_indices,
          p = p,
          num_iter = control$num_iter,
          line_search = control$line_search,
          alpha = control$alpha,
          beta = control$beta
        ) - loglik_const
        
      }
      
    } else {
      
      new_lik <- do.call(
        loglik_func,
        list(
          Y = Y, LL = fit$LL, FF = fit$FF, const = loglik_const
        )
      )
      
    }

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
      fit$progress$max_diff_FF[t + 1] <- max(abs(fit$FF - start_iter_FF))
      
    }
    
    if(inherits(Y, "sparseMatrix") && calc_deriv) {
      
      fit$progress$max_FF_deriv[t + 1] <- max(abs((deriv_product(fit$LL, fit$FF) - fit$LL %*% Y) * FF_mask))
      fit$progress$max_LL_deriv[t + 1] <- max(abs((deriv_product(fit$FF, fit$LL) - Matrix::tcrossprod(fit$FF, Y)) * LL_mask))
      
    } else if (calc_deriv) {
      
      fit$progress$max_FF_deriv[t + 1] <- max(abs(crossprod(exp(crossprod(fit$LL, fit$FF)) - Y, t(fit$LL)) * FF_mask))
      fit$progress$max_LL_deriv[t + 1] <- max(abs(crossprod(exp(crossprod(fit$FF, fit$LL)) - t(Y), t(fit$FF)) * LL_mask))
      
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
  
  colnames(fit$LL) <- rownames(Y)
  colnames(fit$FF) <- colnames(Y)
  
  fit <- postprocess_fit(fit)
  
  return(fit)
  
}

fit_glmpca_control_default <- function() {
  list(
    alpha = .01,
    beta = .5,
    line_search = TRUE,
    num_iter = 5,
    ccd_iter_tol = 0,
    calc_deriv = FALSE,
    calc_max_diff = FALSE
  )
}
