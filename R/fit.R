lik_glmpca_pois_log <- function(Y, LL, FF, const) {
  
  H <- crossprod(LL, FF)
  lik <- sum(Y * H - exp(H)) - const
  if (is.infinite(lik)) {
    browser()
  }
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


lik_glmpca_pois_log1p <- function(Y, LL, FF, const) {
  
  exp_H <- exp(crossprod(LL, FF))
  lik <- sum(Y * log(exp_H - 1) - exp_H) - const
  return(lik)
  
}

#' @title Fit GLM-PCA Model to Count Data
#' 
#' @description Fit a GLM-PCA model to input matrix \code{Y}
#'   by maximum likelihood.
#'   
#' @details In generalized principal component analysis (GLM-PCA)
#' based on a Poisson likelihood (Townes et al, 2019), the counts
#' \eqn{y_{ij}} in the n x p matrix \eqn{Y} are modeled as
#' \deqn{y_{ij} \sim Poisson(\lambda_{ij}).} The logarithm of each
#' Poisson rate is defined as a linear combination of the parameters:
#' \deqn{\log \lambda_{ij} = \sum_{k=1}^K l_{ki} f_{kj} = (L'F)_{ij}.} The model
#' parameters are stored as an K x n matrix \eqn{L} with entries
#' \eqn{l_{jk}} and an K x p matrix \eqn{F} with entries \eqn{f_{ik}}.
#' \eqn{K} is a tuning parameter specifying the rank of the matrices
#' \eqn{L} and \eqn{F}. \code{fit_glmpca} computes maximum-likelihood
#' estimates (MLEs) of \eqn{L} and \eqn{F}.
#' 
#' #' The \code{control} argument is a list in which any of the following
#' named components will override the default optimization algorithm
#' settings (as they are defined by \code{fit_glmpca_control_default}):
#' 
#' \describe{
#'
#' \item{\code{num_iter}}{Number of updates to be made to parameters
#'   at each outer iteration of the algorithm.}
#'
#' \item{\code{line_search}}{Boolean indicating if backtracking line search
#'   should be performed. Only used if \code{algorithm} is set to \code{"ccd"}.}
#'
#' \item{\code{alpha}}{alpha value of line search between 0 and .5. 
#'   Only used if \code{algorithm} is set to \code{"ccd"}.}
#'
#' \item{\code{beta}}{beta value of line search between 0 and .5. 
#'   Only used if \code{algorithm} is set to \code{"ccd"}.}}
#'
#' @param Y The n x p matrix of counts; all entries of Y should be
#'   non-negative. Y may be a sparse matrix from the \code{Matrix}
#'   package. 
#'   
#' @param K An integer 1 or greater giving the matrix rank. This
#'   argument will be ignored if the initial fit,
#'   (\code{fit0}), is provided.
#'   
#' @param fit0 The initial model fit. It should be an object of class
#'   \dQuote{glmpca_fit}, such as an output from \code{init_glmpca}, 
#'   or from a previous call to \code{fit_glmpca}.
#'   
#' @param tol Positive scalar determining relative tolerance for assessing convergence.
#'   Convergence is determined by comparing the log-likelihood at the previous
#'   iteration to the current iteration. 
#'   
#' @param min_iter Minimum number of updates to \code{LL} and \code{FF} to be run.
#' 
#' @param max_iter Maximum number of updates to \code{LL} and \code{FF} to be run.
#' 
#' @param verbose String indicating level of printing to be done during model fitting.
#'   \code{"likelihood"} will print only the log-likelihood at each step, where
#'   \code{"detailed"} will print information about the change in parameter values at each
#'   step.
#'   
#' @param algorithm String determining algorithm to use for updating \code{LL} and \code{FF}
#'
#' @param use_daarem Explain what this input argument is for, and how
#'   to use it.
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
#' @return An object capturing the final state of the model fit.
#'
#' @import Matrix
#' @importFrom utils modifyList
#' @importFrom MatrixExtra mapSparse
#' @importFrom daarem daarem
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
#' data <- generate_glmpca_data(n, p, K)
#' # initialize fit
#' fit0 <- init_glmpca(
#'   Y = data$Y, 
#'   K = K,
#'   fit_col_size_factor = TRUE, 
#'   fit_row_intercept = TRUE
#' )
#' 
#' fit <- fit_glmpca(Y = data$Y, fit0 = fit0)
#'
fit_glmpca <- function(
    Y, 
    K, 
    fit0, 
    warmup = FALSE,
    warmup_steps = 5,
    link = c("log", "log1p"),
    tol = 1e-4,
    min_iter = 1,
    max_iter = 100,
    verbose = c("likelihood", "none"),
    algorithm = c("ccd", "irls"),
    use_daarem = FALSE,
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
  
  verbose <- match.arg(verbose)
  
  algorithm <- match.arg(algorithm)
  control <- modifyList(
    fit_glmpca_control_default(algorithm), 
    control, 
    keep.null = TRUE
  )
  
  if (missing(fit0)) {
    
    if (missing(K)) {
      
      stop("must provide either a \"fit0\" or a value of \"K\"")
      
    } else if (!is.scalar(K) || K < 1) {
      
      stop("\"K\" must be an integer greater than or equal to 1")
      
    }
    
    fit <- init_glmpca(
      Y = Y,
      K = K
    )
    
  } else {
    
    if (!inherits(fit0,"glmpca_fit"))
      stop("Input argument \"fit0\" should be an object of class ",
           "\"glmpca_fit\", such as an output of init_glmpca")
    
    verify.fit(fit0)
    fit <- fit0
    
  }
  
  link <- match.arg(link)
  
  LL_rownames <- rownames(fit$LL)
  FF_rownames <- rownames(fit$FF)
  
  K <- nrow(fit$LL)
    
  # get update indices, subtracting 1 for C++ compatibility
  LL_update_indices <- setdiff(1:K, fit$fixed_loadings) - 1
  FF_update_indices <- setdiff(1:K, fit$fixed_factors) - 1
    
  if (algorithm == "irls") {
    
    # here, must construct constant and update vectors
    LL_update_vec <- rep(1, K)
    FF_update_vec <- rep(1, K)
    
    LL_update_vec[fit$fixed_loadings] <- 0
    FF_update_vec[fit$fixed_factors] <- 0
    
    LL_constant_vec <- as.numeric(!LL_update_vec)
    FF_constant_vec <- as.numeric(!FF_update_vec)
    
  }

  Y_T <- t(Y)
  
  # calculate part of log likelihood that doesn't change
  if (link == "log" && !inherits(Y, "sparseMatrix")) {
    
    loglik_const <- sum(lfactorial(Y))
    loglik_func <- lik_glmpca_pois_log
    
  } else if (link == "log1p" && !inherits(Y, "sparseMatrix")) {
    
    loglik_const <- sum(lfactorial(Y)) - n * p
    loglik_func <- lik_glmpca_pois_log1p
    
  } else if (link == "log" && inherits(Y, "sparseMatrix")) {
    
    loglik_const <- sum(mapSparse(Y, lfactorial))
    loglik_func <- lik_glmpca_pois_log_sp
    
  } else if (link == "log1p" && inherits(Y, "sparseMatrix")) {
    
    loglik_const <- sum(mapSparse(Y, lfactorial)) - n * p
    loglik_func <- lik_glmpca_pois_log
    
  }
  
  current_lik <- do.call(
    loglik_func,
    list(
      Y = Y, LL = fit$LL, FF = fit$FF, const = loglik_const
    )
  )

  converged <- FALSE
  
  t <- 1
  
  if (verbose != "none") {
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
  fit$progress[["max_FF_deriv"]] <- numeric(max_iter)
  fit$progress[["max_LL_deriv"]] <- numeric(max_iter)
  
  LL_mask <- matrix(
    data = 1, nrow = nrow(fit$LL), ncol = ncol(fit$LL)
  )
  
  LL_mask[fit$fixed_loadings, ] <- 0
  LL_mask <- t(LL_mask)
  
  FF_mask <- matrix(
    data = 1, nrow = nrow(fit$FF), ncol = ncol(fit$FF)
  )
  
  FF_mask[fit$fixed_factors, ] <- 0
  FF_mask <- t(FF_mask)
  
  fixed_rows <- union(fit$fixed_factors, fit$fixed_loadings)
  
  fit$progress$iter[1] <- 0
  fit$progress$loglik[1] <- current_lik
  fit$progress$time[1] <- 0
  # note, this may use a lot of memory for sparse matrices
  # figure out a way around this
  # could probably do it in a loop
  fit$progress$max_FF_deriv[1] <- max(abs(crossprod(exp(crossprod(fit$LL, fit$FF)) - Y, t(fit$LL)) * FF_mask))
  fit$progress$max_LL_deriv[1] <- max(abs(crossprod(exp(crossprod(fit$FF, fit$LL)) - Y_T, t(fit$FF)) * LL_mask))
  
  # now, run warmup iterations if desired
  if (warmup) {
    
    fit <- warmup(Y, Y_T, fit, loglik_const, warmup_steps, current_lik, verbose)
    t <- t + warmup_steps
    current_lik <- fit$progress$loglik[t]
    
  }

  if (use_daarem) {
    n  <- ncol(fit$LL)
    m  <- ncol(fit$FF)
    k  <- nrow(fit$LL)
    out <- daarem(c(fit$LL[LL_update_indices + 1,],
                    fit$FF[FF_update_indices + 1,]),
                  glmpca_update,glmpca_objective,
                  Y,Y_T,fit$LL,fit$FF,
                  LL_update_indices + 1,FF_update_indices + 1,
                  glmpca_control = control,
                  control = list(maxiter = max_iter,order = 10,tol = 0,
                                 mon.tol = 0.01,kappa = 20,alpha = 1.2))
    LL <- fit$LL
    FF <- fit$FF
    N  <- n*length(LL_update_indices)
    LL[LL_update_indices + 1,] <- out$par[seq(1,N)]
    FF[FF_update_indices + 1,] <- out$par[seq(N+1,length(out$par))]
    fit$LL <- LL
    fit$FF <- FF
    fit$progress <- list(loglik = out$objfn.track)
  } else {
  
  while (!converged && t <= max_iter) {
    
    start_time <- Sys.time()
    
    if (verbose != "none") {
      
      cat(
        sprintf(
          "Iteration %d: Log-Likelihood = %+0.8e\n", t-1, current_lik
        )
      )
      
    }
    
    FF_T <- t(fit$FF)
      
    if (length(LL_update_indices) > 0) {
      
      if (algorithm == "ccd") {
        
        if (inherits(Y, "sparseMatrix")) {
          
          if (link == "log") {
            
            fit$LL <- update_loadings_sp(
              F_T = FF_T,
              L = fit$LL,
              Y_T = Y_T,
              update_indices = LL_update_indices,
              num_iter = control$num_iter,
              line_search = control$line_search,
              alpha = control$alpha,
              beta = control$beta
            )
            
          } else if (link == "log1p") {
            
            fit$LL <- update_loadings_log1p_sp(
              F_T = FF_T,
              L = fit$LL,
              Y_T = Y_T,
              update_indices = LL_update_indices,
              num_iter = control$num_iter,
              line_search = control$line_search,
              alpha = control$alpha,
              beta = control$beta
            )
            
          }
          
        } else {
          
          if (link == "log") {
            
            fit$LL <- update_loadings(
              F_T = FF_T,
              L = fit$LL,
              Y_T = Y_T,
              deriv_const_mat = -crossprod(FF_T, Y_T),
              update_indices = LL_update_indices,
              num_iter = control$num_iter,
              line_search = control$line_search,
              alpha = control$alpha,
              beta = control$beta,
              ccd_iter_tol = control$ccd_iter_tol
            )
            
          } else if (link == "log1p") {
            
            fit$LL <- update_loadings_log1p(
              F_T = FF_T,
              L = fit$LL,
              Y_T = Y_T,
              update_indices = LL_update_indices,
              num_iter = control$num_iter,
              line_search = control$line_search,
              alpha = control$alpha,
              beta = control$beta
            )
            
          }
          
        }
        
      } else if (algorithm == "irls") {
        
        fit$LL <- update_loadings_irls(
          F_T = FF_T,
          L = fit$LL,
          Y_T = Y_T,
          update_vec = LL_update_vec,
          constant_vec = LL_constant_vec,
          num_iter = control$num_iter
        )
        
      }

    }
    
    LL_T <- t(fit$LL)
    
    if (length(FF_update_indices) > 0) {
      
      if (algorithm == "ccd") {
        
        if (inherits(Y, "sparseMatrix")) {
          
          if (link == "log") {
            
            fit$FF <- update_factors_sp(
              L_T = LL_T,
              FF = fit$FF,
              Y = Y,
              update_indices = FF_update_indices,
              num_iter = control$num_iter,
              line_search = control$line_search,
              alpha = control$alpha,
              beta = control$beta
            )
            
          } else if (link == "log1p") {
            
            fit$FF <- update_factors_log1p_sp(
              L_T = LL_T,
              FF = fit$FF,
              Y = Y,
              update_indices = FF_update_indices,
              num_iter = control$num_iter,
              line_search = control$line_search,
              alpha = control$alpha,
              beta = control$beta
            )
            
          }
          
        } else {
          
          if (link == "log") {
            
            fit$FF <- update_factors(
              L_T = LL_T,
              FF = fit$FF,
              Y = Y,
              deriv_const_mat = -crossprod(LL_T, Y),
              update_indices = FF_update_indices,
              num_iter = control$num_iter,
              line_search = control$line_search,
              alpha = control$alpha,
              beta = control$beta,
              ccd_iter_tol = control$ccd_iter_tol
            )
            
          } else if (link == "log1p") {
            
            fit$FF <- update_factors_log1p(
              L_T = LL_T,
              FF = fit$FF,
              Y = Y,
              update_indices = FF_update_indices,
              num_iter = control$num_iter,
              line_search = control$line_search,
              alpha = control$alpha,
              beta = control$beta
            )
            
          }
          
        }
        
      } else if (algorithm == "irls") {
        
        fit$FF <- update_factors_irls(
          L_T = LL_T,
          FF = fit$FF,
          Y = Y,
          update_vec = FF_update_vec,
          constant_vec = FF_constant_vec,
          num_iter = control$num_iter
        )
        
      }

    }
    
    # rescale loadings and factors for numerical stability
    d <- sqrt(abs(rowMeans(fit$LL)/rowMeans(fit$FF)))
    d[fixed_rows] <- 1
    fit$FF <- fit$FF * d
    fit$LL <- fit$LL / d
    
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
        cat(
          sprintf(
            "Iteration %d: Log-Likelihood = %+0.8e\n", t, new_lik
          )
        )
        
      } 
      
    } 
    
    end_iter_time <- Sys.time()
    time_since_start <- as.numeric(difftime(end_iter_time, start_time, units = "secs"))
    current_lik <- new_lik
    fit$progress$time[t + 1] <- time_since_start
    fit$progress$iter[t + 1] <- t
    fit$progress$loglik[t + 1] <- current_lik
    fit$progress$max_FF_deriv[t + 1] <- max(abs(crossprod(exp(crossprod(fit$LL, fit$FF)) - Y, t(fit$LL)) * FF_mask))
    fit$progress$max_LL_deriv[t + 1] <- max(abs(crossprod(exp(crossprod(fit$FF, fit$LL)) - Y_T, t(fit$FF)) * LL_mask))
    t <- t + 1
    
  }
  
  fit$progress$time <- fit$progress$time[1:t]
  fit$progress$iter <- fit$progress$iter[1:t]
  fit$progress$loglik <- fit$progress$loglik[1:t]
  fit$progress$max_FF_deriv <- fit$progress$max_FF_deriv[1:t]
  fit$progress$max_LL_deriv <- fit$progress$max_LL_deriv[1:t]
  
  rownames(fit$LL) <- LL_rownames
  rownames(fit$FF) <- FF_rownames
  }
  
  return(fit)
  
}

glmpca_update <- function (par, Y, Y_T, LL0, FF0,
                           LL_update_indices, FF_update_indices,
                           glmpca_control) {
  n  <- ncol(LL0)
  m  <- ncol(FF0)
  k  <- nrow(LL0)
  LL <- LL0
  FF <- FF0
  N  <- n*length(LL_update_indices)
  LL[LL_update_indices,] <- par[seq(1,N)]
  FF[FF_update_indices,] <- par[seq(N+1,length(par))]
  FF_T <- t(FF)
  LL_T <- t(LL)
  
  LL <- update_loadings(F_T = FF_T,L = LL,Y_T = Y_T,
                        deriv_const_mat = -crossprod(FF_T, Y_T),
                        update_indices = LL_update_indices - 1,
                        num_iter = glmpca_control$num_iter,
                        line_search = glmpca_control$line_search,
                        alpha = glmpca_control$alpha,
                        beta = glmpca_control$beta,
                        ccd_iter_tol = glmpca_control$ccd_iter_tol)
 
  FF <- update_factors(L_T = LL_T,FF = FF,Y = Y,
                       deriv_const_mat = -crossprod(LL_T, Y),
                       update_indices = FF_update_indices - 1,
                       num_iter = glmpca_control$num_iter,
                       line_search = glmpca_control$line_search,
                       alpha = glmpca_control$alpha,
                       beta = glmpca_control$beta,
                       ccd_iter_tol = glmpca_control$ccd_iter_tol)

  return(c(LL[LL_update_indices,],
           FF[FF_update_indices,]))
}

glmpca_objective <- function (par, Y, Y_T, LL0, FF0,
                              LL_update_indices, FF_update_indices,
                              glmpca_control) {
  loglik_const <- sum(lfactorial(Y))
  n  <- ncol(LL0)
  m  <- ncol(FF0)
  k  <- nrow(LL0)
  LL <- LL0
  FF <- FF0
  N  <- n*length(LL_update_indices)
  LL[LL_update_indices,] <- par[seq(1,N)]
  FF[FF_update_indices,] <- par[seq(N+1,length(par))]
  FF_T <- t(FF)
  LL_T <- t(LL)
  out <- lik_glmpca_pois_log(Y,LL,FF,loglik_const)
  cat(sprintf("loglik = %0.6f\n",out))
  return(out)
}

fit_glmpca_control_default <- function(algorithm) {
  
  if (algorithm == "ccd") {
    
    fit_glmpca_ccd_control_default()
    
  } else if (algorithm == "irls") {
    
    fit_glmpca_irls_control_default()
    
  }
  
}

fit_glmpca_ccd_control_default <- function() {
  list(
    alpha = .25,
    beta = .5,
    line_search = TRUE,
    num_iter = 5,
    ccd_iter_tol = 0
  )
}

fit_glmpca_irls_control_default <- function() {
  list(
    num_iter = 10
  )
}

#' @importFrom Matrix crossprod
warmup <- function(Y, Y_T, fit, loglik_const, n_iter, starting_loglik, verbose) {
  
  # I want to fix a mask for LL and FF
  LL_mask <- matrix(
    data = 1, nrow = nrow(fit$LL), ncol = ncol(fit$LL)
  )
  
  LL_mask[fit$fixed_loadings, ] <- 0
  
  FF_mask <- matrix(
    data = 1, nrow = nrow(fit$FF), ncol = ncol(fit$FF)
  )
  
  FF_mask[fit$fixed_factors, ] <- 0
  
  loglik <- starting_loglik
  
  for (i in 1:n_iter) {
    
    if (verbose != "none") {
      
      cat(
        sprintf(
          "Iteration %d: Log-Likelihood = %+0.8e\n", i - 1, loglik
        )
      )
      
    }
    
    start_time <- Sys.time()
    
    deriv_L_T <- crossprod(exp(crossprod(fit$FF, fit$LL)) - Y_T, t(fit$FF))
    deriv_L_T_2 <- crossprod(exp(crossprod(fit$FF, fit$LL)), t(fit$FF ^ 2))
    
    newton_L <- t(deriv_L_T / deriv_L_T_2) * LL_mask
    
    step <- 1
    converged <- FALSE
    while (!converged) {
      
      proposed_LL <- fit$LL - step * newton_L
      new_loglik <- lik_glmpca_pois_log(Y, proposed_LL, fit$FF, loglik_const)
      if (new_loglik >= loglik) {
        
        converged <- TRUE
        
      } else {
        
        step <- .5 * step
        
      }
      
    }
    
    fit$LL <- proposed_LL
    
    loglik <- new_loglik
    
    deriv_F_T <- crossprod(exp(crossprod(fit$LL, fit$FF)) - Y, t(fit$LL))
    deriv_F_T_2 <- crossprod(exp(crossprod(fit$LL, fit$FF)), t(fit$LL ^ 2))
    
    newton_F <- t(deriv_F_T / deriv_F_T_2) * FF_mask
    
    step <- 1
    converged <- FALSE
    while(!converged) {
      
      proposed_FF <- fit$FF - step * newton_F
      new_loglik <- lik_glmpca_pois_log(Y, fit$LL, proposed_FF, loglik_const)
      
      if (new_loglik >= loglik) {
        
        converged <- TRUE
        
      } else {
        
        step <- .5 * step
        
      }
      
    }
    
    fit$FF <- proposed_FF
    
    loglik <- new_loglik
    
    end_time <- Sys.time()
    time_elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    fit$progress$iter[i + 1] <- i
    fit$progress$loglik[i + 1] <- loglik
    fit$progress$time[i + 1] <- time_elapsed
    
  }
  
  return(fit)
  
}

