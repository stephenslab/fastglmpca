lik_glmpca_pois_log <- function(Y, LL, FF, const) {
  
  H <- crossprod(LL, FF)
  lik <- sum(Y * H - exp(H)) - const
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
    link = c("log", "log1p"),
    tol = 1e-4,
    min_iter = 1,
    max_iter = 100,
    verbose = c("likelihood", "none"),
    algorithm = c("ccd", "irls"),
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

  Y_T <- Matrix::t(Y)
  
  # calculate part of log likelihood that doesn't change
  if (link == "log") {
    
    loglik_const <- sum(lfactorial(Y))
    loglik_func <- lik_glmpca_pois_log
    
  } else if (link == "log1p") {
    
    loglik_const <- sum(lfactorial(Y)) - n * p
    loglik_func <- lik_glmpca_pois_log1p
    
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
  
  iter_vec <- numeric(max_iter)
  loglik_vec <- numeric(max_iter)
  time_vec <- numeric(max_iter)
  
  iter_vec[1] <- 0
  loglik_vec[1] <- current_lik
  
  start_time <- Sys.time()
  time_vec[1] <- 0
  
  while (!converged && t <= max_iter) {
    
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
              update_indices = LL_update_indices,
              num_iter = control$num_iter,
              line_search = control$line_search,
              alpha = control$alpha,
              beta = control$beta
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
              update_indices = FF_update_indices,
              num_iter = control$num_iter,
              line_search = control$line_search,
              alpha = control$alpha,
              beta = control$beta
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
    
    new_lik <- do.call(
      loglik_func,
      list(
        Y = Y, LL = fit$LL, FF = fit$FF, const = loglik_const
      )
    )

    if (new_lik >= current_lik && t >= min_iter) {
      
      rel_improvement <- abs((new_lik - current_lik) / current_lik)
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
    time_vec[t + 1] <- time_since_start
    iter_vec[t + 1] <- t
    loglik_vec[t + 1] <- current_lik
    t <- t + 1
    
  }
  
  fit$progress <- data.frame(
    iteration = iter_vec[1:t], 
    loglik = loglik_vec[1:t],
    time_elapsed = time_vec[1:t]
  )
  
  rownames(fit$LL) <- LL_rownames
  rownames(fit$FF) <- FF_rownames
  
  return(fit)
  
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
    num_iter = 5
  )
}

fit_glmpca_irls_control_default <- function() {
  list(
    num_iter = 10
  )
}

fast_fit <- function(Y, fit0, n_iter = 10) {
  
  fit <- fit0
  loglik_const <- sum(lfactorial(Y))
  loglik <- lik_glmpca_pois_log(Y, fit$LL, fit$FF, loglik_const)
  print(loglik)
  
  Y_T <- t(Y)
  
  for (i in 1:n_iter) {
    
    deriv_L_T <- crossprod(exp(crossprod(fit$FF, fit$LL)) - Y_T, t(fit$FF))
    deriv_L_T_2 <- crossprod(exp(crossprod(fit$FF, fit$LL)), t(fit$FF ^ 2))
    
    newton_L <- t(deriv_L_T / deriv_L_T_2)
    
    step <- 1
    proposed_LL <- fit$LL - step * newton_L
    while(lik_glmpca_pois_log(Y, proposed_LL, fit$FF, loglik_const) < loglik) {
      
      step <- step * .5
      proposed_LL <- fit$LL - step * newton_L
      
    }
    
    fit$LL <- proposed_LL
    
    loglik <- lik_glmpca_pois_log(Y, fit$LL, fit$FF, loglik_const)
    
    print(loglik)
    
    deriv_F_T <- crossprod(exp(crossprod(fit$LL, fit$FF)) - Y, t(fit$LL))
    deriv_F_T_2 <- crossprod(exp(crossprod(fit$LL, fit$FF)), t(fit$LL ^ 2))
    
    newton_F <- t(deriv_F_T / deriv_F_T_2)
    
    step <- 1
    proposed_FF <- fit$FF - step * newton_F
    while(lik_glmpca_pois_log(Y, fit$LL, proposed_FF, loglik_const) < loglik) {
      
      step <- step * .5
      proposed_FF <- fit$FF - step * newton_F
      
    }
    
    fit$FF <- proposed_FF
    
    loglik <- lik_glmpca_pois_log(Y, fit$LL, fit$FF, loglik_const)
    
    print(loglik)
    
  }
  
  return(fit)
  
}

