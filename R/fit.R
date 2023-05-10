lik_glmpca_pois_log <- function(Y, LL, FF, const) {
  
  H <- crossprod(LL, FF)
  lik <- sum(Y * H - exp(H)) - const
  return(lik)
  
}

lik_glmpca_bin <- function(Y, N, LL, FF, ll_const) {
  
  H <- crossprod(LL, FF)
  sum(Y * H - N * log(1 + exp(H))) + ll_const
  
}

lik_glmpca_bin_missing <- function(Y, N, LL, FF, ll_const) {
  
  H <- crossprod(LL, FF)
  exp_H <- pmin(exp(H), .Machine$double.xmax)
  sum(Y * H - N * log(1 + exp_H)) + ll_const
  
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

get_non_missing_indices <- function(Y) {
  
  non_missing_indices <- apply(!is.na(Y), 2, function(x) which(x) - 1)
  return(as.list(non_missing_indices))
  
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
#' @param use_daarem fill in
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
  
  calc_deriv <- control$calc_deriv
  calc_max_diff <- control$calc_max_diff
  
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
  
  if (use_daarem) {
    
    num_updates <- 1
    
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
      
      if(inherits(Y, "sparseMatrix")) {
        
        LL <- update_loadings_sp(F_T = FF_T,L = LL,Y_T = Y_T,
                              deriv_const_mat = -Matrix::crossprod(FF_T, Y_T),
                              update_indices = LL_update_indices - 1,
                              num_iter = glmpca_control$num_iter,
                              line_search = glmpca_control$line_search,
                              alpha = glmpca_control$alpha,
                              beta = glmpca_control$beta,
                              ccd_iter_tol = glmpca_control$ccd_iter_tol)
        
        FF <- update_factors_sp(L_T = LL_T,FF = FF,Y = Y,
                             deriv_const_mat = -Matrix::crossprod(LL_T, Y),
                             update_indices = FF_update_indices - 1,
                             num_iter = glmpca_control$num_iter,
                             line_search = glmpca_control$line_search,
                             alpha = glmpca_control$alpha,
                             beta = glmpca_control$beta,
                             ccd_iter_tol = glmpca_control$ccd_iter_tol)
        
      } else {
        
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
        
      }
      
      num_updates <<- num_updates + 1
      
      if (calc_max_diff) {
        
        fit$progress$max_diff_FF[num_updates] <<- max(abs(LL - LL0))
        fit$progress$max_diff_LL[num_updates] <<- max(abs(FF - FF0))
        
      }
      
      if(inherits(Y, "sparseMatrix") && calc_deriv) {
        
        fit$progress$max_FF_deriv[num_updates] <<- max(abs((deriv_product(LL, FF) - LL %*% Y) * FF_mask))
        fit$progress$max_LL_deriv[num_updates] <<- max(abs((deriv_product(FF, LL) - Matrix::tcrossprod(FF, Y)) * LL_mask))
        
      } else if(calc_deriv) {
        
        fit$progress$max_FF_deriv[num_updates] <<- max(abs(crossprod(exp(crossprod(LL, FF)) - Y, t(LL)) * FF_mask))
        fit$progress$max_LL_deriv[num_updates] <<- max(abs(crossprod(exp(crossprod(FF, LL)) - Y_T, t(FF)) * LL_mask))
        
      }
      
      return(c(LL[LL_update_indices,],
               FF[FF_update_indices,]))
    }
    
    glmpca_objective <- function (par, Y, Y_T, LL0, FF0,
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
      
      if(!inherits(Y, "sparseMatrix")) {
        
        out <- lik_glmpca_pois_log(Y,LL,FF,loglik_const)
        
      } else {
        
        out <- lik_glmpca_pois_log_sp(Y,LL,FF,loglik_const)
        
      }
      
      cat(sprintf("loglik = %0.6f\n",out))
      return(out)
      

    }
    
    n  <- ncol(fit$LL)
    m  <- ncol(fit$FF)
    k  <- nrow(fit$LL)
    
    out <- daarem::daarem(c(fit$LL[LL_update_indices + 1,],
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
    
    fit$progress$time <- fit$progress$time[1:num_updates]
    fit$progress$iter <- fit$progress$iter[1:num_updates]
    fit$progress$loglik <- out$objfn.track
    
    if (calc_deriv) {
      
      fit$progress$max_FF_deriv <- fit$progress$max_FF_deriv[1:num_updates]
      fit$progress$max_LL_deriv <- fit$progress$max_LL_deriv[1:num_updates]
      
    }
    
    if (calc_max_diff) {
      
      fit$progress$max_diff_FF <- fit$progress$max_diff_FF[1:num_updates]
      fit$progress$max_diff_LL <- fit$progress$max_diff_LL[1:num_updates]
      
    }
    
  } else{
  
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
      start_iter_LL <- fit$LL
        
      if (length(LL_update_indices) > 0) {
        
        if (algorithm == "ccd") {
          
          if (inherits(Y, "sparseMatrix")) {
            
            if (link == "log") {
              
              update_loadings_sp(
                F_T = FF_T,
                L = fit$LL,
                Y_T = Y_T,
                deriv_const_mat = -Matrix::crossprod(FF_T, Y_T),
                update_indices = LL_update_indices,
                num_iter = control$num_iter,
                line_search = control$line_search,
                alpha = control$alpha,
                beta = control$beta,
                ccd_iter_tol = control$ccd_iter_tol
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
              
              new_lik <- update_factors_sp(
                L_T = LL_T,
                FF = fit$FF,
                Y = Y,
                deriv_const_mat = -Matrix::crossprod(LL_T, Y),
                update_indices = FF_update_indices,
                num_iter = control$num_iter,
                line_search = control$line_search,
                alpha = control$alpha,
                beta = control$beta,
                ccd_iter_tol = control$ccd_iter_tol
              ) - loglik_const
              
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
      
      # new_lik <- do.call(
      #   loglik_func,
      #   list(
      #     Y = Y, LL = fit$LL, FF = fit$FF, const = loglik_const
      #   )
      # )
  
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
    
  }
  
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
    alpha = .01,
    beta = .5,
    line_search = TRUE,
    num_iter = 5,
    ccd_iter_tol = 0,
    calc_deriv = FALSE,
    calc_max_diff = FALSE
  )
}

fit_glmpca_irls_control_default <- function() {
  list(
    num_iter = 10
  )
}

#' @title Fit Binomial GLM-PCA Model to Count Data
#' 
#' @description Fit a binomial GLM-PCA model to input matrices \code{Y} and \code{N}
#'   by maximum likelihood.
#'   
#' @details In generalized principal component analysis (GLM-PCA)
#' based on a Binomial likelihood (Townes et al, 2019), the counts
#' \eqn{y_{ij}} in the n x p matrix \eqn{Y} are modeled as
#' \deqn{y_{ij} \sim Binomial(n_{ij}, p_{ij}).} The logit of each
#' binomial probability is defined as a linear combination of the parameters:
#' \deqn{logit p_{ij} = \sum_{k=1}^K l_{ki} f_{kj} = (L'F)_{ij}.} The model
#' parameters are stored as an K x n matrix \eqn{L} with entries
#' \eqn{l_{jk}} and an K x p matrix \eqn{F} with entries \eqn{f_{ik}}.
#' \eqn{K} is a tuning parameter specifying the rank of the matrices
#' \eqn{L} and \eqn{F}. \code{fit_glmpca_binom} computes maximum-likelihood
#' estimates (MLEs) of \eqn{L} and \eqn{F}.
#' 
#' The \code{control} argument is a list in which any of the following
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
#'   non-negative. 
#'   
#' @param N The n x p matrix of trials; all entries of N should be
#'   at least 1.
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
#' @import Matrix
#' @importFrom utils modifyList
#' @importFrom MatrixExtra mapSparse
#' 
#' @export
#' 
#' @examples 
#' set.seed(1)
#' 
#' n <- 300
#' p <- 100
#' K <- 2
#' 
#' LL <- matrix(
#'   data = rnorm(K * n, sd = 1),
#'   nrow = K,
#'   ncol = n
#' )
#' 
#' FF <- matrix(
#'   data = rnorm(K * p, sd = 1),
#'   nrow = K,
#'   ncol = p
#' )
#' 
#' H <- crossprod(LL, FF)
#' 
#' prob <- exp(H) / (1 + exp(H))
#' 
#' N <- matrix(
#'   data = 1 + rpois(n * p, lambda = 5),
#'   nrow = n, 
#'   ncol = p
#' )
#' 
#' Y <- matrix(
#' data = rbinom(
#'     n * p,
#'     size = as.vector(N),
#'     prob = as.vector(prob)
#'   ),
#'   nrow = n,
#'   ncol = p
#' )
#' 
#' 
#' fit0 <- init_glmpca_binom(
#'   Y = Y,
#'   K = 2
#' )
#' 
#' final_fit <- fit_glmpca_binom(
#'  Y = Y, 
#'  N = N,
#'  fit0 = fit0,
#'  max_iter = 25
#' )
#'
fit_glmpca_binom <- function(
    Y, 
    N,
    K, 
    fit0, 
    tol = 1e-4,
    min_iter = 1,
    max_iter = 100,
    control = list() 
  ) {
  
  control <- modifyList(
    fit_glmpca_control_default("ccd"), 
    control, 
    keep.null = TRUE
  )
  
  if (missing(fit0)) {
    
    if (missing(K)) {
      
      stop("must provide either a \"fit0\" or a value of \"K\"")
      
    } else if (!is.scalar(K) || K < 1) {
      
      stop("\"K\" must be an integer greater than or equal to 1")
      
    }
    
    fit <- init_glmpca_binom(
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
  
  LL_rownames <- rownames(fit$LL)
  FF_rownames <- rownames(fit$FF)
  
  K <- nrow(fit$LL)
  
  # get update indices, subtracting 1 for C++ compatibility
  LL_update_indices <- setdiff(1:K, fit$fixed_loadings) - 1
  FF_update_indices <- setdiff(1:K, fit$fixed_factors) - 1
  
  fixed_rows <- union(fit$fixed_factors, fit$fixed_loadings)
  
  Y_T <- t(Y)
  N_T <- t(N)
  
  missing_data <- FALSE
  if (any(is.na(Y)) || any(is.na(N))) {
    
    if (!(all(is.na(Y) == is.na(N)))) {
      
      stop("The missing entries of Y and N must have the same indices")
      
    }
    
    missing_data <- TRUE
    nonmissing_idx_Y <- get_non_missing_indices(Y)
    nonmissing_idx_Y_T <- get_non_missing_indices(Y_T)
    
    Y_lik_version <- Y
    Y_lik_version[is.na(Y_lik_version)] <- 0
    
    N_lik_version <- N
    N_lik_version[is.na(N_lik_version)] <- 0
    
  }
  
  if (missing_data) {
    
    ll_const <- sum(lchoose(na.omit(as.vector(N)), na.omit(as.vector(Y))))
    current_lik <- lik_glmpca_bin_missing(
      Y_lik_version, N_lik_version, fit$LL, fit$FF, ll_const
    )
    
  } else {
    
    ll_const <- sum(lchoose(as.vector(N), as.vector(Y)))
    current_lik <- lik_glmpca_bin(Y, N, fit$LL, fit$FF, ll_const)
    
  }
  
  converged <- FALSE
  
  t <- 1
  
  while (!converged && t <= max_iter) {
      
    cat(
      sprintf(
        "Iteration %d: Log-Likelihood = %+0.8e\n", t-1, current_lik
      )
    )
    
    if (length(LL_update_indices) > 0) {
      
      print("Updating L...")
    
      if (missing_data) {
        
        fit$LL <- update_loadings_missing_bin(
          t(fit$FF),
          fit$LL,
          Y_T,
          N_T,
          nonmissing_idx_Y_T,
          LL_update_indices,
          control$num_iter,
          control$line_search,
          control$alpha,
          control$beta,
          control$ccd_iter_tol
        )
        
      } else {
        
        fit$LL <- update_loadings_bin(
          t(fit$FF),
          fit$LL,
          Y_T,
          N_T,
          LL_update_indices,
          control$num_iter,
          control$line_search,
          control$alpha,
          control$beta,
          control$ccd_iter_tol
        )
        
      }
      
    }
    
    if (length(FF_update_indices) > 0) {
       
      print("Updating F...")
      
      if(t == 32) {
        
        print(0)
        
      }
      
      if (missing_data) {
        
        fit$FF <- update_factors_missing_bin(
          t(fit$LL),
          fit$FF,
          Y,
          N,
          nonmissing_idx_Y,
          FF_update_indices,
          control$num_iter,
          control$line_search,
          control$alpha,
          control$beta,
          control$ccd_iter_tol
        )
        
      } else {
        
        fit$FF <- update_factors_bin(
          t(fit$LL),
          fit$FF,
          Y,
          N,
          FF_update_indices,
          control$num_iter,
          control$line_search,
          control$alpha,
          control$beta,
          control$ccd_iter_tol
        )
        
      }
      
    }
    
    # rescale loadings and factors for numerical stability
    d <- sqrt(abs(rowMeans(fit$LL)/rowMeans(fit$FF)))
    d[fixed_rows] <- 1
    fit$FF <- fit$FF * d
    fit$LL <- fit$LL / d
    
    # now, check for convergence
    if (missing_data) {
      
      new_lik <- lik_glmpca_bin_missing(
        Y_lik_version, N_lik_version, fit$LL, fit$FF, ll_const
      )
      
    } else {
      
      new_lik <- lik_glmpca_bin(Y, N, fit$LL, fit$FF, ll_const)
      
    }
    
    
    if (new_lik >= current_lik && t >= min_iter) {
      
      rel_improvement <- new_lik - current_lik
      if (rel_improvement < tol) {
        
        converged <- TRUE
        cat(
          sprintf(
            "Iteration %d: Log-Likelihood = %+0.8e\n", t, new_lik
          )
        )
        
      } else if (t == max_iter) {
        
        warning("Algorithm hit maximum iterations without convergence")
        cat(
          sprintf(
            "Iteration %d: Log-Likelihood = %+0.8e\n", t, new_lik
          )
        )
        
      } 
      
    }
    
    current_lik <- new_lik
    t <- t + 1

  }
  
  return(fit)
  
}

