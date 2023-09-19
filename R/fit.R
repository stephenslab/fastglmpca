#' @title Fit Poisson GLM-PCA Model to Count Data
#' 
#' @description Fit a Poisson GLM-PCA model by maximum-likelihood.
#'   
#' @details In generalized principal component analysis (GLM-PCA)
#' based on a Poisson likelihood, the counts
#' \eqn{y_{ij}} in the n x p matrix \eqn{Y} are modeled as
#' \deqn{y_{ij} \sim Poisson(\lambda_{ij}).} The logarithm of each
#' Poisson rate is defined as a linear combination of the parameters:
#' \deqn{\log \lambda_{ij} = \sum_{k=1}^K d_{k} u_{ik} v_{jk} + 
#' \sum_{m=1}^{n_{x}} x_{im} b_{jm} +
#' \sum_{l=1}^{n_{z}} w_{il} z_{jl} = (UDV' + XB' + WZ')_{ij},} where
#' \eqn{U} is an orthonormal \eqn{n \times K} matrix, \eqn{D} is a diagonal \eqn{K \times K} matrix,
#' \eqn{V} is an orthonormal \eqn{p \times K} matrix, \eqn{X} is a fixed \eqn{n \times n_{x}}
#' matrix of row specific covariates, \eqn{B} is an \eqn{p \times n_{x}} matrix
#' of coefficients of the row specific covariates, \eqn{Z} is a fixed 
#' \eqn{p \times n_{z}} matrix of column specific covariates, and
#' \eqn{W} is an \eqn{n \times n_{z}} matrix of coefficients of column
#' specific covariates. \eqn{K} is a tuning parameter specifying the rank of the matrices
#' \eqn{U} and \eqn{V}. \code{fit_glmpca_pois} computes maximum-likelihood
#' estimates (MLEs) of \eqn{U}, \eqn{D}, \eqn{V}, and if any row or column
#' covariates are present, respectively, \eqn{B} and \eqn{W}.
#' 
#' The algorithm works by repeatedly alternating between updating \eqn{U} and \eqn{W} with
#' \eqn{V} and \eqn{Z} fixed and updating
#' \eqn{V} and \eqn{Z} with \eqn{U} and \eqn{W} fixed. Each update takes the 
#' form of a series of Poisson regressions solved using cyclic co-ordinate
#' descent (CCD). When the algorithm
#' terminates, we rotate the columns of \eqn{U} and \eqn{V} 
#' so that they form an orthonormal set, and then we calculate \eqn{D} appropriately.
#' This rotation does not change the final log-likelihood.
#' 
#' The \code{control} argument is a list in which any of the following
#' named components will override the default optimization algorithm
#' settings (as they are defined by \code{fit_glmpca_pois_control_default}):
#' 
#' \describe{
#'
#' \item{\code{num_ccd_iter}}{Number of ccd updates to be made to parameters
#'   at each iteration of the algorithm.}
#'
#' \item{\code{line_search}}{Boolean indicating if backtracking line search
#'   should be performed at each ccd iteration.}
#'
#' \item{\code{alpha}}{alpha value of line search between 0 and 0.5.}
#'
#' \item{\code{beta}}{beta value of line search between 0 and 0.5.}
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
#' @param Y The n x m matrix of counts; all entries of \code{Y} should
#'   be non-negative. It can be a sparse matrix (class
#'   \code{"dgCMatrix"}) or dense matrix (class \code{"matrix"}).
#'   
#' @param K Integer 1 or greater specifying the rank of the matrix
#'   factorization. This should only be provided if the initial fit
#'   (\code{fit0}) is not.
#'   
#' @param fit0 Initial model fit. It should be an object of class
#'   \dQuote{glmpca_fit_pois}, such as an output from
#'   \code{init_glmpca_pois} or a previous call to
#'   \code{fit_glmpca_pois}.
#'   
#' @param tol The optimization steps when the change in the
#'   log-likelihood between two successive iterations is less than this
#'   amount.
#'   
#' @param min_iter Minimum nummber of updates to be performed.
#' 
#' @param max_iter Maximum number of updates to be performed.
#' 
#' @param verbose Boolean indicating if likelihood should be printed
#' at each step.
#'   
#' @param control List of control parameters to modify behavior of
#' \code{algorithm}.
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
#' @return An object capturing the final state of the model fit. It
#'   will contain the final values of \eqn{U}, \eqn{D} and \eqn{V}, as
#'   well as a list \code{progress} that records information about the
#'   algorithm progress at each iteration.
#'
#' @importFrom utils modifyList
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
    fit0 = init_glmpca_pois(Y,K), 
    tol = 1e-4,
    min_iter = 1,
    max_iter = 100,
    verbose = TRUE,
    control = list()) {

  # Verify and prepare input argument "Y".
  verify.count.matrix(Y)
  n <- nrow(Y)
  m <- ncol(Y)
  if (is.integer(Y))
    storage.mode(Y) <- "double"

  # Check and prepare input arguments "K" and "fit0"
  if (!((missing(K) & !missing(fit0)) |
        (!missing(K) & missing(fit0))))
    stop("Provide either \"K\" or \"fit0\" but not both")
  if (missing(fit0)) {
    if (!(is.scalar(K) & all(K >= 1)))
      stop("Input argument \"K\" should be an integer 1 or greater")
    force(fit0)
  }
  if (!inherits(fit0,"glmpca_pois_fit"))
    stop("Input argument \"fit0\" should be an object of class ",
         "\"glmpca_fit_pois\", such as an output of init_glmpca_pois")
  verify.fit(fit0)
  K <- ncol(fit0$U)
  
  # Check input arguments "tol", "min_iter" and "max_iter".
  if (!(is.scalar(tol) & all(tol > 0)))
    stop("Input argument \"tol\" should be a positive number")
  if (any(min_iter > max_iter))
    stop("\"min_iter\" should be less than or equal to \"max_iter\"")
  
  # Check and process input argument "control".
  control <- modifyList(fit_glmpca_pois_control_default(),control, 
                        keep.null = TRUE)

  # Set up the internal "fit" object.
  fit <- list(LL = t(cbind(fit0$U %*% fit0$D,fit0$X,fit0$W)),
              FF = t(cbind(fit0$V,fit0$B,fit0$Z)),
              fixed_l = numeric(0)
              fixed_f = numeric(0))

  # Determine which rows of LL, FF are "clamped".
  nx <- ifelse(length(fit0$X) > 0,ncol(fit0$X),0)
  nz <- ifelse(length(fit0$Z) > 0,ncol(fit0$Z),0)
  if (nx > 0)
    fit$fixed_l <- c(fit$fixed_l,K + seq(1,nx))
  fit$fixed_l <- c(fit$fixed_l,K + nx + fit0$fixed_w_cols)
  fit$fixed_f <- c(fit$fixed_f,K + fit0$fixed_b_cols)
  if (nz > 0)
    fit$fixed_f <- c(fit$fixed_f,K + nx + seq(1,nz))

  # Perform the updates.
  fit <- fit_glmpca_pois_main_loop(fit,Y,min_iter,max_iter,tol,verbose,control)
  
  # POSTPROCESSIING.
  fit$progress$time <- fit$progress$time[1:t]
  fit$progress$iter <- fit$progress$iter[1:t]
  fit$progress$loglik <- fit$progress$loglik[1:t]
    
  if (control$calc_deriv) {
    fit$progress$max_FF_deriv <- fit$progress$max_FF_deriv[1:t]
    fit$progress$max_LL_deriv <- fit$progress$max_LL_deriv[1:t]
  }
  if (control$calc_max_diff) {
    fit$progress$max_diff_FF <- fit$progress$max_diff_FF[1:t]
    fit$progress$max_diff_LL <- fit$progress$max_diff_LL[1:t] 
  }
  
  colnames(fit$LL) <- rownames(Y)
  colnames(fit$FF) <- colnames(Y)
  
  fit <- postprocess_fit(fit,nx,nz,K)
  return(fit)
}

# This implements the core part of fit_glmpca_pois.
#
#' @importFrom Matrix t
#' @importFrom Matrix tcrossprod
#' @importFrom MatrixExtra tcrossprod
#' @importFrom MatrixExtra mapSparse
fit_glmpca_pois_main_loop <- function (fit, Y, min_iter, max_iter, tol,
                                       verbose, control) {
  n <- nrow(Y)
  m <- ncol(Y)
  K <- nrow(fit$LL)
    
  # Get the rows to update, subtracting 1 for C/C++.
  LL_update_indices_R    <- setdiff(1:K,fit$fixed_l)
  FF_update_indices_R    <- setdiff(1:K,fit$fixed_f)
  joint_update_indices_R <- sort(intersect(LL_update_indices_R,
                                           FF_update_indices_R))
  LL_update_indices      <- LL_update_indices_R - 1
  FF_update_indices      <- FF_update_indices_R - 1

  # These variables are used to compute the log-likelihood below.
  if (inherits(Y,"sparseMatrix")) {
    loglik_const <- sum(mapSparse(Y,lfactorial))
    loglik_func  <- lik_glmpca_pois_log_sp
  } else {
    loglik_const <- sum(lfactorial(Y))
    loglik_func  <- lik_glmpca_pois_log
  } 

  # Set up the data structure for recording the algorithm's progress.
  progress <- data.frame(iter         = 1:max_iter,
                         loglik       = rep(0,max_iter),
                         time         = rep(0,max_iter),
                         max_FF_deriv = rep(as.numeric(NA),max_iter),
                         max_LL_deriv = rep(as.numeric(NA),max_iter),
                         max_diff_FF  = rep(as.numeric(NA),max_iter),
                         max_diff_LL  = rep(as.numeric(NA),max_iter))

  # Set up other data structures used in the calculations below.
  if (control$calc_deriv) {
    LL_mask <- matrix(1,K,n)
    LL_mask[fit$fixed_u_cols,] <- 0
    if (!inherits(Y,"sparseMatrix"))
      LL_mask <- t(LL_mask)
    FF_mask <- matrix(1,K,m)
    FF_mask[fit$fixed_v_cols,] <- 0
    if (!inherits(Y,"sparseMatrix")) 
      FF_mask <- t(FF_mask)
  }
  
  # Some other housekeeping.
  current_lik <- fit0$loglik
  Y_T <- Matrix::t(Y)
  
  converged <- FALSE
  iter <- 1
  if (verbose)
    cat(sprintf("Fitting GLM-PCA model to a %d x %d count matrix.\n",n,m))
  while (!converged && iter <= max_iter) {
    start_time <- proc.time()
    if (control$calc_max_diff) {
      start_iter_LL <- fit$LL
      start_iter_FF <- fit$FF
    }

    browser()
    
    if (length(LL_update_indices) > 0) {
      
      # Orthonormalize rows of FF that are not fixed.
      if (length(joint_update_indices_R) > 1) {
        svd_out <- svd(t(fit$FF[joint_update_indices_R,]))
        fit$FF[joint_update_indices_R,] <- t(svd_out$u)
        fit$LL[joint_update_indices_R,] <-
          diag(svd_out$d) %*% t(svd_out$v) %*% fit$LL[joint_update_indices_R,]
      }

      update_loadings_faster_parallel(
        F_T = t(fit$FF),
        L = fit$LL,
        M = as.matrix(MatrixExtra::tcrossprod(fit$FF[LL_update_indices_R,],Y)),
        update_indices = LL_update_indices,
        num_iter = control$num_ccd_iter,
        line_search = control$line_search,
        alpha = control$alpha,
        beta = control$beta
      )
    }
    
    if (length(FF_update_indices) > 0) {
      if (length(joint_update_indices_R) > 1) {
        
        # Orthonormalize rows of LL that are not fixed.
        #
        # TO DO: Add option to turn this on or off.
        #
        svd_out <- svd(t(fit$LL[joint_update_indices_R,]))
        fit$LL[joint_update_indices_R,] <- t(svd_out$u)
        fit$FF[joint_update_indices_R,] <-
          diag(svd_out$d) %*% t(svd_out$v) %*% fit$FF[joint_update_indices_R,]
      }

      if (length(fit$fixed_v_cols) > 0) {
        update_factors_faster_parallel(
          L_T = t(fit$LL),
          FF = fit$FF,
          M = as.matrix(MatrixExtra::tcrossprod(fit$LL[FF_update_indices_R, ], Y_T)),
          update_indices = FF_update_indices,
          num_iter = control$num_ccd_iter,
          line_search = control$line_search,
          alpha = control$alpha,
          beta = control$beta
        ) 
        
      } else {
        update_factors_faster_parallel(
          L_T = t(fit$LL),
          FF = fit$FF,
          M = as.matrix(MatrixExtra::tcrossprod(fit$LL,Y_T)),
          update_indices = FF_update_indices,
          num_iter = control$num_ccd_iter,
          line_search = control$line_search,
          alpha = control$alpha,
          beta = control$beta
        ) 
      }
    } 
    
    new_lik <- loglik_func(Y,fit$LL,fit$FF,loglik_const)
    # TO DO: Add check for decreasing log-likelihood.
    if (new_lik >= current_lik && iter >= min_iter) {
      rel_improvement <- new_lik - current_lik
      if (rel_improvement < tol) {
        converged <- TRUE
        if (verbose)
          cat(sprintf("Iteration %d: Log-Likelihood = %+0.8e\n",iter,new_lik))
      } 
    }

    end_iter_time <- proc.time()
    time_since_start <-
      as.numeric(difftime(end_iter_time,start_time,units = "secs"))
    current_lik <- new_lik
    fit$progress[iter,"time"] <- time_since_start
    fit$progress[iter,"iter"] <- iter
    fit$progress$loglik[iter] <- current_lik
    
    if (calc_max_diff) {
      fit$progress[iter,"max_diff_LL"] <- max(abs(fit$LL - start_iter_LL))
      fit$progress[iter,"max_diff_FF"] <- max(abs(fit$FF - start_iter_FF))
    }
    
    if(inherits(Y,"sparseMatrix") && calc_deriv) {
      fit$progress$max_FF_deriv[iter] <-
        max(abs((deriv_product(fit$LL, fit$FF) - fit$LL %*% Y) * FF_mask))
      fit$progress$max_LL_deriv[iter] <-
        max(abs((deriv_product(fit$FF, fit$LL) -
                 Matrix::tcrossprod(fit$FF, Y)) * LL_mask))
    } else if (calc_deriv) {
      fit$progress$max_FF_deriv[iter] <- max(abs(crossprod(exp(crossprod(fit$LL, fit$FF)) - Y, t(fit$LL)) * FF_mask))
      fit$progress$max_LL_deriv[iter] <- max(abs(crossprod(exp(crossprod(fit$FF, fit$LL)) - t(Y), t(fit$FF)) * LL_mask))
    }
    iter <- iter + 1
  }
  
  if (verbose)
    cat(sprintf("Iteration %d: Log-Likelihood = %+0.8e\n",
                iter,current_lik))
  if (t > max_iter && !converged) {
    if (verbose)
      cat(sprintf("Iteration %d: Log-Likelihood = %+0.8e\n",t - 1,new_lik))
    warning("Algorithm reached maximum iterations without convergence")
  }
}

#' @rdname fit_glmpca_pois
#'
#' @export
#' 
fit_glmpca_pois_control_default <- function() {
  list(
    alpha = .01,
    beta = .5,
    line_search = TRUE,
    num_iter = 3,
    ccd_iter_tol = 0,
    calc_deriv = FALSE,
    calc_max_diff = FALSE
  )
}
