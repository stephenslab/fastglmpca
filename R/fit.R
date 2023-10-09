#' @title Fit Poisson GLM-PCA Model to Count Data
#' 
#' @description Fit a Poisson GLM-PCA model by maximum-likelihood.
#'   
#' @details In generalized principal component analysis (GLM-PCA)
#' based on a Poisson likelihood, the counts \eqn{y_{ij}} stored in an
#' \eqn{n \times m}{n x m} matrix \eqn{Y} are modeled as \deqn{y_{ij}
#' \sim Pois(\lambda_{ij}),} in which the logarithm of each rate
#' parameter \eqn{\lambda_{ij}} is defined as a linear combination of
#' rank-K matrices to be estimated from the data: \deqn{\log
#' \lambda_{ij} = (UDV')_{ij},} where \eqn{U} and \eqn{V} are
#' orthogonal matrices of dimension \eqn{n \times K}{n x K} and \eqn{m
#' \times K}{m x K}, respectively, and \eqn{D} is a diagonal \eqn{K
#' \times K}{K x K} matrix in which the entries along its diagonal are
#' positive and decreasing. \eqn{K} is a tuning parameter specifying
#' the rank of the matrix factorization. This is the same as the
#' low-rank matrix decomposition underlying PCA (that is, the singular
#' value decomposition), but because we are not using a linear
#' (Gaussian) model, this is called \dQuote{generalized PCA} or
#' \dQuote{GLM PCA}.
#'
#' To allow for additional components that may be fixed,
#' \code{fit_glmpca_pois} can also fit the more general model
#' \deqn{\log \lambda_{ij} = (UDV' + XB' + WZ')_{ij},} in which
#' \eqn{X}, \eqn{Z} are fixed matrices of dimension \eqn{n \times
#' n_x}{n x n_x} and \eqn{m \times n_z}{m x n_z}, respectively, and
#' \eqn{B}, \eqn{W} are matrices of dimension \eqn{m \times n_x}{m x
#' n_x} and \eqn{n \times n_z}{n x n_z} to be estimated from the data.
#' 
#' \code{fit_glmpca_pois} computes maximum-likelihood estimates (MLEs)
#' of \eqn{U}, \eqn{V}, \eqn{D}, \eqn{B} and \eqn{W} satistifying the
#' orthogonality constraints for \eqn{U} and \eqn{V} and the
#' additional constraints on \eqn{D} that the entries are positive and
#' decreasing. This is accomplished by iteratively fitting a series of
#' Poisson GLMs, where each of these individual Poissons GLMs is fitted
#' using a fast cyclic co-ordinate descent (CCD) algorithm.
#' 
#' The \code{control} argument is a list in which any of the following
#' named components will override the default optimization algorithm
#' settings (as they are defined by \code{fit_glmpca_pois_control_default}):
#' 
#' \describe{
#'
#' \item{\code{num_ccd_iter}}{Number of co-ordinate descent
#'  updates to be made to parameters at each iteration of 
#'  the algorithm.}
#'
#' \item{\code{line_search}}{If \code{line_search = TRUE}, a
#'   backtracking line search is performed at each iteration of CCD to
#'   guarantee improvement in the objective (the log-likelihood).}
#'
#' \item{\code{alpha}}{alpha parameter for backtracking line search.
#'   (Should be a number between 0 and 0.5, typically a number near
#'   zero.)}
#'
#' \item{\code{beta}}{beta parameter for backtracking line search
#'   controlling the rate at which the step size is decreased.
#'   (Should be a number between 0 and 0.5.)}
#'   
#' \item{\code{calc_deriv}}{If \code{calc_deriv = TRUE}, the maximum
#'   gradient of \eqn{U} and \eqn{V} is calculated and stored after each
#'   update. This may be useful for assessing convergence of the
#'   optimization, though increases overhead.}
#'   
#' \item{\code{calc_max_diff}}{If \code{calc_max_diff = TRUE}, the
#'   largest change in \eqn{U} and \eqn{V} after each update is
#'   calculated and stored. This may be useful for monitoring progress
#'   of the optimization algorithm.}
#' 
#' \item{\code{orthonormalize}}{If \code{orthonormalize = TRUE}, the
#'   matrices \eqn{U} and \eqn{V} are made to be orthogonal after each
#'   update step. This generally improves the speed of convergence
#'   while incurring minimal overhead.}
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
#' @param tol The optimization stops when the change in the
#'   log-likelihood between two successive iterations is less than this
#'   amount.
#'   
#' @param min_iter Minimum number of updates to be performed.
#' 
#' @param max_iter Maximum number of updates to be performed.
#' 
#' @param verbose If \code{verbose = TRUE}, information about the
#'   algorithm's progress is printed after each update.
#'   
#' @param control List of control parameters to modify behavior of
#'   the optimization algorithm; see \dQuote{Details}.
#' 
#' @references
#' Townes, F. W., Hicks, S. C., Aryee, M. J. and Irizarry,
#'   R. A. (2019). Feature selection and dimension reduction for
#'   single-cell RNA-Seq based on a multinomial model. \emph{Genome Biology}
#'   \bold{20}, 295. \url{https://doi.org/10.1186/s13059-019-1861-6}
#'
#'   Collins, M., Dasgupta, S. and Schapire, R. E. (2002). A
#'   generalization of principal components analysis to the exponential
#'   family. In \emph{Advances in Neural Information Processing Systems} 14.
#'
#' @return An object capturing the final state of the model fit. It
#'   contains the final estimates of \eqn{U}, \eqn{D} and \eqn{V}, and
#'   the other parameters (\eqn{X}, \eqn{B}, \eqn{Z}, \eqn{W}), as well
#'   as the log-likelihood achieved (\code{loglik}) and a data frame
#'   \code{progress} storing information about the algorithm's progress
#'   after each update.
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
  if (n != nrow(fit0$U))
    stop("Input \"Y\" should have same number of rows as fit0$U. ",
         "Did fit0 come from a different Y?")
  if (m != nrow(fit0$V))
    stop("Input \"Y\" should have same number of rows as fit0$V. ",
         "Did fit0 come from a different Y?")
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
  fit <- list(LL = t(cbind(fit0$U %*% sqrt(fit0$D),fit0$X,fit0$W)),
              FF = t(cbind(fit0$V %*% sqrt(fit0$D),fit0$B,fit0$Z)),
              fixed_l = numeric(0),
              fixed_f = numeric(0),
              loglik = fit0$loglik)

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
  out <- fit_glmpca_pois_main_loop(fit,Y,min_iter,max_iter,tol,verbose,control)
  fit <- out$fit
  
  # Prepare the final output.
  out$progress$iter <- max(fit0$progress$iter) + out$progress$iter 
  fit <- list(U = t(fit$LL),
              V = t(fit$FF),
              fixed_b_cols = fit0$fixed_b_cols,
              fixed_w_cols = fit0$fixed_w_cols,
              loglik       = fit$loglik,
              progress     = rbind(fit0$progress,out$progress))
  if (nx > 0) {
    fit$X <- fit$U[,K + seq(1,nx),drop = FALSE]
    fit$B <- fit$V[,K + seq(1,nx),drop = FALSE]
  } else {
    fit$X <- numeric(0)
    fit$B <- numeric(0)
  }
  if (nz > 0) {
    fit$Z <- fit$V[,K + nx + seq(1,nz),drop = FALSE]
    fit$W <- fit$U[,K + nx + seq(1,nz),drop = FALSE]
  } else {
    fit$Z <- numeric(0)
    fit$W <- numeric(0)
  }
  fit$U <- fit$U[,1:K,drop = FALSE]
  fit$V <- fit$V[,1:K,drop = FALSE]
  fit <- orthonormalize_fit(fit)
  dimnames(fit$U) <- dimnames(fit0$U)
  dimnames(fit$V) <- dimnames(fit0$V)
  if (length(fit$X) > 0) {
    dimnames(fit$X) <- dimnames(fit0$X)
    dimnames(fit$B) <- dimnames(fit0$B)
  }
  if (length(fit$Z) > 0) {
    dimnames(fit$Z) <- dimnames(fit0$Z)
    dimnames(fit$W) <- dimnames(fit0$W)
  }
  class(fit) <- c("glmpca_pois_fit","list")
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
  update_indices_l <- sort(setdiff(1:K,fit$fixed_l))
  update_indices_f <- sort(setdiff(1:K,fit$fixed_f))

  # These variables are used to compute the log-likelihood below.
  if (inherits(Y,"sparseMatrix")) {
    loglik_const <- sum(mapSparse(Y,lfactorial))
    loglik_func  <- lik_glmpca_pois_log_sp
  } else {
    loglik_const <- sum(lfactorial(Y))
    loglik_func  <- lik_glmpca_pois_log
  } 

  # Set up the data structure for recording the algorithm's progress.
  progress <- data.frame(iter        = 1:max_iter,
                         loglik      = rep(0,max_iter),
                         time        = rep(0,max_iter),
                         max_deriv_f = rep(as.numeric(NA),max_iter),
                         max_deriv_l = rep(as.numeric(NA),max_iter),
                         max_diff_f  = rep(as.numeric(NA),max_iter),
                         max_diff_l  = rep(as.numeric(NA),max_iter))

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
  
  converged <- FALSE
  iter <- 0
  Y_T <- Matrix::t(Y)
  if (verbose)
    cat(sprintf("Fitting GLM-PCA model to %d x %d count matrix.\n",n,m))
  while (!converged & iter < max_iter) {
    iter <- iter + 1
    start_time <- proc.time()
    if (control$calc_max_diff) {
      start_iter_LL <- fit$LL
      start_iter_FF <- fit$FF
    }

    # Perform a single update of L and F.
    fit <- update_glmpca_pois(Y,Y_T,fit,update_indices_l,
                              update_indices_f,control)

    # Update the "progress" data frame.
    new_lik                 <- loglik_func(Y,fit$LL,fit$FF,loglik_const)
    end_iter_time           <- proc.time()
    progress[iter,"loglik"] <- new_lik
    progress[iter,"time"]   <- (end_iter_time - start_time)["elapsed"]
    if (control$calc_max_diff) {
      progress[iter,"max_diff_l"] <- max(abs(fit$LL - start_iter_LL))
      progress[iter,"max_diff_f"] <- max(abs(fit$FF - start_iter_FF))
    }
    if (control$calc_deriv) {
      if(inherits(Y,"sparseMatrix")) {
        progress[iter,"max_deriv_f"] <-
          max(abs((deriv_product(fit$LL,fit$FF) - fit$LL %*% Y) * FF_mask))
        progress[iter,"max_deriv_l"] <-
          max(abs((deriv_product(fit$FF,fit$LL) -
                   Matrix::tcrossprod(fit$FF,Y)) * LL_mask))
      } else {
        progress[iter,"max_deriv_f"] <-
          max(abs(crossprod(exp(crossprod(fit$LL,fit$FF)) - Y,
                            t(fit$LL)) * FF_mask))
        progress[iter,"max_deriv_l"] <-
          max(abs(crossprod(exp(crossprod(fit$FF,fit$LL)) - t(Y),
                            t(fit$FF)) * LL_mask))
      }
    }
    
    # Check whether the stopping criterion is met.
    if (verbose)
      cat(sprintf("Iteration %d: log-likelihood = %+0.12e\n",iter,new_lik))
    if (new_lik < fit$loglik)
      warning("Detected decrease in the log-likelihood")
    if (new_lik >= fit$loglik & iter >= min_iter) {
      if (new_lik - fit$loglik < tol)
        converged <- TRUE
    }
    fit$loglik <- new_lik
  }
  if (!converged)
    warning("Algorithm did not meet convergence criterion")
  return(list(fit = fit,progress = progress[1:iter,]))
}

#' @rdname fit_glmpca_pois
#'
#' @export
#' 
fit_glmpca_pois_control_default <- function()
  list(alpha = 0.01,
       beta = 0.5,
       line_search = TRUE,
       num_ccd_iter = 3,
       ccd_iter_tol = 0,
       calc_deriv = FALSE,
       calc_max_diff = FALSE,
       orthonormalize = TRUE)

# This implements a single update of L and F in the GLM-PCA model.
update_glmpca_pois <- function (Y, Y_T, fit, update_indices_l,
                                update_indices_f, control) {
  n <- nrow(Y)
  m <- ncol(Y)
  K <- nrow(fit$FF)
  k <- sort(intersect(update_indices_l,update_indices_f))
  if (length(update_indices_l) > 0) {
      
    # If requested, orthogonalize rows of FF that are not fixed.
    if (length(k) > 1 & control$orthonormalize) {
      out        <- svd(t(fit$FF[k,,drop = FALSE]))
      fit$FF[k,] <- t(out$u)
      fit$LL[k,] <- diag(out$d) %*% t(out$v) %*% fit$LL[k,,drop = FALSE]
    }

    # Update the LL matrix.
    LLnew  <- matrix(fit$LL,K,n)
    i      <- update_indices_l - 1
    update_factors_faster_parallel(
        L_T = t(fit$FF),
        FF = LLnew,
        M = as.matrix(MatrixExtra::tcrossprod(fit$FF[update_indices_l,,
                                                     drop = FALSE],Y)),
        update_indices = i,
        num_iter = control$num_ccd_iter,
        line_search = control$line_search,
        alpha = control$alpha,
        beta = control$beta)
    fit$LL <- LLnew
  }

  if (length(update_indices_f) > 0) {

    # If requested, orthogonalize rows of LL that are not fixed.
    if (length(k) > 1 & control$orthonormalize) {
      out        <- svd(t(fit$LL[k,,drop = FALSE]))
      fit$LL[k,] <- t(out$u)
      fit$FF[k,] <- diag(out$d) %*% t(out$v) %*% fit$FF[k,,drop = FALSE]
    }

    # Update the FF matrix.
    FFnew <- matrix(fit$FF,K,m)
    i     <- update_indices_f - 1
    update_factors_faster_parallel(
      L_T = t(fit$LL),
      FF = FFnew,
      M = as.matrix(MatrixExtra::tcrossprod(fit$LL[update_indices_f,,
                                                   drop = FALSE],Y_T)),
      update_indices = i,
      num_iter = control$num_ccd_iter,
      line_search = control$line_search,
      alpha = control$alpha,
      beta = control$beta)
    fit$FF <- FFnew
  }

  return(fit)
}
