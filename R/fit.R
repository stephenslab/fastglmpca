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
#' using a fast \dQuote{cyclic co-ordinate descent} (CCD) algorithm.
#'
#' The \code{control} argument is a list in which any of the following
#' named components will override the default optimization algorithm
#' settings (as they are defined by
#' \code{fit_glmpca_pois_control_default}). Additional control
#' arguments not listed here can be used to control the behaviour of
#' \code{\link[daarem]{fpiter}} or \code{\link[daarem]{daarem}}; see
#' the help accompanying these functions for details.
#'
#' \describe{
#'
#' \item{\code{use_daarem}}{If \code{use_daarem = TRUE}, the updates
#'   are accelerated using DAAREM; see \code{\link[daarem]{daarem}} for
#'   details.}
#'
#' \item{\code{tol}}{This is the value of the \dQuote{tol} control
#'   argument for \code{\link[daarem]{fpiter}} or
#'   \code{\link[daarem]{daarem}} that determines when to stop the
#'   optimization. In brief, the optimization stops when the change in
#'   the estimates or in the log-likelihood between two successive
#'   updates is less than \dQuote{tol}.}
#'
#' \item{\code{maxiter}}{This is the value of the \dQuote{maxiter}
#'   control argument for \code{\link[daarem]{fpiter}} or
#'   \code{\link[daarem]{daarem}}. In brief, it sets the upper limit on
#'   the number of CCD updates.}
#'
#' \item{\code{convtype}}{This is the value of the \dQuote{convtype}
#'   control argument for \code{\link[daarem]{daarem}}. It determines
#'   whether the stopping criterion is based on the change in the
#'   estimates or the change in the log-likelihood between two
#'   successive updates.}
#'
#' \item{\code{mon.tol}}{This is the value of the \dQuote{mon.tol}
#'   control argument for \code{\link[daarem]{daarem}}. This setting
#'   determines to what extent the monotonicity condition can be
#'   violated.}
#'   
#' \item{\code{training_frac}}{Fraction of the columns of input data \code{Y}
#'   to fit initial model on. If set to \code{1} (default), the model is fit
#'   by optimizing the parameters on the entire dataset. If set between \code{0}
#'   and \code{1}, the model is optimized by first fitting a model on a randomly
#'   selected fraction of the columns of \code{Y}, and then projecting the 
#'   remaining columns of \code{Y} onto the solution. Setting this to a smaller
#'   value will increase speed but decrease accuracy.
#' }
#' 
#' \item{\code{num_projection_ccd_iter}}{Number of co-ordinate descent updates
#'   be made to elements of \code{V} if and when a subset of \code{Y} is 
#'   projected onto \code{U}. Only used if \code{training_frac} is less than 
#'   \code{1}.
#' }
#'
#' \item{\code{num_ccd_iter}}{Number of co-ordinate descent updates to
#'   be made to parameters at each iteration of the algorithm.}
#'
#' \item{\code{line_search}}{If \code{line_search = TRUE}, a
#'   backtracking line search is performed at each iteration of CCD to
#'   guarantee improvement in the objective (the log-likelihood).}
#'
#' \item{\code{ls_alpha}}{alpha parameter for backtracking line search.
#'   (Should be a number between 0 and 0.5, typically a number near
#'   zero.)}
#'
#' \item{\code{ls_beta}}{beta parameter for backtracking line search
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
#'   update step. This improves the speed of convergence without the
#'   DAAREM acceleration; however, should not be used when
#'   \code{use_daarem = TRUE}.}}
#'
#' You may use function \code{\link{set_fastglmpca_threads}} to adjust
#' the number of threads used in performing the updates.
#'
#' @param Y The n x m matrix of counts; all entries of \code{Y} should
#'   be non-negative. It can be a sparse matrix (class
#'   \code{"dsparseMatrix"}) or dense matrix (class \code{"matrix"}).
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
#' @param verbose If \code{verbose = TRUE}, information about the
#'   algorithm's progress is printed after each update.
#'
#' @param control List of control parameters to modify behavior of
#'   the optimization algorithm; see \dQuote{Details}.
#'
#' @references
#'   Townes, F. W., Hicks, S. C., Aryee, M. J. and Irizarry,
#'   R. A. (2019). Feature selection and dimension reduction for
#'   single-cell RNA-Seq based on a multinomial model. \emph{Genome Biology}
#'   \bold{20}, 295. \doi{10.1186/s13059-019-1861-6}
#'
#'   Collins, M., Dasgupta, S. and Schapire, R. E. (2002). A
#'   generalization of principal components analysis to the exponential
#'   family. In \emph{Advances in Neural Information Processing Systems} 14.
#'
#' @return An object capturing the state of the model fit. It contains
#'   estimates of \eqn{U}, \eqn{V} and \eqn{D} (stored as matrices
#'   \code{U}, \code{V} and a vector of diagonal entries \code{d},
#'   analogous to the \code{\link{svd}} return value); the other
#'   parameters (\eqn{X}, \eqn{B}, \eqn{Z}, \eqn{W}); the log-likelihood
#'   achieved (\code{loglik}); information about which columns of
#'   \eqn{B} and \eqn{W} are fixed (\code{fixed_b_cols},
#'   \code{fixed_w_cols}); and a data frame \code{progress} storing
#'   information about the algorithm's progress after each update.
#'
#' @importFrom utils modifyList
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' p <- 100
#' K <- 3
#' dat  <- generate_glmpca_data_pois(n,p,K)
#' fit0 <- init_glmpca_pois(dat$Y,K)
#' fit  <- fit_glmpca_pois(dat$Y,fit0 = fit0)
#'
fit_glmpca_pois <- function(
    Y,
    K,
    fit0 = init_glmpca_pois(Y,K),
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

  # Check and process input argument "control".
  control <- modifyList(fit_glmpca_pois_control_default(),
                        control,keep.null = TRUE)
  
  # Set up the internal fit.
  D <- sqrt(fit0$d)
  if (K == 1)
    D <- matrix(D)
  else
    D <- diag(D)
  LL <- t(cbind(fit0$U %*% D,fit0$X,fit0$W))
  FF <- t(cbind(fit0$V %*% D,fit0$B,fit0$Z))
  
  # Determine which rows of LL and FF are "clamped".
  fixed_l <- numeric(0)
  fixed_f <- numeric(0)
  nx <- ifelse(length(fit0$X) > 0,ncol(fit0$X),0)
  nz <- ifelse(length(fit0$Z) > 0,ncol(fit0$Z),0)
  if (nx > 0)
    fixed_l <- c(fixed_l,K + seq(1,nx))
  fixed_l <- c(fixed_l,K + nx + fit0$fixed_w_cols)
  fixed_f <- c(fixed_f,K + fit0$fixed_b_cols)
  if (nz > 0)
    fixed_f <- c(fixed_f,K + nx + seq(1,nz))
  
  if (control$training_frac == 1) {
    
    # Perform the updates.
    res <- fit_glmpca_pois_main_loop(LL,FF,Y,fixed_l,fixed_f,verbose,control)
    
  } else {
    
    if (control$training_frac <= 0 || control$training_frac > 1)
      stop("control argument \"training_frac\" should be between 0 and 1")
    
    train_idx <- sample(
      1:ncol(Y), 
      size = ceiling(ncol(Y) * control$training_frac)
    )
    
    Y_train <- Y[, train_idx]
    
    if (any(Matrix::rowSums(Y_train) == 0) || any(Matrix::colSums(Y_train) == 0)) {
      
      stop(
        "After subsetting, the remaining values of \"Y\" ",
        "contain a row or a column where all counts are 0. This can cause ",
        "problems with optimization. Please either remove rows / columns ",
        "with few non-zero counts from \"Y\", or set \"training_frac\" to ",
        "a larger value."
        )
      
    }
    
    FF_train <- FF[, train_idx, drop = FALSE]
    FF_test <- FF[, -train_idx, drop = FALSE]
    Y_test <- Y[, -train_idx, drop = FALSE]
    
    test_idx <- 1:ncol(Y)
    test_idx <- test_idx[-train_idx]
    
    # Perform the updates.
    res <- fit_glmpca_pois_main_loop(
      LL,
      FF_train,
      Y_train,
      fixed_l,
      fixed_f,
      verbose,
      control
    )
    
    update_indices_f <- sort(setdiff(1:K,fixed_f))
    
    # now, I just need to project the results back
    update_factors_faster_parallel(
      L_T = t(res$fit$LL),
      FF = FF_test,
      M = as.matrix(res$fit$LL[update_indices_f,,drop = FALSE] %*% Y_test),
      update_indices = update_indices_f - 1,
      num_iter = control$num_projection_ccd_iter,
      line_search = control$line_search,
      alpha = control$ls_alpha,
      beta = control$ls_beta
    )
    
    # now, I need to reconstruct FF, and hopefully compute the log-likelihood
    FF[, train_idx] <- res$fit$FF
    FF[, test_idx] <- FF_test
    res$fit$FF <- FF
    
    if (inherits(Y,"sparseMatrix")) {
      test_loglik_const <- sum(lfactorial(Y_test@x))
      loglik_func  <- lik_glmpca_pois_log_sp
    } else {
      test_loglik_const <- sum(lfactorial(Y_test))
      loglik_func  <- lik_glmpca_pois_log
    }
    
    test_loglik <- loglik_func(Y_test,res$fit$LL,FF_test,test_loglik_const)
    res$loglik <- res$loglik + test_loglik
    
  }
  
  # Prepare the final output.
  res$progress$iter <- max(fit0$progress$iter) + res$progress$iter
  fit <- list(U            = t(res$fit$LL),
              V            = t(res$fit$FF),
              fixed_b_cols = fit0$fixed_b_cols,
              fixed_w_cols = fit0$fixed_w_cols,
              loglik       = res$loglik,
              progress     = rbind(fit0$progress,res$progress))
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
  names(fit$d)    <- paste("k",1:K,sep = "_")
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
#' @importFrom daarem fpiter
#' @importFrom daarem daarem
fit_glmpca_pois_main_loop <- function (LL, FF, Y, fixed_l, fixed_f,
                                       verbose, control) {
  n <- nrow(Y)
  m <- ncol(Y)
  K <- nrow(LL)

  # Get the rows of LL and FF to update.
  update_indices_l <- sort(setdiff(1:K,fixed_l))
  update_indices_f <- sort(setdiff(1:K,fixed_f))

  # These variables are used to compute the log-likelihood below.
  if (inherits(Y,"sparseMatrix")) {
    loglik_const <- sum(lfactorial(Y@x))
    loglik_func  <- lik_glmpca_pois_log_sp
  } else {
    loglik_const <- sum(lfactorial(Y))
    loglik_func  <- lik_glmpca_pois_log
  }

  # Set up the data structure for recording the algorithm's progress.
  fastglmpca_internal$progress <-
    data.frame(iter        = 1:control$maxiter,
               loglik      = rep(0,control$maxiter),
               time        = rep(0,control$maxiter),
               max_deriv_f = rep(as.numeric(NA),control$maxiter),
               max_deriv_l = rep(as.numeric(NA),control$maxiter),
               max_diff_f  = rep(as.numeric(NA),control$maxiter),
               max_diff_l  = rep(as.numeric(NA),control$maxiter))

  # Set up other data structures used in the calculations below.
  LL_mask <- matrix(1,K,n)
  if (!inherits(Y,"sparseMatrix"))
    LL_mask <- t(LL_mask)
  FF_mask <- matrix(1,K,m)
  if (!inherits(Y,"sparseMatrix"))
    FF_mask <- t(FF_mask)

  # Perform the updates using fpiter or daarem.
  fastglmpca_internal$main_loop_iter <- 0
  Y_T <- Matrix::t(Y)
  if (verbose)
    cat(sprintf("Fitting GLM-PCA model to %d x %d count matrix.\n",n,m))
  if (control$use_daarem) {
    updater <- daarem
    control_settings <- c("maxiter","order","tol","mon.tol","cycl.mon.tol",
                          "alpha","kappa","resid.tol","convtype")
  } else {
    updater <- fpiter
    control_settings <- c("tol","maxiter","trace")
  }
  control_settings <- intersect(control_settings,names(control))
  control_daarem <- control[control_settings]
  res <- updater(# These are the inputs needed to run fpiter or daarem.
                 par = fit2par(list(LL = LL,FF = FF),
                               update_indices_l,
                               update_indices_f),
                 fixptfn = fpiter_update,
                 objfn   = fpiter_objective,
                 control = control_daarem,

                 # These arguments are passed along to
                 # fpiter_objective and fpiter_update.
                 control_glmpca_pois = control,
                 LL = LL,FF = FF,
                 LL_mask = LL_mask,FF_mask = FF_mask,
                 Y = Y,Y_T = Y_T,
                 update_indices_l = update_indices_l,
                 update_indices_f = update_indices_f,
                 loglik_func = loglik_func,
                 loglik_const = loglik_const,
                 verbose = verbose)
  if (!res$convergence)
    message(sprintf(paste("Note that the specified convergence criteria was not",
                          "met within %d iterations."),
                    control$maxiter))

  # Prepare the output.
  return(list(fit = par2fit(res$par,LL,FF,update_indices_l,update_indices_f),
              progress = fastglmpca_internal$progress[
                1:fastglmpca_internal$main_loop_iter,],
              loglik = res$value.objfn))
}

#' @rdname fit_glmpca_pois
#'
#' @export
#'
fit_glmpca_pois_control_default <- function()
  list(use_daarem = FALSE,
       maxiter = 100,
       tol = 1e-4,
       training_frac = 1,
       num_projection_ccd_iter = 10,
       mon.tol = 0.05,
       convtype = "objfn",
       line_search = TRUE,
       ls_alpha = 0.001,
       ls_beta = 0.5,
       num_ccd_iter = 3,
       ccd_iter_tol = 0,
       calc_deriv = FALSE,
       calc_max_diff = FALSE,
       orthonormalize = TRUE)

# This implements "objfn" in fpiter or daarem.
fpiter_objective <- function (par, LL, FF, LL_mask, FF_mask, Y, Y_T,
                              update_indices_l, update_indices_f,
                              loglik_func, loglik_const,
                              control_glmpca_pois, verbose) {
  fit <- par2fit(par,LL,FF,update_indices_l,update_indices_f)
  return(loglik_func(Y,fit$LL,fit$FF,loglik_const))
}

# This implements "fixptfn" in fpiter or daarem.
#
# @importFrom Matrix tcrossprod
fpiter_update <- function (par, LL, FF, LL_mask, FF_mask, Y, Y_T,
                           update_indices_l, update_indices_f,
                           loglik_func, loglik_const,
                           control_glmpca_pois, verbose) {
  fastglmpca_internal$main_loop_iter <- fastglmpca_internal$main_loop_iter + 1
  start_time <- proc.time()

  # Set up the internal "fit" object.
  fit <- par2fit(par,LL,FF,update_indices_l,update_indices_f)

  # Keep track of the current estimates.
  fit0 <- fit

  # Perform a single update of LL and FF.
  fit <- update_glmpca_pois(fit$LL,fit$FF,Y,Y_T,update_indices_l,
                            update_indices_f,control_glmpca_pois)

  # Update the "progress" data frame.
  new_lik <- loglik_func(Y,fit$LL,fit$FF,loglik_const)
  fastglmpca_internal$progress[fastglmpca_internal$main_loop_iter,"loglik"] <-
    new_lik
  if (control_glmpca_pois$calc_max_diff) {
    fastglmpca_internal$progress[
      fastglmpca_internal$main_loop_iter,"max_diff_l"] <-
        max(abs(fit$LL - fit0$LL))
    fastglmpca_internal$progress[
      fastglmpca_internal$main_loop_iter,"max_diff_f"] <-
        max(abs(fit$FF - fit0$FF))
  }
  if (control_glmpca_pois$calc_deriv) {
    if (inherits(Y,"sparseMatrix")) {
      fastglmpca_internal$progress[
        fastglmpca_internal$main_loop_iter,"max_deriv_f"] <-
          max(abs((deriv_prod(fit$LL,fit$FF) - fit$LL %*% Y) * FF_mask))
      fastglmpca_internal$progress[
        fastglmpca_internal$main_loop_iter,"max_deriv_l"] <-
          max(abs((deriv_prod(fit$FF,fit$LL) -
                   Matrix::tcrossprod(fit$FF,Y)) * LL_mask))
    } else {
      fastglmpca_internal$progress[
        fastglmpca_internal$main_loop_iter,"max_deriv_f"] <-
          max(abs(crossprod(exp(crossprod(fit$LL,fit$FF)) - Y,
                            t(fit$LL)) * FF_mask))
      fastglmpca_internal$progress[
        fastglmpca_internal$main_loop_iter,"max_deriv_l"] <-
          max(abs(crossprod(exp(crossprod(fit$FF,fit$LL)) - t(Y),
                            t(fit$FF)) * LL_mask))
    }
  }
  stop_time <- proc.time()
  fastglmpca_internal$progress[
    fastglmpca_internal$main_loop_iter,"time"] <-
      (stop_time - start_time)["elapsed"]
  if (verbose)
    cat(sprintf("Iteration %d: log-likelihood = %+0.12e\n",
                fastglmpca_internal$main_loop_iter,new_lik))
  return(fit2par(fit,update_indices_l,update_indices_f))
}

# Extract the model fit from the value of "par" provided by fpiter or
# daarem.
par2fit <- function (par, LL, FF, update_indices_l, update_indices_f) {
  n <- ncol(LL)
  i <- seq(1,n * length(update_indices_l))
  fit <- list(LL = LL,FF = FF)
  fit$LL[update_indices_l,] <- par[i]
  fit$FF[update_indices_f,] <- par[-i]
  return(fit)
}

# Convert the model fit to a "par" value accepted by fpiter or daarem.
fit2par <- function (fit, update_indices_l, update_indices_f)
  c(fit$LL[update_indices_l,],fit$FF[update_indices_f,])

# This implements a single update of LL and FF.
#
#' @importMethodsFrom Matrix tcrossprod
update_glmpca_pois <- function (LL, FF, Y, Y_T, update_indices_l,
                                update_indices_f, control) {
  n <- nrow(Y)
  m <- ncol(Y)
  K <- nrow(FF)
  k <- sort(intersect(update_indices_l,update_indices_f))
  if (length(update_indices_l) > 0) {

    # If requested, orthogonalize rows of FF that are not fixed.
    if (length(k) > 1 & control$orthonormalize) {
      out    <- svd(t(FF[k,,drop = FALSE]))
      FF[k,] <- t(out$u)
      LL[k,] <- diag(out$d) %*% t(out$v) %*% LL[k,,drop = FALSE]
    }

    # Update the LL matrix.
    LLnew  <- matrix(LL,K,n)
    i      <- update_indices_l - 1
    update_factors_faster_parallel(
        L_T = t(FF),
        FF = LLnew,
        M = as.matrix(tcrossprod(FF[update_indices_l,,drop = FALSE],Y)),
        update_indices = i,
        num_iter = control$num_ccd_iter,
        line_search = control$line_search,
        alpha = control$ls_alpha,
        beta = control$ls_beta)
    LL <- LLnew
  }

  if (length(update_indices_f) > 0) {

    # If requested, orthogonalize rows of LL that are not fixed.
    if (length(k) > 1 & control$orthonormalize) {
      out    <- svd(t(LL[k,,drop = FALSE]))
      LL[k,] <- t(out$u)
      FF[k,] <- diag(out$d) %*% t(out$v) %*% FF[k,,drop = FALSE]
    }

    # Update the FF matrix.
    FFnew <- matrix(FF,K,m)
    i     <- update_indices_f - 1
    update_factors_faster_parallel(
      L_T = t(LL),
      FF = FFnew,
      M = as.matrix(tcrossprod(LL[update_indices_f,,drop = FALSE],Y_T)),
      update_indices = i,
      num_iter = control$num_ccd_iter,
      line_search = control$line_search,
      alpha = control$ls_alpha,
      beta = control$ls_beta)
    FF <- FFnew
  }

  return(list(LL = LL,FF = FF))
}
