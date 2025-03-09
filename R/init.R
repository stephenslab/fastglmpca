#' @rdname fit_glmpca_pois
#'
#' @param U An optional argument giving the initial estimate of the
#'   loadings matrix. It should be an n x K matrix, where n is the
#'   number of rows in the counts matrix \code{Y}, and K > 0 is the rank
#'   of the matrix factorization. When \code{U} and \code{V} are not
#'   provided, input argument \code{K} should be specified instead.
#'
#' @param V An optional argument giving is the initial estimate of the
#'   factors matrix. It should be a m x K matrix, where m is the number
#'   of columns in the counts matrix \code{Y}, and K > 0 is the rank of
#'   the matrix factorization. When \code{U} and \code{V} are not
#'   provided, input argument \code{K} should be specified instead.
#'
#' @param X Optional argument giving row covariates of the count
#'   matrix \code{Y}. It should be an n x nx matrix, where nx is
#'   the number of row covariates.
#'
#' @param B Optional argument giving the initial estimates for the
#'   coefficients of the row covariates. It should be an m x nx matrix,
#'   where nx is the number of row covariates.
#'   This argument is ignored if X is not provided.
#'
#' @param Z Optional argument giving column covariates of the count
#'   matrix \code{Y}. It should be an m x nz matrix, where nz is the
#'   number of column covariates.
#'
#' @param W Optional argument giving the initial estimates for the
#'   coefficients of the column covariates.  It should be an n x nz matrix,
#'   where nz is the number of column covariates.
#'   This argument is ignored if Z is not provided.
#'
#' @param fixed_b_cols Optional numeric vector specifying which
#'   columns of \code{B} (if any) should be fixed during
#'   optimization. This argument is ignored if X is not provided.
#'
#' @param fixed_w_cols Optional numeric vector specifying which
#'   columns of \code{W} (if any) should be fixed during
#'   optimization. This argument is ignored if Z is not provided.
#'
#' @param col_size_factor If \code{col_size_factor = TRUE}, add a
#'   fixed factor accounting for average differences in Poisson rates
#'   across columns of \code{Y}. Setting \code{col_size_factor = TRUE}
#'   and \code{row_intercept = TRUE} is intended to replicate the
#'   default behavior of \code{glmpca}.
#'
#' @param row_intercept If \code{row_intercept = TRUE}, add a factor
#'   accounting for average differences in Poisson rates across rows of
#'   \code{Y}. Setting \code{col_size_factor = TRUE} and
#'   \code{row_intercept = TRUE} is intended to replicate the default
#'   behavior of \code{glmpca}.
#'
#' @seealso \code{\link{fit_glmpca_pois}}
#'
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom stats rnorm
#'
#' @export
#'
init_glmpca_pois <- function(
    Y,
    K,
    U,
    V,
    X = numeric(0),
    Z = numeric(0),
    B = numeric(0),
    W = numeric(0),
    fixed_b_cols = numeric(0),
    fixed_w_cols = numeric(0),
    col_size_factor = TRUE,
    row_intercept = TRUE) {

  # Check and prepare input argument Y.
  verify.count.matrix(Y)
  n <- nrow(Y)
  m <- ncol(Y)
  if (is.integer(Y))
    storage.mode(Y) <- "double"

  # Only one of K or (U, V) should be provided.
  if (!((missing(K) & !missing(U) & !missing(V)) |
        (!missing(K) & missing(U) & missing(V))))
    stop("Provide a rank (K) or an initialization of U and V, but not both")

  # Check and prepare input arguments K, U and V.
  if (missing(K)) {

    # Check the provided U and V.
    if (!(is.matrix(U) & is.numeric(U)))
      stop("Input argument \"U\" should be a numeric matrix (is.matrix(U) ",
           "should return TRUE)")
    if (is.integer(U))
      storage.mode(U) <- "double"
    if (!(is.matrix(V) & is.numeric(V)))
      stop("Input argument \"V\" should be a numeric matrix (is.matrix(V) ",
           "should return TRUE)")
    if (is.integer(V))
      storage.mode(V) <- "double"
    if (nrow(U) != nrow(Y))
      stop("Input argument \"U\" should have same number of rows as \"Y\"")
    if (nrow(V) != ncol(Y))
      stop("Input argument \"V\" should have same number of columns as \"Y\"")
    if (ncol(U) != ncol(V))
      stop("Inputs \"U\" and \"V\" should have same number of columns")
    K <- ncol(U)
  } else {
    if (K >= min(m,n))
      stop("Input \"K\" should be less than the number of rows and columns ",
           "of Y.")

    # Initialize U and V.
    U <- matrix(rnorm(n*K,sd = 1e-4),n,K)
    V <- matrix(rnorm(m*K,sd = 1e-4),m,K)
  }

  # Check and prepare input arguments X and B.
  if (!missing(X) && length(X) > 0) {
    verify.matrix(X)
    nx <- ncol(X)
    if (nrow(X) != nrow(Y))
      stop("Inputs \"X\" and \"Y\" should have same number of rows")
    if (missing(B))
      B <- matrix(rnorm(m*nx,sd = 0.1),m,nx)
    if (nrow(B) != ncol(Y))
      stop("Input \"B\" should have same number of rows as \"Y\" has columns")
    if (ncol(B) != ncol(X))
        stop("Inputs \"B\" and \"X\" should have same number of columns")
    if (is.integer(X))
      storage.mode(X) <- "double"
    if (is.integer(B))
      storage.mode(B) <- "double"
    if (is.null(colnames(X)))
      colnames(X) <- paste0("x_",1:nx)
  } else {
    B <- numeric(0)
    nx <- 0
    fixed_b_cols <- numeric(0)
  }

  # Check and prepare input arguments Z and W.
  if (!missing(Z) && length(Z) > 0) {
    verify.matrix(Z)
    nz <- ncol(Z)
    if (nrow(Z) != ncol(Y))
      stop("Input \"Z\" should have as many rows as columns of \"Y\"")
    if (missing(W))
      W <- matrix(rnorm(n*nz,sd = 0.1),n,nz)
    if (nrow(W) != nrow(Y))
      stop("Input \"W\" should have same number of rows as \"Y\"")
    if (ncol(W) != ncol(Z))
      stop("Inputs \"W\" and \"Z\" should have same number of columns")
    if (is.integer(Z))
      storage.mode(Z) <- "double"
    if (is.integer(W))
      storage.mode(W) <- "double"
    if (is.null(colnames(Z)))
      colnames(Z) <- paste0("z_",1:nz)
  } else {
    W <- numeric(0)
    nz <- 0
    fixed_w_cols <- numeric(0)
  }

  # Add the fixed intercept and/or size factor if requested.
  if (col_size_factor) {
    x <- matrix(rep(1,nrow(Y)))
    b <- matrix(log(colMeans(Y)))
    colnames(x) <- "col_size_factor"
    colnames(b) <- "col_size_factor"
    X <- cbind(X,x)
    B <- cbind(B,b)
    nx <- nx + 1
    fixed_b_cols <- c(fixed_b_cols,nx)
  }
  if (row_intercept) {
    z <- matrix(rep(1,ncol(Y)))
    w <- matrix(log(rowSums(Y)/sum(colMeans(Y))))
    colnames(z) <- "row_intercept"
    colnames(w) <- "row_intercept"
    Z <- cbind(Z,z)
    W <- cbind(W,w)
    nz <- nz + 1
  }

  # Compute the log-likelihood.
  LL <- t(cbind(U,X,W))
  FF <- t(cbind(V,B,Z))
  if (inherits(Y,"sparseMatrix"))
    loglik <-
      lik_glmpca_pois_log_sp(Y,LL,FF,sum(lfactorial(Y@x)))
  else
    loglik <- lik_glmpca_pois_log(Y,LL,FF,sum(lfactorial(Y)))

  # Prepare the final output.
  fit <- list(U = U,V = V,X = X,B = B,Z = Z,W = W,
              fixed_b_cols = fixed_b_cols,
              fixed_w_cols = fixed_w_cols,
              loglik = loglik,
              progress = data.frame(iter        = 0,
                                    loglik      = loglik,
                                    time        = 0,
                                    max_deriv_f = as.numeric(NA),
                                    max_deriv_l = as.numeric(NA),
                                    max_diff_f  = as.numeric(NA),
                                    max_diff_l  = as.numeric(NA)))
  fit <- orthonormalize_fit(fit)
  rownames(fit$U) <- rownames(Y)
  rownames(fit$V) <- colnames(Y)
  colnames(fit$U) <- paste("k",1:K,sep = "_")
  colnames(fit$V) <- paste("k",1:K,sep = "_")
  names(fit$d)    <- paste("k",1:K,sep = "_")
  if (length(fit$X) > 0) {
    rownames(fit$X) <- rownames(Y)
    rownames(fit$B) <- colnames(Y)
    colnames(fit$B) <- colnames(fit$X)
  }
  if (length(fit$Z) > 0) {
    rownames(fit$Z) <- colnames(Y)
    rownames(fit$W) <- rownames(Y)
    colnames(fit$W) <- colnames(fit$Z)
  }
  class(fit) <- c("glmpca_pois_fit","list")
  return(fit)
}
