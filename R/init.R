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
#'   matrix \code{Y}. It should be an n x nx matrix, where n is the
#'   number of rows in \code{Y} and nx is the number of row
#'   covariates.
#' 
#' @param B An optional argument giving the initial estimates
#'   for the coefficients of the row specific covariates of the 
#'   count matrix \code{Y}. It should be p x n_x matrix, where p
#'   is the number of columns of the count matrix \code{Y}, and
#'   n_x is the number of row covariates. This argument is ignored
#'   if X is not provided.
#' 
#' @param Z An optional argument giving column specific covariates of the 
#'   count matrix \code{Y}. It should be p x n_z matrix, where p
#'   is the number of columns of the count matrix \code{Y}, and
#'   n_z is the number of column covariates.
#' 
#' @param W An optional argument giving the initial estimates
#'   for the coefficients of the column specific covariates of the 
#'   count matrix \code{Y}. It should be n x n_z matrix, where n
#'   is the number of rows of the count matrix \code{Y}, and
#'   n_z is the number of column covariates.
#'   
#' @param fixed_b_cols An optional argument giving which, if any, columns
#'   of \code{B} should be fixed during optimization.
#'   
#' @param fixed_w_cols An optional argument giving which, if any, columns
#'   of \code{W} should be fixed during optimization.
#'   
#' @param col_size_factor Boolean indicating if a size factor should be
#'   used to normalize the likelihood across columns. This may be
#'   useful when \code{Y} is a matrix from a scRNA experiment where rows 
#'   represent genes and columns represent cells.
#'   
#' @param row_size_factor Boolean indicating if a size factor should be
#'   used to normalize the likelihood across rows. This may be
#'   useful when \code{Y} is a matrix from a scRNA experiment where rows 
#'   represent cells and columns represent genes.
#' 
#' @param row_intercept Boolean indicating if intercept term should be
#' fit for each row of \code{Y}. This may be useful when \code{Y} is a
#' matrix from a scRNA experiment where rows represent genes and
#' columns represent cells, and one wants to regress out mean
#' differences between genes.
#'   
#' @param col_intercept Boolean indicating if intercept term should be
#' fit for each column of \code{Y}. This may be useful when \code{Y}
#' is a matrix from a scRNA experiment where rows represent cells and
#' columns represent genes, and one wants to regress out mean
#' differences between genes.
#'
#' @return An object capturing the initial state of the model fit. See
#'   \code{\link{fit_glmpca_pois}} for details.
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
    row_intercept = TRUE,
    row_size_factor = FALSE,
    col_intercept = FALSE
) {

  # Check and prepare input argument Y.
  verify.count.matrix(Y)
  n <- nrow(Y)
  m <- ncol(Y)
  
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

    # Initialize U and V.
    U <- matrix(rnorm(n*K,sd = 0.1),n,K)
    V <- matrix(rnorm(m*K,sd = 0.1),m,K)
  }

  # Check and prepare input arguments X and B.
  if (!missing(X) && length(X) > 0) {
    if (nrow(X) != nrow(Y))
      stop("Inputs \"X\" and \"Y\" must have same number of rows")
    nx <- ncol(X)
    if (missing(B))
      B <- matrix(rnorm(m*n_x,sd = 0.1),m,nx)      
    if (nrow(B) != ncol(Y))
      stop("Input \"B\" must have same number of rows as \"Y\" has columns")
    if (ncol(B) != ncol(X))
        stop("Inputs \"B\" and \"X\" must have same number of columns")
  } else
    B <- numeric(0)

  # Check and prepare input arguments Z and W.
  if (!missing(X) && length(X) > 0) {

  }
  
  # Prepare the final output.
  fit <- list(U = U,V = V,X = X,Z = Z,B = B,W = W)
  fit <- orthonormalize_fit(fit)
  class(fit) <- c("glmpca_pois_fit","list")
  rownames(fit$U) <- rownames(Y)
  rownames(fit$V) <- colnames(Y)
  colnames(fit$U) <- paste("k",1:K,sep = "_")
  colnames(fit$V) <- paste("k",1:K,sep = "_")
  rownames(fit$D) <- paste("k",1:K,sep = "_")
  colnames(fit$D) <- paste("k",1:K,sep = "_")
  if (length(fit$X) > 0) {
    rownames(fit$X) <- rownames(Y)
    if (is.null(colnames(fit$X)))
      colnames(fit$X) <- paste0("x_",1:nx)
    rownames(fit$B) <- colnames(Y)
    colnames(fit$B) <- colnames(fit$X)
  }
  return(fit)
  
  if (missing(Z)) {
    
    fit$Z <- numeric(0)
    fit$W <- numeric(0)
    
    n_z <- 0
    
  }
  else {
    
    if (nrow(Z) != ncol(Y)) {
      
      stop("Input \"Z\" must have same number of rows as \"Y\" has columns")
      
    }
    
    nz <- ncol(Z)
    rownames(fit$Z) <- colnames(Y)
    if(is.null(colnames(fit$Z))) {
      
      colnames(fit$Z) <- paste0("z_", 1:n_z)
      
    }
    
    if(missing(W)) {
      
      W <- matrix(
        data = rnorm(n = n * n_z, sd = 1e-5),
        nrow = n,
        ncol = n_z
      )
      
    } else {
      
      if (nrow(W) != nrow(Y)) {
        
        stop("Inputs \"W\" and \"Y\" must have same number of rows")
        
      }
      
      if (ncol(W) != ncol(Z)) {
        
        stop("Inputs \"W\" and \"Z\" must have same number of columns")
        
      }
      
    }
    
    fit$W <- W
    rownames(W) <- rownames(Y)
    if(is.null(colnames(W))) {
      
      colnames(W) <- paste0("w_", 1:n_z)
      
    }
    
  }
  
  # Now, want to add various row or column intercepts / size factors
  if (fit_row_intercept) {
    
    Z_int <- rep(1, p)
    names(Z_int) <- colnames(Y)
    
    W_int <- log(Matrix::rowSums(Y) / sum(Matrix::colMeans(Y)))
    names(W_int) <- rownames(Y)
    
    fit$Z <- cbind(fit$Z, Z_int)
    fit$W <- cbind(fit$W, W_int)
    n_z <- n_z + 1
    
    colnames(fit$Z)[ncol(fit$Z)] <- "row_intercept"
    colnames(fit$W)[ncol(fit$W)] <- "row_intercept"
    
  }
  
  if (fit_col_intercept) {
    
    X_int <- rep(1, n)
    names(X_int) <- rownames(Y)
    
    B_int <- log(Matrix::colSums(Y) / sum(Matrix::rowMeans(Y)))
    names(B_int) <- colnames(Y)
    
    fit$X <- cbind(fit$X, X_int)
    fit$B <- cbind(fit$B, B_int)
    n_x <- n_x + 1
    
    colnames(fit$X)[ncol(fit$X)] <- "col_intercept"
    colnames(fit$B)[ncol(fit$B)] <- "col_intercept"
    
  }
  
  if (fit_row_size_factor) {
    
    Z_size <- rep(1, p)
    names(Z_size) <- colnames(Y)
    
    W_size <- log(Matrix::rowMeans(Y))
    names(W_size) <- rownames(Y)
    
    fit$Z <- cbind(fit$Z, Z_size)
    fit$W <- cbind(fit$W, W_size)
    n_z <- n_z + 1
    fit$fixed_w_cols <- c(fit$fixed_w_cols, n_z)
    
    colnames(fit$Z)[ncol(fit$Z)] <- "row_size_factor"
    colnames(fit$W)[ncol(fit$W)] <- "row_size_factor"
    
  } 
  
  if (fit_col_size_factor) {
    
    X_size <- rep(1, n)
    names(X_size) <- rownames(Y)
    
    B_size <- log(Matrix::colMeans(Y))
    names(B_size) <- colnames(Y)
    
    fit$X <- cbind(fit$X, X_size)
    fit$B <- cbind(fit$B, B_size)
    n_x <- n_x + 1
    fit$fixed_b_cols <- c(fit$fixed_b_cols, n_x)
    
    colnames(fit$X)[ncol(fit$X)] <- "col_size_factor"
    colnames(fit$B)[ncol(fit$B)] <- "col_size_factor"
    
  } 
  
  # here, want to calculate the loglik to add to the fit
  if (!inherits(Y, "sparseMatrix")) {
    
    H <- tcrossprod(fit$U %*% fit$D, fit$V)
    
    if(!identical(fit$X, numeric(0))) {
      
      H <- H + tcrossprod(fit$X, fit$B)
      
    }
    
    if(!identical(fit$Z, numeric(0))) {
      
      H <- H + tcrossprod(fit$W, fit$Z)
      
    }
    
    loglik <- sum(Y * H - exp(H)) - sum(lfactorial(Y))
    
  } else {
    
    loglik <- lik_glmpca_pois_log_sp(
      Y = Y, 
      LL = t(
        cbind(
          fit$U %*% fit$D,
          fit$X,
          fit$W
        )
      ),
      FF = t(
        cbind(
          fit$V,
          fit$B,
          fit$Z
        )
      ),
      const = sum(mapSparse(Y, lfactorial))
    )
    
  }
  
  fit$loglik <- loglik
  fit$progress <- data.frame(
    iter = 0,
    loglik = loglik,
    time = 0
  )
  
  return(fit)
  
}
