#' @title Initialize GLM-PCA Poisson Fit
#' 
#' @description Initialize a GLM-PCA Poisson model fit.
#'
#' @param Y n x p matrix of counts; all entries should be
#'   non-negative.
#'   
#' @param K An integer 1 or greater giving the matrix rank. This
#'   argument should only be specified if the initial fit (\code{U, V})
#'   is not provided.
#'   
#' @param U An optional argument giving the initial estimate of the
#'   loadings matrix. It should be an n x K matrix, where n is the
#'   number of rows in the counts matrix \code{Y}, and K >= 1 is the rank
#'   of the matrix factorization. When \code{U} and \code{V} are not
#'   provided, input argument \code{K} should be specified instead.
#'   
#' @param V An optional argument giving is the initial estimate of the
#'   factors matrix. It should be a p x K matrix, where p is the number
#'   of columns in the counts matrix \code{Y}, and K >= 1 is the rank of
#'   the matrix factorization. When \code{U} and \code{V} are not
#'   provided, input argument \code{K} should be specified instead.
#'   
#' @param X An optional argument giving row specific covariates of the 
#'   count matrix \code{Y}. It should be n x n_x matrix, where n
#'   is the number of rows of the count matrix \code{Y}, and
#'   n_x is the number of row covariates.
#' 
#' @param B An optional argument giving the initial estimates
#'   for the coefficients of the row specific covariates of the 
#'   count matrix \code{Y}. It should be p x n_x matrix, where p
#'   is the number of columns of the count matrix \code{Y}, and
#'   n_x is the number of row covariates.
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
#' @param fit_col_size_factor Boolean indicating if a size factor should be
#'   used to normalize the likelihood across columns. This may be
#'   useful when \code{Y} is a matrix from a scRNA experiment where rows 
#'   represent genes and columns represent cells.
#'   
#' @param fit_row_size_factor Boolean indicating if a size factor should be
#'   used to normalize the likelihood across rows. This may be
#'   useful when \code{Y} is a matrix from a scRNA experiment where rows 
#'   represent cells and columns represent genes.
#' 
#' @param fit_row_intercept Boolean indicating if intercept term should be fit
#'   for each row of \code{Y}. This may be useful when \code{Y} is a matrix from a
#'   scRNA experiment where rows represent genes and columns represent cells,
#'   and one wants to regress out mean differences between genes.
#'   
#' @param fit_col_intercept Boolean indicating if intercept term should be fit
#'   for each column of \code{Y}. This may be useful when \code{Y} is a matrix from a
#'   scRNA experiment where rows represent cells and columns represent genes,
#'   and one wants to regress out mean differences between genes.
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
    X,
    Z,
    B,
    W,
    fixed_b_cols,
    fixed_w_cols,
    fit_col_size_factor = TRUE,
    fit_row_intercept = TRUE,
    fit_row_size_factor = FALSE,
    fit_col_intercept = FALSE
) {
    
  verify.count.matrix(Y)
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  fit <- list()
  
  if (missing(fixed_w_cols)) {
    
    fit$fixed_w_cols <- numeric(0)
    
  } else {
    
    fit$fixed_w_cols
    
  }
  
  if (missing(fixed_b_cols)) {
    
    fit$fixed_b_cols <- numeric(0)
    
  } else {
    
    fit$fixed_b_cols
    
  }
  
  if (missing(U)) {
    
    if (missing(K)) {
      
      stop("if \"U\" is missing, must provide \"K\" ")
      
    }
      
    fit$U <- matrix(
      data = rnorm(K, 0, sd=1e-5),
      ncol = K, 
      nrow = n
    )
    
  } else {
    
    if (nrow(U) != nrow(Y)) {
      
      stop("Input \"U\" must have same number of rows as \"Y\" ")
      
    }
    
    fit$U <- U
    
  }
  
  if (missing(V)) {
    
    if (missing(K)) {
      
      stop("if \"V\" is missing, must provide \"K\" ")
      
    }
      
    fit$V <- matrix(
      data = rnorm(K * p, 0, sd = 1e-5), 
      ncol = K, 
      nrow = p
    )
    
  } else {
    
    if (nrow(V) != ncol(Y)) {
      
      stop("Input \"V\" must have same number of rows as there are columns of \"Y\"")
      
    }
    
    fit$V <- V
    
  }
  
  if (!missing(U) && !missing(V)) {
    
    if (ncol(U) != ncol(V)) {
      
      stop("Inputs \"U\" and \"V\" must have same number of columns")
      
    }
    
  }
  
  if (missing(X)) {
    
    fit$X <- numeric(0)
    fit$B <- numeric(0)
    
    n_x <- 0
    
  }
  else {
    
    if (nrow(X) != nrow(Y)) {
      
      stop("Inputs \"X\" and \"Y\" must have same number of rows")
      
    }
    
    n_x <- ncol(X)
    
    fit$X <- X
    
    rownames(fit$X) <- rownames(Y)
    
    if(is.null(colnames(fit$X))) {
      
      colnames(fit$X) <- paste0("x_", 1:n_x)
      
    }
    
    if(missing(B)) {
      
      B <- matrix(
        data = rnorm(n = p * n_x, sd = 1e-5),
        nrow = p,
        ncol = n_x
      )
      
    } else {
      
      if (nrow(B) != ncol(Y)) {
        
        stop("Input \"B\" must have same number of rows as \"Y\" has columns")
        
      }
      
      if (ncol(B) != ncol(X)) {
        
        stop("Inputs \"B\" and \"X\" must have same number of columns")
        
      }
      
    }
    
    fit$B <- B
    
    rownames(fit$B) <- colnames(Y)
    if(is.null(colnames(fit$B))) {
      
      colnames(fit$B) <- paste0("b_", 1:n_x)
      
    }
    
  }
  
  if (missing(Z)) {
    
    fit$Z <- numeric(0)
    fit$W <- numeric(0)
    
    n_z <- 0
    
  }
  else {
    
    if (nrow(Z) != ncol(Y)) {
      
      stop("Input \"Z\" must have same number of rows as \"Y\" has columns")
      
    }
    
    n_z <- ncol(Z)
    
    fit$Z <- Z
    
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
  
  fit <- orthonormalize_fit_qr(fit)
  
  class(fit) <- c("glmpca_pois_fit", "list")
  
  rownames(fit$U) <- rownames(Y)
  colnames(fit$U) <- paste0("k_", 1:ncol(fit$U))
  rownames(fit$V) <- colnames(Y)
  colnames(fit$V) <- paste0("k_", 1:ncol(fit$V))
  rownames(fit$D) <- paste0("k_", 1:nrow(fit$D))
  colnames(fit$D) <- paste0("k_", 1:ncol(fit$D))
  
  return(fit)
  
}
