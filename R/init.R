#' @title Initialize GLM-PCA Fit
#' 
#' @description Initialize a GLM-PCA model fit.
#'
#' @param Y The n x p matrix of counts; all entries of Y should be
#'   non-negative.
#'   
#' @param K An integer 1 or greater giving the matrix rank. This
#'   argument should only be specified if the initial fit (\code{LL, FF})
#'   is not provided.
#'   
#' @param LL An optional argument giving the initial estimate of the
#'   loadings matrix. It should be an K x n matrix, where n is the
#'   number of rows in the counts matrix \code{Y}, and K >= 1 is the rank
#'   of the matrix factorization. When \code{LL} and \code{FF} are not
#'   provided, input argument \code{K} should be specified instead, and
#'   \code{LL} and \code{FF} are initialized as a constant matrix.
#'   
#' @param FF An optional argument giving is the initial estimate of the
#'   factors matrix. It should be a K x p matrix, where p is the number
#'   of columns in the counts matrix \code{Y}, and K >= 1 is the rank of
#'   the matrix factorization. When \code{LL} and \code{FF} are not
#'   provided, input argument \code{K} should be specified instead, and
#'   \code{LL} and \code{FF} are initialized randomly
#'   
#' @param fit_col_size_factor Boolean indicating if a size factor should be
#'   used to normalize the likelihood across columns. This is done by fixing
#'   an element of each loading to 1 and fixing and element of each factor
#'   as the log the mean value of it's corresponding column. This may be
#'   useful when \code{Y} is a matrix from a scRNA experiment where rows 
#'   represent genes and columns represent cells, and one wants to control
#'   for different sequencing depths between cells.
#' 
#' @param fit_row_intercept Boolean indicating if intercept term should be fit
#'   for each row of \code{Y}. This is done by fixing an element of each
#'   factor to 1. This may be useful when \code{Y} is a matrix from a
#'   scRNA experiment where rows represent genes and columns represent cells.
#'
#' @param fixed_loadings Vector of integers indicating, which, if any, loadings
#'   should be fixed at their initial values. This argument will be ignored if
#'   \code{LL} is not provided.
#'
#' @param fixed_factors Vector of integers indicating which, if any, factors
#'   should be fixed at their initial values. This argument will be ignored
#'   if \code{FF} is not provided.
#'
#' @return An object capturing the initial state of the model fit. See
#'   \code{\link{fit_glmpca}} for details.
#'
#' @seealso \code{\link{fit_glmpca}}
#' 
#' @export
#' 
init_glmpca <- function(
    Y,
    K,
    LL,
    FF,
    fit_col_size_factor = FALSE,
    fit_row_intercept = FALSE,
    fixed_loadings = NULL,
    fixed_factors = NULL
) {
  
  if (!missing(LL) && !missing(FF)) {
    
    if (nrow(LL) != nrow(FF)) {
      
      stop("Inputs \"LL\" and \"FF\" must have same number of rows")
      
    }
    
  }
  
  if (!missing(Y) && !missing(LL)) {
    
    if (ncol(LL) != nrow(Y)) {
      
      stop("Input \"LL\" must have same number of columns as \"Y\" has rows")
      
    }
    
  }
  
  if (!missing(Y) && !missing(FF)) {
    if (ncol(FF) != ncol(Y)) {
      
      stop("Input \"FF\" must have same number of columns as \"Y\"")
      
    }
    
  }
  
  if (!missing(Y)) {
    
    n <- nrow(Y)
    p <- ncol(Y)
    
  }
  
  fit <- list()
  
  if(!missing(K)) {
    
    if(!is.scalar(K) || K < 1) 
      stop("Input argument \"K\" must be a positive integer")
    
    K_total <- K + fit_col_size_factor + fit_row_intercept
    
  }
  
  
  if (missing(LL)) {
    
    if (missing(K) || missing(Y)) {
      
      stop("if \"LL\" is missing, must provide \"K\" and \"Y\"")
      
    }
    
    fit$LL <- matrix(
      data = runif(K_total * n, .Machine$double.eps, 1),
      nrow = K_total, 
      ncol = n
    )
    
    if (fit_col_size_factor) {
      
      fit$LL[1, ] <- 1
      fit$fixed_loadings <- c(1)
      
      rownames(fit$LL) <- c(
        "size_factor", 
        paste0("loading_", c(1:(K + fit_row_intercept)))
      )
      
    } else {
      
      fit$fixed_loadings <- NULL
      
    }
    
  } else {
    
    fit$LL <- LL
    fit$fixed_loadings <- fixed_loadings
    
  }
  
  if (missing(FF)) {
    
    if (missing(K) || missing(Y)) {
      
      stop("if \"FF\" is missing, must provide \"K\" and \"Y\"")
      
    }
    
    fit$FF <- matrix(
      data = runif(K_total * p, .Machine$double.eps, 1), nrow = K_total, ncol = p
    )
    
    if (fit_col_size_factor && fit_row_intercept) {
      
      if (missing(Y)) {
        
        stop("if \"fit_col_size_factor\" is TRUE, must provide \"Y\"")
        
      }
      
      fit$FF[1, ] <- log(Matrix::colMeans(Y))
      
      # Intercept
      fit$FF[2, ] <- 1
      
      fit$fixed_factors <- c(1, 2)
      rownames(fit$FF) <- c("size_factor", "intercept", paste0("factor_", c(1:K)))
      
    } else if (fit_col_size_factor) {
      
      if (missing(Y)) {
        
        stop("if \"fit_col_size_factor\" is TRUE, must provide \"Y\"")
        
      }
      
      fit$FF[1, ] <- log(Matrix::colMeans(Y))
      fit$fixed_factors <- c(1)
      rownames(fit$FF) <- c("size_factor", paste0("factor_", c(1:K)))
      
    } else if (fit_row_intercept) {
      
      fit$FF[1, ] <- 1
      fit$fixed_factors <- c(1)
      rownames(fit$FF) <- c("intercept", paste0("factor_", c(1:K)))
      
    } else {
      
      fit$fixed_factors <- NULL
      
    }
    
  } else {
    
    fit$FF <- init_FF
    fit$fixed_factors <- fixed_factors
    
  }
  
  class(fit) <- c("glmpca_fit", "list")
  
  return(fit)
  
}