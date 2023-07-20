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
#' @param fit_col_size_factor Boolean indicating if a size factor should be
#'   used to normalize the likelihood across columns. This is done by fixing
#'   an element of each loading to 1 and fixing and element of each factor
#'   as the log the mean value of it's corresponding column. This may be
#'   useful when \code{Y} is a matrix from a scRNA experiment where rows 
#'   represent genes and columns represent cells.
#' 
#' @param fit_row_intercept Boolean indicating if intercept term should be fit
#'   for each row of \code{Y}. This is done by fixing an element of each
#'   factor to 1. This may be useful when \code{Y} is a matrix from a
#'   scRNA experiment where rows represent genes and columns represent cells,
#'   and one wants to regress out mean differences between genes.
#'
#' @param fixed_loadings Vector of integers indicating, which, if any, loadings
#'   (i.e. columns of \code{U}) should be fixed at their initial values. 
#'   This argument will be ignored if \code{U} is not provided.
#'
#' @param fixed_factors Vector of integers indicating which, if any, factors
#'   (i.e. columns of \code{V}) should be fixed at their initial values. 
#'   This argument will be ignored if \code{V} is not provided.
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
    fit_col_size_factor = FALSE,
    fit_row_intercept = FALSE,
    fixed_loadings = NULL,
    fixed_factors = NULL
) {
  
  if (!missing(U) && !missing(V)) {
    
    if (ncol(U) != ncol(V)) {
      
      stop("Inputs \"U\" and \"V\" must have same number of columns")
      
    }
    
  }
  
  if (!missing(Y) && !missing(U)) {
    
    if (nrow(U) != nrow(Y)) {
      
      stop("Input \"U\" must have same number of rows as \"Y\" ")
      
    }
    
  }
  
  if (!missing(Y) && !missing(V)) {
    if (nrow(V) != ncol(Y)) {
      
      stop("Input \"V\" must have same number of rows as there are columns of \"Y\"")
      
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
  
  if (missing(U)) {
    
    if (missing(K)) {
      
      stop("if \"U\" is missing, must provide \"K\" ")
      
    }
      
    fit$U <- matrix(
      data = rnorm(K_total * n, 0, sd=1e-5/(K + as.numeric(fit_row_intercept))),
      ncol = K_total, 
      nrow = n
    )
    
    if (fit_col_size_factor) {
      
      if (missing(Y)) {
        
        stop("if \"fit_col_size_factor\" is true, must provide \"Y\" ")
        
      }
      
      fit$U[, 1] <- 1
      fit$fixed_loadings <- c(1)
      
      if (fit_row_intercept) {
        
        fit$U[, 2] <- log(rowSums(Y) / sum(colMeans(Y)))
        
      }
      
    } else {
      
      fit$fixed_loadings <- NULL
      
      if (fit_row_intercept) {
        
        if (missing(Y)) {
          
          stop("if \"fit_row_intercept\" is true, must provide \"Y\" ")
          
        }
        
        fit$U[ ,1] <- log(rowSums(Y) / sum(colMeans(Y)))
        
      }
      
    }
    
  } else {
    
    fit$U <- U
    fit$fixed_loadings <- fixed_loadings
    
  }
  
  if (missing(V)) {
    
    if (missing(K)) {
      
      stop("if \"V\" is missing, must provide \"K\" ")
      
    }
      
    fit$V <- matrix(
      data = rnorm(K_total * p, 0, sd = 1e-5 / K), 
      ncol = K_total, 
      nrow = p
    )
    
    if (fit_col_size_factor && fit_row_intercept) {
      
      if (missing(Y)) {
        
        stop("if \"fit_col_size_factor\" is TRUE, must provide \"Y\"")
        
      }
      
      fit$V[, 1] <- log(colMeans(Y))
      
      # Intercept
      fit$V[, 2] <- 1
      
      fit$fixed_factors <- c(1, 2)
      
    } else if (fit_col_size_factor) {
      
      if (missing(Y)) {
        
        stop("if \"fit_col_size_factor\" is TRUE, must provide \"Y\"")
        
      }
      
      fit$V[ ,1] <- log(colMeans(Y))
      fit$fixed_factors <- c(1)
      
    } else if (fit_row_intercept) {
      
      fit$V[ ,1] <- 1
      fit$fixed_factors <- c(1)
      
    } else {
      
      fit$fixed_factors <- NULL
      
    }
    
  } else {
    
    fit$V <- V
    fit$fixed_factors <- fixed_factors
    
  }
  
  if (missing(U)) {
    
    colnames_U <- paste0("loading_", c(1:K))
    if (fit_row_intercept) {
      
      colnames_U <- c("intercept", colnames_U)
      
    }
    
    if (fit_col_size_factor) {
      
      colnames_U <- c("size_factor", colnames_U)
      
    }
    
    colnames(fit$U) <- colnames_U
    
  }
  
  if (missing(V)) {
    
    colnames_V <- paste0("factor_", c(1:K))
    if (fit_row_intercept) {
      
      colnames_V <- c("intercept", colnames_V)
      
    }
    
    if (fit_col_size_factor) {
      
      colnames_V <- c("size_factor", colnames_V)
      
    }
    
    colnames(fit$V) <- colnames_V
    
  }
  
  fit <- orthonormalize_fit(fit)
  
  class(fit) <- c("glmpca_pois_fit", "list")
  
  return(fit)
  
}

#' @title Initialize Binomial GLM-PCA Fit
#' 
#' @description Initialize a binomial GLM-PCA model fit.
#'
#' @param Y The n x p matrix of counts; all entries of Y should be
#'   non-negative.
#'   
#' @param N The n x p matrix of trials; all entries of N should be
#'   at least 1.
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
#' @seealso \code{\link{fit_glmpca_binom}}
#'
#' @importFrom stats rnorm
#' 
#' @export
#' 
init_glmpca_binom <- function(
    Y,
    K,
    LL,
    FF,
    fixed_loadings = NULL,
    fixed_factors = NULL
) {
  
  fit <- list()
  n <- nrow(Y)
  p <- ncol(Y)
  
  if(missing(LL)) {
    
    fit$LL <- matrix(
      data = rnorm(K * n, sd = .1),
      nrow = K,
      ncol = n
    )
    
  } else {
    
    fit$LL <- LL
    
  }
  
  if(missing(FF)) {
    
    fit$FF <- matrix(
      data = rnorm(K * p, sd = .1),
      nrow = K,
      ncol = p
    )
    
  } else {
    
    fit$FF <- FF
    
  }
  
  fit$fixed_loadings <- fixed_loadings
  fit$fixed_factors <- fixed_factors
  
  class(fit) <- c("glmpca_fit", "list")
  
  return(fit)
  
}
