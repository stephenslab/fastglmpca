#' @title Generate Data from a GLMCPA Model
#' 
#' @description Generate data from a GLMPCA model with a specified
#'   rank.
#'   
#' @details This function assumes that each column of the data is generated
#'   from a multinomial distribution. Let \deqn{Y_j} denote column j of the 
#'   generated data matrix. First, we set \deqn{sum(Y_j)} equal to a value 
#'   generated from a \deqn{Uniform(50, 5000)} distribution. Then,
#'   we generate \deqn{L} and \deqn{F} from mixture distributions,
#'   and calculate \deqn{H = exp(L'F)}. Then, we generate the individual
#'   elements of \deqn{Y_j} from a multinomial model where the probability
#'   for each individual element is just \deqn{H_j} normalized.
#'
#' @param n Number of rows (genes).
#' @param p Number of columns (cells).
#' @param K Rank of the underlying mean structure
#'
#' @return list with the following components
#' \itemize{
#'   \item LL - loadings of underlying mean structure. A K x n matrix
#'   \item FF - factors of underlying mean structure. A K x p matrix
#'   \item Y - n x p matrix of generated data.
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' sim_data <- generate_glmpca_data(1000, 500, 1)
#' 
generate_glmpca_data <- function(n, p, K) {
  
  if (!is.scalar(K) || K < 1) {
    
    stop("\"K\" must be an integer greater than or equal to 1")
    
  }
  
  if (!is.scalar(n) || n < 1) {
    
    stop("\"n\" must be an integer greater than or equal to 1")
    
  }
  
  if (!is.scalar(p) || p < 1) {
    
    stop("\"p\" must be an integer greater than or equal to 1")
    
  }
  
  LL <- matrix(nrow = K, ncol = n)
  FF <- matrix(nrow = K, ncol = p)
  
  l_dist <- distr::UnivarMixingDistribution(
    distr::Unif(0, .01),
    distr::Dirac(1.35),
    distr::Dirac(2.35),
    mixCoeff = rep((1/3), 3)
  )
  
  f_dist <- distr::UnivarMixingDistribution(
    distr::Unif(0, .01),
    distr::Dirac(1.35),
    distr::Dirac(2.35),
    mixCoeff = rep((1/3), 3)
  )
  
  l_sampler <- distr::r(l_dist)
  f_sampler <- distr::r(f_dist)
  
  cell_totals <- runif(n = p, min = 50, max = 5000)
  
  for (k in 1:K) {
    
    LL[k, ] <- l_sampler(n)
    FF[k, ] <- f_sampler(p)
    
  }
  
  c_dist <- distr::UnivarMixingDistribution(
    distr::Dirac(3),
    distr::Dirac(0),
    distr::Unif(0, 3),
    mixCoeff = c(.6, .2, .2)
  )
  
  Lambda <- exp(crossprod(LL, FF)) 
  
  Pi <- sweep(Lambda, 2, colSums(Lambda), `/`)
  
  Y <- matrix(nrow = n, ncol = p)
  
  for (j in 1:p) {
    
    Y[, j] <- rmultinom(n = 1, size = cell_totals[j], prob = Pi[, j])
    
  }
  
  return(
    list(
      Y = Y, LL = LL, FF = FF
    )
  )
  
}

