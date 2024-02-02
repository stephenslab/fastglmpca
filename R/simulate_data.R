#' @title Generate Data from a GLM-PCA Model
#' 
#' @description Generate data from a GLM-PCA model with a specified
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
#' 
#' @param p Number of columns (cells).
#' 
#' @param K Rank of the underlying mean structure.
#' 
#' @param link Character vector describing the link between the product 
#'   of the loading and factors and the mean of the data.
#'
#' @return list with the following components
#' \itemize{
#'   \item LL - loadings of underlying mean structure. A K x n matrix
#'   \item FF - factors of underlying mean structure. A K x p matrix
#'   \item Y - n x p matrix of generated data.
#' }
#'
#' @examples
#' set.seed(1)
#' sim_data <- generate_glmpca_data_pois(1000, 500, 1)
#' 
#' @importFrom stats runif
#' @importFrom stats rmultinom
#' @importFrom distr UnivarMixingDistribution
#' @importFrom distr Unif
#' @importFrom distr Dirac
#' @importFrom distr r
#' 
#' @export
#'
generate_glmpca_data_pois <- function(n, p, K, link = c("log", "log1p")) {
  if (!is.scalar(K) || K < 1)
    stop("\"K\" must be an integer greater than or equal to 1")
  if (!is.scalar(n) || n < 1)
    stop("\"n\" must be an integer greater than or equal to 1")
  if (!is.scalar(p) || p < 1)
    stop("\"p\" must be an integer greater than or equal to 1")
  
  link <- match.arg(link)
  
  LL <- matrix(nrow = K,ncol = n)
  FF <- matrix(nrow = K,ncol = p)
  
  l_dist <- UnivarMixingDistribution(
    Unif(0,0.01),
    Dirac(1.35),
    Dirac(2.35),
    mixCoeff = rep(1/3,3)
  )
  
  f_dist <- UnivarMixingDistribution(
    Unif(0,0.01),
    Dirac(1.35),
    Dirac(2.35),
    mixCoeff = rep(1/3,3)
  )
  
  l_sampler <- distr::r(l_dist)
  f_sampler <- distr::r(f_dist)
  cell_totals <- runif(n = p,min = 50,max = 5000)
  
  for (k in 1:K) {
    LL[k,] <- l_sampler(n)
    FF[k,] <- f_sampler(p)
  }
  
  if (link == "log")
    Lambda <- exp(crossprod(LL,FF)) 
  else if (link == "log1p")
    Lambda <- exp(crossprod(LL,FF)) - 1
  Pi <- sweep(Lambda,2,colSums(Lambda),`/`)
  Y <- matrix(nrow = n,ncol = p)
  
  for (j in 1:p)
    Y[,j] <- rmultinom(n = 1,size = cell_totals[j],prob = Pi[,j])

  rownames(Y) <- paste0("s_",1:n)
  colnames(Y) <- paste0("g_",1:p)
  return(list(Y = Y, LL = LL, FF = FF))
}

#' @importFrom stats rpois
#' @importFrom distr UnivarMixingDistribution
#' @importFrom distr Unif
#' @importFrom distr r
generate_data_simple <- function(n, p, K, link = c("log", "log1p")) {
  if (!is.scalar(K) || K < 1)
    stop("\"K\" must be an integer greater than or equal to 1")
  if (!is.scalar(n) || n < 1)
    stop("\"n\" must be an integer greater than or equal to 1")
  if (!is.scalar(p) || p < 1)    
    stop("\"p\" must be an integer greater than or equal to 1")
  
  link <- match.arg(link)
  LL <- matrix(nrow = K, ncol = n)
  FF <- matrix(nrow = K, ncol = p)
  
  l_dist <- UnivarMixingDistribution(
    Unif(0,0.01),
    Unif(.25,0.5),
    Unif(.5,0.75),
    mixCoeff = rep(1/3,3)
  )
  
  f_dist <- UnivarMixingDistribution(
    Unif(0,0.01),
    Unif(.25,0.5),
    Unif(.5,0.75),
    mixCoeff = rep(1/3,3)
  )
  
  l_sampler <- distr::r(l_dist)
  f_sampler <- distr::r(f_dist)
  
  for (k in 1:K) {
    LL[k, ] <- l_sampler(n)
    FF[k, ] <- f_sampler(p)
  }
  
  if (link == "log")
    Lambda <- exp(crossprod(LL,FF)) 
  else if (link == "log1p")
    Lambda <- exp(crossprod(LL,FF)) - 1
  
  Y_dat <- rpois(n = n * p,lambda = as.vector(Lambda))
  Y <- matrix(data = Y_dat, nrow = n, ncol = p)
  return(list(Y = Y,LL = LL,FF = FF))
}

