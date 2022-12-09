# code to simulate plash data here
# I'd like a realistic model that is sparse and represents 
# scRNA data vaguely

generate_plash_data <- function(n, p, K) {
  
  LL <- matrix(nrow = K, ncol = n)
  FF <- matrix(nrow = K, ncol = p)
  
  l_dist <- distr::UnivarMixingDistribution(
    distr::Dirac(0),
    distr::Dirac(1),
    distr::Dirac(2),
    mixCoeff = rep((1/3), 3)
  )
  
  f_dist <- distr::UnivarMixingDistribution(
    distr::Dirac(0),
    distr::Dirac(1),
    distr::Dirac(2),
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
  
  c_sampler <- distr::r(c_dist)
  
  cc <- c_sampler(n)
  
  cc <- 0
  
  Lambda <- exp(crossprod(LL, FF)) - cc
  
  Pi <- sweep(Lambda, 2, colSums(Lambda), `/`)
  
  Y <- matrix(nrow = n, ncol = p)
  
  for (j in 1:p) {
    
    Y[, j] <- rmultinom(n = 1, size = cell_totals[j], prob = Pi[, j])
    
  }
  
  return(
    list(
      Y = Y, LL = LL, FF = FF, cc = cc
    )
  )
  
}

