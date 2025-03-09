library(fastglmpca)

set.seed(1)
dat <- generate_glmpca_data_pois(500, 250, 2)

set.seed(1)
fit1 <- fit_glmpca_pois(
  Y = dat$Y, 
  K = 2,
  control = list(training_frac = 0.99, maxiter = 100)
)

set.seed(1)
fit2 <- fit_glmpca_pois(
  Y = dat$Y, 
  K = 2,
  control = list(training_frac = 1, maxiter = 100)
)
