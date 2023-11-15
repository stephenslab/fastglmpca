# library(fastglmpca)
library(ggplot2)
library(cowplot)
set.seed(1)
data(pbmc_facs)
Y <- pbmc_facs$counts
n <- nrow(Y)
m <- ncol(Y)
X <- matrix(rnorm(2*n),n,2)
Z <- matrix(rnorm(m),m,1)
set.seed(1)
fit0 <- init_glmpca_pois(Y,X = X,Z = Z,K = 3)
fit0_init <- init_glmpca_pois(Y,X = X,Z = Z,
                              U = matrix(rnorm(3*n),n,3),
                              V = matrix(rnorm(3*m),m,3))
fit0_rank1 <- init_glmpca_pois(Y,X = X,Z = Z,K = 1)
fit1 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(use_daarem = FALSE,
                                       maxiter = 50,
                                       calc_max_diff = TRUE,
                                       calc_deriv = TRUE,
                                       orthonormalize = TRUE))
fit2 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(use_daarem = TRUE,
                                       maxiter = 50,
                                       calc_max_diff = TRUE,
                                       calc_deriv = TRUE))

# TO DO: Add plot comparing progress of fit1 and fit2.
