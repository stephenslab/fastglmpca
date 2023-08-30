# library(fastglmpca)
library(ggplot2)
library(cowplot)
set.seed(1)
data(pbmc_facs)
Y <- as.matrix(pbmc_facs$counts)
n <- nrow(Y)
m <- ncol(Y)
X <- matrix(rnorm(2*n),n,2)
Z <- matrix(rnorm(m),m,1)
fit0 <- init_glmpca_pois(Y,X = X,Z = Z,K = 3)
fit0_init <- init_glmpca_pois(Y,X = X,Z = Z,
                              U = matrix(rnorm(3*n),n,3),
                              V = matrix(rnorm(3*m),m,3))
fit0_rank1 <- init_glmpca_pois(Y,X = X,Z = Z,K = 1)
# fit <- fit_glmpca_pois(Y,K = 3,max_iter = 4)

