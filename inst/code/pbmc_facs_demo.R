# library(fastglmpca)
library(ggplot2)
library(cowplot)
set.seed(1)
data(pbmc_facs)
Y <- as.matrix(pbmc_facs$counts)
n <- nrow(Y)
m <- ncol(Y)
fit0_rank1 <- init_glmpca_pois(Y,K = 1)
fit0_init <- init_glmpca_pois(Y,
                              U = matrix(rnorm(3*n),n,3),
                              V = matrix(rnorm(3*m),m,3))
fit0 <- init_glmpca_pois(Y,K = 3)
# fit <- fit_glmpca_pois(Y,K = 3,max_iter = 4)

