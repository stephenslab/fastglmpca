# library(fastglmpca)
library(ggplot2)
library(cowplot)
set.seed(1)
data(pbmc_facs)
Y_sparse <- pbmc_facs$counts
Y <- as.matrix(Y_sparse)
n <- nrow(Y)
m <- ncol(Y)
X <- matrix(rnorm(2*n),n,2)
Z <- matrix(rnorm(m),m,1)
set.seed(1)
fit0 <- init_glmpca_pois(Y,X = X,Z = Z,K = 3)
set.seed(1)
fit0_sparse <- init_glmpca_pois(Y_sparse,X = X,Z = Z,K = 3)
print(fit0$loglik - fit0_sparse$loglik)
fit0_init <- init_glmpca_pois(Y,X = X,Z = Z,
                              U = matrix(rnorm(3*n),n,3),
                              V = matrix(rnorm(3*m),m,3))
fit0_rank1 <- init_glmpca_pois(Y,X = X,Z = Z,K = 1)

# Using a single thread.
set_fastglmpca_threads(1)
t0 <- proc.time()
fit1 <- fit_glmpca_pois(Y,fit0 = fit0,max_iter = 10,
                        control = list(calc_max_diff = FALSE,
                                       calc_deriv = FALSE))
t1 <- proc.time()
print(t1 - t0)                

# Using 2 threads.
set_fastglmpca_threads(2)
t0 <- proc.time()
fit2 <- fit_glmpca_pois(Y,fit0 = fit0,max_iter = 10,
                        control = list(calc_max_diff = FALSE,
                                       calc_deriv = FALSE))
t1 <- proc.time()
print(t1 - t0)                
