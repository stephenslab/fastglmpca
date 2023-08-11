# library(fastglmpca)
library(ggplot2)
library(cowplot)
set.seed(1)
data(pbmc_facs)
Y <- as.matrix(pbmc_facs$counts)
fit0 <- init_glmpca_pois(Y,K = 3)
stop()
crossprod(fit0$U)
crossprod(fit0$V)
fit0$D
fit <- fit_glmpca_pois(Y,K = 3,max_iter = 4)

