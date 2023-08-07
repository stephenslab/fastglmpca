# library(fastglmpca)
library(ggplot2)
library(cowplot)
set.seed(1)
data(pbmc_facs)
Y <- pbmc_facs$counts
fit <- fit_glmpca_pois(Y,K = 3,max_iter = 4)
range(with(fit,crossprod(LL,FF)) - with(fit2,tcrossprod(U %*% D,V)))
