# library(fastglmpca)
library(ggplot2)
library(cowplot)
set.seed(1)
data(pbmc_facs)
Y <- as.matrix(pbmc_facs$counts)
<<<<<<< HEAD
fit0 <- init_glmpca_pois(Y,K = 3)
stop()
=======
fit0 <- init_glmpca_pois(Y,K = 3,
                         fit_col_size_factor = TRUE,
                         fit_row_intercept = TRUE)
>>>>>>> 85424f4d3e05969f6e6d0897c3ce0cc2304e3474
crossprod(fit0$U)
crossprod(fit0$V)
fit0$D
fit <- fit_glmpca_pois(Y,K = 3,max_iter = 4)

