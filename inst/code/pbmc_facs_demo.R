library(ggplot2)
library(cowplot)
set.seed(1)
data(pbmc_facs)
Y <- pbmc_facs$counts
fit <- fit_glmpca_pois(Y,K = 2)
