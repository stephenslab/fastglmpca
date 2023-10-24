# TO DO: Explain here what this script is for, and how to use it.
# library(fastglmpca)
devtools::load_all()
library(Matrix)
load("../datasets/newsgroups.RData")
counts <- counts[1:4000,1:1e4]
topics <- topics[1:4000]
# counts <- as.matrix(counts)
t0 <- proc.time()
fit <- fit_glmpca_pois(counts,K = 3,max_iter = 10,
                       control = list(calc_max_diff = TRUE,
                                       calc_deriv = TRUE))
t1 <- proc.time()
print(t1 - t0)
