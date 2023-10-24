# TO DO: Explain here what this script is for, and how to use it.
#
#   source activate base
#   export MEM_CHECK_INTERVAL=0.01
#   python3 monitor_memory.py Rscript sparse_benchmark.R
#
# library(fastglmpca)
devtools::load_all()
library(Matrix)
load("../datasets/newsgroups.RData")
counts <- counts[1:6000,1:2e4]
topics <- topics[1:6000]
# counts <- as.matrix(counts)
print(format(object.size(counts),unit = "Mb"))
t0 <- proc.time()
fit <- fit_glmpca_pois(counts,K = 3,max_iter = 10,
                       control = list(calc_max_diff = TRUE,calc_deriv = TRUE))
t1 <- proc.time()
print(t1 - t0)
