# sinteractive -p mstephens --account=pi-mstephens --time=24:00:00 \
#   -c 10 --mem=16G
line_search    <- TRUE
use_daarem     <- FALSE
orthonormalize <- TRUE
library(Matrix)
library(RhpcBLASctl)
library(fastglmpca)
load(file.path("/project2/mstephens/pcarbo/git/fastTopics-experiments/data",
               "droplet.RData"))
omp_set_num_threads(10)
blas_set_num_threads(1)
set.seed(1)
counts <- t(counts)
cat("line_search    =",line_search,"\n")
cat("use_daarem     =",use_daarem,"\n")
cat("orthonormalize =",orthonormalize,"\n")
fit0 <- fit_glmpca_pois(counts,K = 5,
                        control = list(maxiter = 10,tol = 0,
                                       line_search = TRUE,
                                       use_daarem = FALSE,
                                       orthonormalize = TRUE))
fit <- fit_glmpca_pois(counts,fit0 = fit0,
                       control = list(maxiter = 10,tol = 0,
                                      line_search = line_search,
                                      use_daarem = use_daarem,
                                      orthonormalize = orthonormalize))
save(list = "fit",file = "droplet_k=5_daarem=no.RData")
