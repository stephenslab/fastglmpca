# sinteractive -p mstephens --account=pi-mstephens --time=24:00:00 \
#   -c 10 --mem=16G
# use_daarem     <- FALSE
# orthonormalize <- TRUE
# outfile        <- "droplet_k=5_daarem=no.RData"
use_daarem     <- TRUE
orthonormalize <- FALSE
outfile        <- "droplet_k=5_daarem=yes.RData"
library(Matrix)
library(RhpcBLASctl)
library(fastglmpca)
load(file.path("/project2/mstephens/pcarbo/git/fastTopics-experiments/data",
               "droplet.RData"))
omp_set_num_threads(10)
blas_set_num_threads(1)
set.seed(1)
counts <- t(counts)
cat("use_daarem     =",use_daarem,"\n")
cat("orthonormalize =",orthonormalize,"\n")
cat("outfile        =",outfile,"\n")
fit0 <- fit_glmpca_pois(counts,K = 5,
                        control = list(maxiter = 100,tol = 0,
                                       use_daarem = FALSE,
                                       orthonormalize = TRUE))
fit <- fit_glmpca_pois(counts,fit0 = fit0,
                       control = list(maxiter = 100,tol = 0,
                                      use_daarem = use_daarem,
                                      orthonormalize = orthonormalize))
save(list = "fit",file = outfile)
