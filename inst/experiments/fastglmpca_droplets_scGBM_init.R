command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])
n_iter = as.integer(command_args[2])
n_cores = as.integer(command_args[3])

library(Matrix)
library(RhpcBLASctl)
library(scGBM)
blas_set_num_threads(1)
omp_set_num_threads(n_cores)

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")

set.seed(1)

# first, get an initial fit from scGBM
scGBM_fit0 <- scGBM::gbm.sc(
  Y = as.matrix(Matrix::t(counts)),
  M = n_factor,
  max.iter = 1,
  min.iter = 1,
  return.W = FALSE,
  infer.beta = FALSE
)

fit0 <- fastglmpca::init_glmpca_pois(
  Y = Matrix::t(counts), K = n_factor, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
)

# add intercepts
fit0$LL[2, ] <- scGBM_fit0$alpha
fit0$FF[1, ] <- scGBM_fit0$beta

for (k in 1:n_factor) {
  
  fit0$LL[k+2, ] <- scGBM_fit0$U[,k]
  fit0$FF[k+2, ] <- scGBM_fit0$V[,k]
  
}

rm(scGBM_fit0)
gc()

library(tictoc)
tic()
fit <- fastglmpca::fit_glmpca_pois(
  Y = Matrix::t(counts),
  fit0 = fit0,
  control = list(line_search = TRUE, num_iter = 3, alpha = .01),
  max_iter = n_iter,
  tol = .Machine$double.eps
)
toc()

readr::write_rds(
  fit,
  glue::glue("droplets_fastglmpca_fit_{n_cores}_cores_{n_factor}_factors_{n_iter}_iter_scGBM_init.rds")
)