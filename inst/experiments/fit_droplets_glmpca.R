command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])
n_iter = as.integer(command_args[2])
n_cores = as.integer(command_args[3])

library(Matrix)
library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(n_cores)
load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplets.RData")


scGBM_fit0 <- readr::read_rds(
  glue::glue("/home/ericweine/glmpca_experiments/droplets_scGBM_fit_{n_factor}_factors_init_1_iter_no_beta_infer.rds")
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
fit <- glmpca::glmpca(
  Y = Matrix::t(counts),
  L = n_factor,
  fam = "poi",
  optimizer = "avagrad",
  minibatch = "stochastic",
  ctl = list(minIter = n_iter - 1, maxIter = n_iter, verbose = TRUE, tol = .Machine$double.eps),
  init = list(factors = t(fit0$FF[-c(1,2),]), loadings = t(fit0$LL[-c(1,2),])),
  sz = exp(fit0$FF[1, ]),
  init_coefX = matrix(data = fit0$LL[2, ], ncol = 1)
)

toc()

readr::write_rds(
  fit,
  glue::glue(
    "droplets_glmpca_fit_{n_factor}_factors_{n_iter}_iter_avagrad_optimizer_minibatch_stochastic.rds"
  )
)