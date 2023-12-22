command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])
n_iter = as.integer(command_args[2])

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")

library(RhpcBLASctl)
omp_set_num_threads(28)
blas_set_num_threads(1)

set.seed(1)
fit0 <- fastglmpca::init_glmpca_pois(
  Y = Matrix::t(counts),
  K = n_factor
)

fit <- fastglmpca::fit_glmpca_pois(
  Y = Matrix::t(counts),
  fit0 = fit0,
  control = list(
    use_daarem = FALSE,
    maxiter = n_iter,
    tol = 0,
    mon.tol = 0.05,
    convtype = "objfn",
    line_search = TRUE,
    ls_alpha = 0.001,
    ls_beta = 0.5,
    num_ccd_iter = 3,
    ccd_iter_tol = 0,
    calc_deriv = FALSE,
    calc_max_diff = FALSE,
    orthonormalize = TRUE
  )
)

readr::write_rds(
  fit,
  glue::glue("droplets_fastglmpca_fit_{n_factor}_factors_{n_iter}_iter_28_cores_dec_23.rds")
)
