command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])
n_iter = as.integer(command_args[2])

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")

data <- as.matrix(Matrix::t(counts))

set.seed(1)

library(tictoc)

tic()
fit <- scGBM::gbm.sc(
  Y = data,
  M = n_factor,
  max.iter = n_iter,
  tol = .Machine$double.eps,
  time.by.iter = TRUE,
  infer.beta = FALSE,
  return.W = FALSE,
  min.iter = n_iter - 1
)
toc()

readr::write_rds(
  fit,
  glue::glue("droplets_scGBM_fit_{n_factor}_factors_{n_iter}_iter_no_beta_infer.rds")
)