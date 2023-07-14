command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplets.RData")

library(scGBM)

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

readr::write_rds(
  scGBM_fit0,
  glue::glue("droplets_scGBM_fit_{n_factor}_factors_init_1_iter_no_beta_infer.rds")
)