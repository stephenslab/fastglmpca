command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])
n_iter = as.integer(command_args[2])
n_cores = as.integer(command_args[3])

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")
#load("~/Documents/plash/scratch/droplet.RData")

data <- as.matrix(counts)

set.seed(1)
fit0 <- plash::init_glmpca(
  Y = data, K = n_factor, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
)

library(tictoc)
tic()
fit_init <- plash::fit_glmpca(
  Y = data, 
  fit0 = fit0, 
  algorithm = "ccd", 
  link = "log",
  control = list(line_search = TRUE, num_iter = 3, alpha = .01),
  warmup = FALSE, 
  tol = .Machine$double.eps,
  use_daarem = FALSE,
  max_iter = 3
)
toc()

readr::write_rds(
  fit, 
  glue::glue("droplets_plash_fit_init_{n_cores}_cores_{n_factor}_factors_{n_iter}_iter.rds")
)

tic()
fit_daarem <- plash::fit_glmpca(
  Y = data, 
  fit0 = fit_init, 
  algorithm = "ccd", 
  link = "log",
  control = list(line_search = TRUE, num_iter = 3, alpha = .01),
  warmup = FALSE, 
  tol = .Machine$double.eps,
  use_daarem = TRUE,
  max_iter = n_iter - 3
)
toc()

readr::write_rds(
  fit, 
  glue::glue("droplets_plash_fit_daarem_{n_cores}_cores_{n_factor}_factors_{n_iter}_iter.rds")
)

