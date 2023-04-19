command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])
n_iter = as.integer(command_args[2])
n_cores = as.integer(command_args[3])

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")

data <- as.matrix(counts)

set.seed(1)
fit0 <- plash::init_glmpca(
  Y = data, K = n_factor, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
)

tic()
fit <- plash::fit_glmpca(
  Y = data, 
  fit0 = fit0, 
  algorithm = "ccd", 
  link = "log",
  control = list(line_search = TRUE, num_iter = 3, alpha = .1),
  warmup = TRUE, 
  max_iter = n_iter,
  tol = .Machine$double.eps
)
toc()

readr::write_rds(
  fit, 
  glue::glue("droplets_plash_fit_{n_cores}_cores_{n_factor}_factors_{n_iter}_iter.rds")
)

