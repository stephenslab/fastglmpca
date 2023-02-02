# I think perhaps what is happening is that there is an issue when you fit an intercept
# I want to go through a few more cases without an intercept just to make sure things work
# then, once I am confident in this, I can move on to the case where there is an intercept
set.seed(1)
data <- plash:::generate_plash_data_simple(2500, 1250, 5)

library(tictoc)
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)

# note that before I do speed comparisons
# it may make sense to have a control parameter for how many iterations of bfgs we do
# because I may be "over-optimizing"
tic()
plash_fit_bfgs <- plash:::plash_glmpca(
  Y = data$Y, K = 5, offset = FALSE, intercept = FALSE, tol = 1e-5,
  update_L = TRUE, update_F = TRUE, init_FF = NULL, algorithm = "bfgs",
  parallel = TRUE, bfgs_ctl = list(maxit = 50)
)
toc()

tic()
plash_fit_ccd <- plash:::plash_glmpca(
  Y = data$Y, K = 5, offset = TRUE, intercept = TRUE, tol = 1e-5,
  update_L = TRUE, update_F = TRUE, init_FF = NULL, algorithm = "ccd",
  parallel = FALSE, ccd_ctl = list(num_iter = 10, alpha = .25, beta = .5, line_search = TRUE)
)
toc()
