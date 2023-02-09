# I think perhaps what is happening is that there is an issue when you fit an intercept
# I want to go through a few more cases without an intercept just to make sure things work
# then, once I am confident in this, I can move on to the case where there is an intercept
set.seed(1)
data <- plash:::generate_plash_data_simple(5000, 2500, 5)

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

LL <- matrix(
  data = .1, nrow = 5, ncol = 5000
)

FF <- matrix(
  data = .1, nrow = 5, ncol = 2500
)

#   # Offset
LL[1, ] <- 1
FF[1, ] <- log(colMeans(data$Y))
#   
#   # Intercept
FF[2, ] <- 1

tictoc::tic()
plash_fit_ccd <- plash:::fit_glmpca(
  Y = data$Y, K = 5, LL_update_indices = c(2:5), FF_update_indices = c(3:5), tol = 1e-5,
  update_L = TRUE, update_F = TRUE, init_FF = FF, init_LL = LL,
  parallel = FALSE, ccd_ctl = list(num_iter = 10, alpha = .25, beta = .5, line_search = TRUE)
)
tictoc::toc()

# I want to think about how to design an initialization method here
# I think that the init function should be delicate.
# I think one option should be just to have an intercept and cell specific offset
# and then if the user doesn't want that then they can just specify their own loadings and factors
