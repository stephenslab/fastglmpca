data <- plash:::generate_plash_data_simple(8000, 4000, 15)

# Now, I want to see how well the algorithms I've built work on solving
# some very simple problems. 

# First, I'd like to see if I get the same solution for fitting L with F fixed

# Now, see if I can fit FF in addition to LL

# It's pretty clear to me that this is at least a correct implementation
# The issue is now can I make the algorithm faster by using a line search
# this will hopefully prevent the initial crazy increases in the objective
# which will increase the cost per iteration but likely massively
# decrease the total number of iterations required

tic()
plash_fit2 <- plash:::fit_plash_ccd(
  Y = data$Y, K = 15, update_LL = FALSE, update_FF = TRUE, init_LL = data$LL, 
  init_FF = NULL, min_iter = 5, max_iter = 30, tol = 1e-5, line_search = TRUE
)
toc()

n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)

tic()
plash_fit_glmp <- plash:::plash_glmpca(
  Y = data$Y, K = 15, offset = FALSE, intercept = FALSE, tol = 1e-5,
  update_L = FALSE, update_F = TRUE, init_LL = data$LL
)
toc()

tic()
glmpca_fit <- glmpca::glmpca(
  Y = data$Y, L = 15, optimizer = "fisher", ctl = list(verbose = TRUE)
)
toc()
