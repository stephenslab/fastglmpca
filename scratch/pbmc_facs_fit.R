# #n.cores <- parallel::detectCores() - 1
# set.seed(1)
# n.cores <- 11
# my.cluster <- parallel::makeCluster(
#   n.cores,
#   type = "PSOCK"
# )
# 
# # first, obtain a glmpca fit
# glmpca_out <- glmpca(Matrix::t(fastTopics::pbmc_facs$counts),L = 6,optimizer = "fisher",
#                      ctl = list(minIter = 2,maxIter = 35,tol = 1e-5,
#                                 verbose = TRUE))
# 
# count_mat <- t(as.matrix(fastTopics::pbmc_facs$counts))

n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

pbmc_plash_fit_L_c_sep_c_last <- plash::plash_omni(
  Y = Matrix::t(fastTopics::pbmc_facs$counts), K = 6, offset = TRUE, intercept = TRUE, tol = 1e-4, parallel = TRUE, update_c = TRUE,
  update_c_first = FALSE
)

pbmc_plash_fit_L_c_simul <- plash::plash_omni(
  Y = count_mat, K = 6, offset = TRUE, intercept = TRUE, tol = 1e-4, parallel = TRUE, update_c = TRUE,
  update_L_c_simul = TRUE
)

pbmc_plash_fit_L_c_sep_c_first <- plash::plash_omni(
  Y = count_mat, K = 6, offset = TRUE, intercept = TRUE, tol = 1e-4, parallel = TRUE, update_c = TRUE,
  update_c_first = TRUE
)

plash_fit_c_v2 <- plash::plash_parallel_optim_v2(
  Y = count_mat, cluster = my.cluster, K = 6, tol = 1e-4
)

plash_fit_no_c <- plash::plash_parallel(
  Y = count_mat, cluster = my.cluster, K = 6, true_cc = rep(0, nrow(count_mat)), update_c = FALSE, tol = 1e-3
)

plash_fit_c_1 <- plash::plash_parallel(
  Y = count_mat, cluster = my.cluster, K = 6, true_cc = rep(1, nrow(count_mat)), update_c = FALSE, tol = 1e-3
)

plash_fit_no_c <- plash::plash_parallel(
  Y = count_mat, cluster = my.cluster, update_c = F, true_cc = rep(0, nrow(count_mat)), K = 6
)
