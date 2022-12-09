# Here, I want to test some of the code I have built

generate_pois_flash <- function(n, p, K, l_dist_list, f_dist_list) {

  LL <- matrix(nrow = K, ncol = n)
  FF <- matrix(nrow = K, ncol = p)

  for (k in 1:K) {

    l_sampler <- distr::r(l_dist_list[[k]])
    f_sampler <- distr::r(f_dist_list[[k]])

    LL[k, ] <- l_sampler(n)
    FF[k, ] <- f_sampler(p)

  }

  Z <- t(LL) %*% FF

  Lambda <- exp(Z)

  Y <- matrix(
    data = rpois(n * p, as.vector(Lambda)), nrow = n, ncol = p
  )

  return(
    list(
      LL = LL, FF = FF, Z = Z, Lambda = Lambda, Y = Y
    )
  )

}

l_dist_list <- list()
f_dist_list <- list()

for (i in 1:5) {

  l_dist_list[[i]] <- distr::UnivarMixingDistribution(
    distr::Dirac(0),
    distr::Unif(0, 1),
    mixCoeff = c(.2, .8)
  )

  f_dist_list[[i]] <- distr::UnivarMixingDistribution(
    distr::Dirac(0),
    distr::Dirac(i/2.5),
    mixCoeff = c((1 / 2), (1 / 2))
  )

}

p_out <- generate_pois_flash(2500, 1000, 5, l_dist_list, f_dist_list)

glmpca_out <- glmpca(p_out$Y,L = 5,optimizer = "fisher",
              ctl = list(minIter = 2,maxIter = 15,tol = 1e-8,
                         verbose = TRUE))

n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

# this should give essentially the same fit as glmpca
plash_fit_L_c_sep_c_first <- plash::plash_omni(
  Y = p_out$Y, K = 5, offset = TRUE, intercept = TRUE, tol = 5e-5, parallel = TRUE,
  init_cc = rep(0, nrow(p_out$Y)), update_c = TRUE, update_c_first = TRUE
)

plash_fit_L_c_sep_c_last <- plash::plash_omni(
  Y = p_out$Y, K = 5, offset = TRUE, intercept = TRUE, tol = 5e-5, parallel = TRUE,
  init_cc = rep(0, nrow(p_out$Y)), update_c = TRUE, update_c_first = FALSE
)

plash_fit_L_c_siuml <- plash::plash_omni(
  Y = p_out$Y, K = 5, offset = TRUE, intercept = TRUE, tol = 5e-5, parallel = TRUE,
  init_cc = rep(0, nrow(p_out$Y)), update_c = TRUE, update_L_c_simul = TRUE
)

# Now, I will attempt to see if the code works when I update L and c simultaneously
# I have the least confidence in this algorithm, because I feel like the construction of the 
# jacobian and gradient are a bit sketchy

plash_fit_L_c_simul <- plash::plash_omni(
  Y = p_out$Y, K = 5, offset = TRUE, intercept = TRUE, tol = 1e-5, parallel = TRUE,
  update_L_c_simul = TRUE
)

# Now, fit the model with my code
plash_fit <- plash::plash_parallel_glmpca(
  Y = p_out$Y, K = 5, cluster = my.cluster, offset = TRUE, offset_vec = glmpca_out$offsets,
  intercept = TRUE, tol = 1e-6
)

p_v2_test <- plash::plash_parallel_optim_v2(
  Y = p_out$Y, K = 5, update_L = T, update_F = T, update_c = T, tol = 1e-4, cluster = my.cluster
)

tic()
p_test <- plash::plash(
  Y = p_out$Y, K = 5, true_LL = p_out$LL, true_FF = p_out$FF, true_cc = rep(0, 1000),
  update_L = T, update_F = T, update_c = T, tol = 1e-4
)
toc()

tic()
# note, before running this I must transform this into a matrix
# and then I also want to take the transpose
p_test_par <- plash::plash_parallel(
  Y = p_out$Y, cluster = my.cluster, K = 5, update_L = T, update_F = T, update_c = F, tol = 1e-6,
  true_cc = rep(0, 2500),
)
toc()

tic()
glmpca_test <- fit_glmpca(
  X = p_out$Y, k = 5, numiter = 6, method = "fastglm"
)
toc()



fl_out <- flashier::flash(
  data = log(p_out$Y + 1), greedy.Kmax = 5, backfit = T
)

# Here, need to test plash omni
# Once I do this I will have all of the algorithms in one place, and it will be much easier to test 
# different options
# this will also help me clean up the package a bit, FWIW
