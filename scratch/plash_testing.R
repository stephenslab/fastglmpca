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

p_out <- generate_pois_flash(2500, 2500, 5, l_dist_list, f_dist_list)

n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

tic()
p_test <- plash::plash(
  Y = p_out$Y, K = 5, true_LL = p_out$LL, true_FF = p_out$FF, true_cc = rep(0, 1000),
  update_L = T, update_F = T, update_c = T, tol = 1e-4
)
toc()

tic()
p_test_par <- plash::plash_parallel(
  Y = p_out$Y, cluster = my.cluster, K = 5, true_LL = p_out$LL, true_FF = p_out$FF, true_cc = rep(0, 1000),
  update_L = T, update_F = T, update_c = T, tol = 1e-4
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
