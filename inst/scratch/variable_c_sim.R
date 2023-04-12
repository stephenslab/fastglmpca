# Here, want to run a simulation where I fit glmpca on a grid of c values
set.seed(2020)

get_lik_nmf <- function(Y, ft_fit) {
  
  Lambda <- tcrossprod(ft_fit$L, ft_fit$F)
  lik <- dpois(as.vector(Y), as.vector(Lambda), log = TRUE)
  return(mean(lik))
  
}

get_lik_glmpca <- function(Y, glmpca_fit) {
  
  Lambda <- exp(crossprod(glmpca_fit$LL, glmpca_fit$FF)) - outer(glmpca_fit$cc, glmpca_fit$size)
  lik <- dpois(as.vector(Y), as.vector(Lambda), log = TRUE)
  return(mean(lik))
  
}

n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

nmf_1f_data <- plash:::generate_1_factor_nmf_data(750, 1000)

nmf_fit <- fastTopics:::fit_pnmf_rank1(nmf_1f_data$Y_train)
nmf_ll <- get_lik_nmf(nmf_1f_data$Y_train, nmf_fit)

c_vec <- seq(0, 2.5, .1)
plash_ll_vec <- numeric(length(c_vec))

init <- plash:::get_feasible_init(
  Y = t(nmf_1f_data$Y_train),
  K = 1,
  offset = FALSE,
  intercept = FALSE,
  cc = rep(2.6, 1000)
)

for (i in 1:length(c_vec)) {
  
  print(glue::glue("Running simulations for c = {c_vec[i]}"))
  
  plash_mod <- plash::plash_omni(
    Y = t(nmf_1f_data$Y_train),
    K = 1,
    init_cc = rep(c_vec[i], 1000),
    update_c = FALSE,
    offset = FALSE,
    intercept = FALSE,
    tol = 5e-5,
    init_LL = init$LL,
    init_FF = init$FF,
    max_iter = 20,
    parallel = TRUE
  )
  
  plash_ll <- get_lik_glmpca(t(nmf_1f_data$Y_train), plash_mod)
  plash_ll_vec[i] <- plash_ll
  
}

sim_df <- data.frame(cc = c_vec, plash_ll = plash_ll_vec)
#sim_df <- data.frame(cc = c_vec[-26], plash_ll = plash_ll_vec[-26])

ggplot(data = sim_df, aes(x = cc, y = plash_ll)) +
  geom_point() +
  xlab("c") +
  ylab("GLMPCA LogLik") +
  geom_hline(yintercept=nmf_ll, linetype="dashed", color = "red") +
  annotate("text", x = .375, y = nmf_ll + .05, label = "fastTopics LogLik", size = 3.5, colour = "red") +
  ylim(-1.05, nmf_ll + .075)

plash_mod_fin <- plash::plash_omni(
  Y = t(nmf_1f_data$Y_train),
  K = 1,
  init_cc = rep(0, 1000),
  update_c = T,
  offset = FALSE,
  intercept = FALSE,
  tol = 5e-5,
  max_iter = 25,
  parallel = TRUE,
  update_L_c_simul = T
)
