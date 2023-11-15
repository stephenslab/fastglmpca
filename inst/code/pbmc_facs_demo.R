# library(fastglmpca)
library(ggplot2)
library(cowplot)
set.seed(1)
data(pbmc_facs)
Y <- as.matrix(pbmc_facs$counts)
# n <- nrow(Y)
# m <- ncol(Y)
# X <- matrix(rnorm(2*n),n,2)
# Z <- matrix(rnorm(m),m,1)
set.seed(1)
fit0 <- init_glmpca_pois(Y,K = 3)
# fit0_init <- init_glmpca_pois(Y,X = X,Z = Z,
#                               U = matrix(rnorm(3*n),n,3),
#                               V = matrix(rnorm(3*m),m,3))
# fit0_rank1 <- init_glmpca_pois(Y,X = X,Z = Z,K = 1)
fit0 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(use_daarem = FALSE,
                                       maxiter = 4,
                                       calc_max_diff = TRUE,
                                       calc_deriv = TRUE,
                                       orthonormalize = TRUE))
fit1 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(use_daarem = FALSE,
                                       maxiter = 40,
                                       calc_max_diff = TRUE,
                                       calc_deriv = TRUE,
                                       orthonormalize = FALSE))
fit2 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(use_daarem = TRUE,
                                       maxiter = 40,
                                       calc_max_diff = TRUE,
                                       calc_deriv = TRUE,
                                       orthonormalize = FALSE))
pdat <- rbind(data.frame(method = "fpiter",
                         iter   = fit1$progress$iter,
                         loglik = fit1$progress$loglik),
              data.frame(method = "daarem",
                         iter   = fit2$progress$iter,
                         loglik = fit2$progress$loglik))
best_loglik <- max(pdat$loglik)
pdat <- transform(pdat,loglik = best_loglik - loglik + 0.01)
ggplot(pdat,aes(x = iter,y = loglik,color = method)) +
  geom_line(size = 0.5) +
  scale_color_manual(values = c("darkblue","darkorange"))+
  scale_y_continuous(trans = "log10") +
  theme_cowplot(font_size = 12)
