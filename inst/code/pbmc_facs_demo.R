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
fit0 <- init_glmpca_pois(Y,K = 3,col_size_factor = FALSE,row_intercept = FALSE)
fit0 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(maxiter = 4,
                                       orthonormalize = TRUE,
                                       daarem = FALSE))
# fit0_init <- init_glmpca_pois(Y,X = X,Z = Z,
#                               U = matrix(rnorm(3*n),n,3),
#                               V = matrix(rnorm(3*m),m,3))
# fit0_rank1 <- init_glmpca_pois(Y,X = X,Z = Z,K = 1)
fit1 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(use_daarem = FALSE,
                                       maxiter = 100,
                                       calc_max_diff = FALSE,
                                       calc_deriv = FALSE,
                                       orthonormalize = FALSE))
fit2 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(use_daarem = TRUE,
                                       maxiter = 100,
                                       calc_max_diff = FALSE,
                                       calc_deriv = FALSE,
                                       orthonormalize = FALSE))
fit3 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(use_daarem = TRUE,
                                       maxiter = 100,
                                       calc_max_diff = FALSE,
                                       calc_deriv = FALSE,
                                       orthonormalize = TRUE))
pdat <- rbind(data.frame(method = "fpiter",
                         iter   = fit1$progress$iter,
                         loglik = fit1$progress$loglik),
              data.frame(method = "daarem",
                         iter   = fit2$progress$iter,
                         loglik = fit2$progress$loglik),
              data.frame(method = "daarem+ortho",
                         iter   = fit3$progress$iter,
                         loglik = fit3$progress$loglik))
best_loglik <- max(pdat$loglik)
pdat <- transform(pdat,loglik = best_loglik - loglik + 1)
p <- ggplot(pdat,aes(x = iter,y = loglik,color = method)) +
  geom_line(size = 0.5) +
  scale_color_manual(values = c("darkblue","darkorange","magenta"))+
  scale_y_continuous(trans = "log10") +
  theme_cowplot(font_size = 12)
print(p)
