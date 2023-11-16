# library(fastglmpca)
library(ggplot2)
library(cowplot)
set.seed(1)
data(pbmc_facs)
Y <- pbmc_facs$counts
set.seed(1)
fit0 <- init_glmpca_pois(Y,K = 3)
fit0 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(maxiter = 4,
                                       orthonormalize = TRUE,
                                       daarem = FALSE))
fit1 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(use_daarem = FALSE,
                                       maxiter = 100,
                                       calc_max_diff = TRUE,
                                       calc_deriv = TRUE,
                                       orthonormalize = TRUE))
fit2 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(use_daarem = TRUE,
                                       maxiter = 100,
                                       calc_max_diff = TRUE,
                                       calc_deriv = TRUE,
                                       orthonormalize = FALSE))
fit3 <- fit_glmpca_pois(Y,fit0 = fit0,
                        control = list(use_daarem = TRUE,
                                       maxiter = 100,
                                       calc_max_diff = TRUE,
                                       calc_deriv = TRUE,
                                       orthonormalize = TRUE))
pdat <- rbind(data.frame(method = "fpiter",
                         iter   = fit1$progress$iter,
                         loglik = fit1$progress$loglik),
              data.frame(method = "daarem-ortho",
                         iter   = fit2$progress$iter,
                         loglik = fit2$progress$loglik),
              data.frame(method = "daarem+ortho",
                         iter   = fit3$progress$iter,
                         loglik = fit3$progress$loglik))
best_loglik <- max(pdat$loglik)
pdat <- transform(pdat,loglik = best_loglik - loglik + 100)
p <- ggplot(pdat,aes(x = iter,y = loglik,color = method)) +
  geom_line(size = 0.5) +
  scale_color_manual(values = c("darkblue","darkorange","magenta"))+
  scale_y_continuous(trans = "log10",breaks = c(1,10,100,1e3,1e4,1e5,1e6)) +
  theme_cowplot(font_size = 12)
print(p)
