# library(fastglmpca)
library(ggplot2)
library(cowplot)
set.seed(1)
data(pbmc_facs)
Y <- as.matrix(pbmc_facs$counts)
set.seed(1)
fit0 <- fit_glmpca_pois(Y,K = 3,max_iter = 4,
                        control = list(orthonormalize = TRUE,daarem = FALSE))
t0 <- proc.time()
fit1 <- fit_glmpca_pois(Y,fit0 = fit0,max_iter = 40,
                        control = list(orthonormalize = TRUE,daarem = FALSE))
t1 <- proc.time()
print(t1 - t0)
t0 <- proc.time()
fit2 <- fit_glmpca_pois(Y,fit0 = fit0,max_iter = 40,
                        control = list(orthonormalize = FALSE,daarem = TRUE))
t1 <- proc.time()
print(t1 - t0)
pdat <- rbind(data.frame(method = "fpiter",
                         iter   = fit1$progress$iter,
                         loglik = fit1$progress$loglik),
              data.frame(method = "daarem",
                         iter   = fit2$progress$iter,
                         loglik = fit2$progress$loglik))
bestloglik <- max(fit1$loglik,fit2$loglik)
pdat <- subset(pdat,iter > 4)
pdat <- transform(pdat,loglik = bestloglik - loglik + 1)
p <- ggplot(pdat,aes(x = iter,y = loglik,color = method)) +
  geom_line(size = 0.75) +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = c("darkblue","darkorange")) +
  labs(y = "distance to best loglik") +
  theme_cowplot(font_size = 12)
print(p)
