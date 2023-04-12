# TO DO: Describe here briefly what this script is for.

# remotes::install_github("eweine/glmpca",upgrade = "never")
library(glmpca)
library(ggplot2)
library(cowplot)

# Simulate a 200 x 500 counts matrix.
set.seed(1)
dat <- generate_glmpca_data(200,500,K = 3)
Y <- dat$Y

# Generate an initial fit of the Poisson GLM-PCA by running 8
# iterations of the glmpca "Fisher scoring" algorithm.
set.seed(1)
out <- glmpca(Y,L = 3,optimizer = "fisher",
              ctl = list(minIter = 2,maxIter = 8,tol = 1e-8,penalty = 0))
U <- as.matrix(cbind(out$X,out$offsets,out$factors))
V <- as.matrix(cbind(out$coefX,1,out$loadings))
colnames(U) <- c("intercept","offset",paste0("d",1:3))
colnames(V) <- c("intercept","offset",paste0("d",1:3))
fit0 <- init_glmpca(Y,LL = t(V),FF = t(U),fixed_factors = 1:2,
                    fixed_loadings = 2)

# Fit the Poisson GLM-PCA model by running 40 iterations of the
# glmpca Fisher scoring algorithm.
set.seed(1)
out <- glmpca(Y,L = 3,optimizer = "fisher",
              ctl = list(maxIter = 80,tol = 1e-8,penalty = 0))

# Fit the Poisson GLM-PCA model by running 40 ccd updates.
fit <- fit_glmpca(Y,fit0 = fit0,max_iter = 40,algorithm = "ccd",tol = 1e-12)

# Plot the progress of the two optimization algorithms.
pdat <- rbind(data.frame(method = "glmpca (fisher)",
                         iter   = seq(1,80),
                         loglik = out$lik),
              data.frame(method = "fastglmpca (ccd)",
                         iter   = 1:41,
                         loglik = fit$progress$loglik))
bestloglik <- max(pdat$loglik)
pdat <- transform(pdat,loglik = bestloglik - loglik + 1e-4)
p <- ggplot(pdat,aes(x = iter,y = loglik,color = method)) +
  geom_line(size = 0.75) +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = c("darkblue","darkorange")) +
  labs(y = "distance to best loglik") +
  theme_cowplot(font_size = 12)
print(p)
U <- as.matrix(cbind(out$X,out$offsets,out$factors))
V <- as.matrix(cbind(out$coefX,1,out$loadings))
plot(U,t(fit$FF),pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")
plot(V,t(fit$LL),pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")
