# This is a script used to verify implementation of fit_glmpca in this
# package against glmpca.

# remotes::install_github("eweine/glmpca",upgrade = "never")
library(glmpca)
library(ggplot2)
library(cowplot)

# Simulate a 200 x 500 counts matrix.
set.seed(1)
dat <- generate_glmpca_data(200,500,K = 3)
Y <- dat$Y

# Generate an initial fit of the Poisson GLM-PCA by running 10
# iterations of the glmpca "Fisher scoring" algorithm.
set.seed(1)
out <- glmpca(Y,L = 3,optimizer = "fisher",
              ctl = list(minIter = 2,maxIter = 10,tol = 1e-8,penalty = 0))
loglik0 <- out$lik
U <- as.matrix(cbind(out$X,out$offsets,out$factors))
V <- as.matrix(cbind(out$coefX,1,out$loadings))
colnames(U) <- c("intercept","offset",paste0("d",1:3))
colnames(V) <- c("intercept","offset",paste0("d",1:3))
fit0 <- init_glmpca(Y,LL = t(V),FF = t(U),fixed_factors = 1:2,
                    fixed_loadings = 2)

# Fit the Poisson GLM-PCA model by running 80 iterations of the
# glmpca Fisher scoring algorithm.
set.seed(1)
out <- glmpca(Y,L = 3,optimizer = "fisher",
              ctl = list(maxIter = 80,tol = 1e-8,penalty = 0))

# Fit the Poisson GLM-PCA model by running 70 ccd updates.
fit <- fit_glmpca(Y,fit0 = fit0,max_iter = 70,algorithm = "ccd",tol = 1e-15)

# Check that the glmpca log-likelihood and deviance calculations agree.
k0 <- -sum(lgamma(Y + 1))
k1 <- 2*sum(Y*(log(Y + 1e-15) - 1))
ll <- k0 - (out$dev - k1)/2
print(range(ll - out$lik))

# Check that the estimates returned by glmpca and fit_glmpca mostly
# agree.
U <- as.matrix(cbind(out$X,out$offsets,out$factors))
V <- as.matrix(cbind(out$coefX,1,out$loadings))
plot(U,t(fit$FF),pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")
plot(V,t(fit$LL),pch = 20)
abline(a = 0,b = 1,lty = "dashed",col = "magenta")

# Plot the progress of the two optimization algorithms.
pdat <- rbind(data.frame(method = "glmpca (fisher)",
                         iter   = seq(1,80),
                         loglik = out$lik),
              data.frame(method = "fastglmpca (ccd)",
                         iter   = seq(1,81),
                         loglik = c(loglik0,fit$progress$loglik)))
bestloglik <- max(pdat$loglik)
pdat <- transform(pdat,loglik = bestloglik - loglik + 1e-6)
p <- ggplot(pdat,aes(x = iter,y = loglik,color = method)) +
  geom_line(size = 0.75) +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = c("darkblue","darkorange")) +
  labs(y = "distance to best loglik") +
  theme_cowplot(font_size = 12)
print(p)

