library(Matrix)
library(matrixStats)
set.seed(1)
load("~/git/single-cell-topics/data/droplet.RData")
Y <- as.matrix(counts)
s <- colSds(Y)
rows <- sample(7193,2000)
cols <- which(s > 1)
Y <- Y[rows,cols]

# First initialize without extrapolation.
fit0 <- init_glmpca(Y = Y,K = 5,fit_col_size_factor = TRUE,
                    fit_row_intercept = TRUE)
fit0 <- fit_glmpca(Y = Y,fit0 = fit0,tol = 1e-4,algorithm = "ccd",
                   link = "log",max_iter = 4,warmup = FALSE,
                   use_extrapolation = FALSE,
                   control = list(line_search = TRUE,num_iter = 4,
                                  alpha = 1e-4,beta = 0.5))

# Fit the GLM-PCA model *without* extrapolation.
fit <- fit_glmpca(Y = Y,fit0 = fit0,tol = 1e-8,algorithm = "ccd",link = "log",
                  max_iter = 100,warmup = FALSE,use_extrapolation = FALSE,
                  control = list(line_search = TRUE,num_iter = 4,alpha = 1e-4,
                                 beta = 0.5))

# Fit the GLM-PCA model *with* extrapolation.
fit_extra <- fit_glmpca(Y = Y,fit0 = fit0,tol = 1e-4,algorithm = "ccd",
                         link = "log",max_iter = 100,warmup = FALSE,
                         use_extrapolation = TRUE,
                         control = list(line_search = TRUE,num_iter = 4,
                                        alpha = 1e-4,beta = 0.5))

y <- max(c(fit$progress$loglik,fit_extra$progress$loglik))
plot(y - fit$progress$loglik + 1,type = "l",lwd = 2,col = "darkblue",
     xlab = "iteration",ylab = "loglik",log = "y",ylim = c(1,1e6))
lines(y - fit_extra$progress$loglik + 1,lwd = 2,col = "darkorange",
      lty = "dashed")

# plot(fit0$FF,fit$FF,pch = 20)
# plot(fit0$LL,fit$LL,pch = 20)
