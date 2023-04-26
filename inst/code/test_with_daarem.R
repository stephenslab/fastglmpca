library(Matrix)
library(matrixStats)
set.seed(1)
load("~/git/single-cell-topics/data/droplet.RData")
Y <- as.matrix(counts)
s <- colSds(Y)
cols <- which(s > 1)
Y <- Y[,cols]

# First initialize without daarem.
fit0 <- init_glmpca(Y = Y,K = 3,fit_col_size_factor = TRUE,
                    fit_row_intercept = TRUE)
fit0 <- fit_glmpca(Y = Y,fit0 = fit0,tol = 1e-4,algorithm = "ccd",
                   link = "log",max_iter = 4,warmup = FALSE,use_daarem = FALSE,
                   control = list(line_search = TRUE,num_iter = 4,
                                  alpha = 1e-4,beta = 0.5))

# Fit the GLM-PCA model *without* daarem.
fit <- fit_glmpca(Y = Y,fit0 = fit0,tol = 1e-8,algorithm = "ccd",link = "log",
                  max_iter = 100,warmup = FALSE,use_daarem = FALSE,
                  control = list(line_search = TRUE,num_iter = 4,alpha = 1e-4,
                                 beta = 0.5))

# Fit the GLM-PCA model *with* daarem.
fit_daarem <- fit_glmpca(Y = Y,fit0 = fit0,tol = 1e-4,algorithm = "ccd",
                         link = "log",max_iter = 100,warmup = FALSE,
                         use_daarem = TRUE,
                         control = list(line_search = TRUE,num_iter = 4,
                                        alpha = 1e-4,beta = 0.5))

plot(fit$progress$loglik,type = "l",lwd = 2,col = "darkblue",
     xlab = "iteration",ylab = "loglik")
lines(fit_daarem$progress$loglik,lwd = 2,col = "darkorange",
      lty = "dashed")
