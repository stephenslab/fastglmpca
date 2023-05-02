library(Matrix)
library(matrixStats)
set.seed(1)
load("~/git/single-cell-topics/data/droplet.RData")
Y <- as.matrix(counts)
s <- colSds(Y)
cols <- which(s > 1)
rows <- sample(7193,2000)
# cols <- which(s > 1)
Y <- Y[,cols]
Y <- as.matrix(Y)

# First initialize without daarem.
fit0 <- init_glmpca(Y = Y,K = 2,fit_col_size_factor = TRUE,
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

y <- max(c(fit$progress$loglik,fit_daarem$progress$loglik))
plot(y - fit$progress$loglik + 0.01,type = "l",lwd = 2,col = "darkblue",
     xlab = "iteration",ylab = "loglik",log = "y",ylim = c(0.01,1e6))
lines(y - fit_daarem$progress$loglik + 0.01,lwd = 2,col = "darkorange",
      lty = "dashed")
