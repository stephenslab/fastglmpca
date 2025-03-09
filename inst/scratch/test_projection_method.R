library(fastglmpca)
cc <- pbmc_facs$counts[Matrix::rowSums(pbmc_facs$counts) > 10, ]

set.seed(1)
# leave out a single sample from fitting
fit1 <- fit_glmpca_pois(
  Y = cc, 
  K = 2,
  control = list(training_frac = 0.75, maxiter = 100)
)

set.seed(1)
# leave out a single sample from fitting
fit2 <- fit_glmpca_pois(
  Y = cc, 
  K = 2,
  control = list(training_frac = 1, maxiter = 100)
)

# for some reason the calculated log-likelihood and the expected
# are not matching up
set.seed(1)
fit2 <- fit_glmpca_pois(
  Y = pbmc_facs$counts, 
  K = 2,
  control = list(training_frac = 0.25, maxiter = 10, num_projection_ccd_iter = 25)
)

set.seed(1)
fit3 <- fit_glmpca_pois(
  Y = pbmc_facs$counts, 
  K = 2,
  control = list(training_frac = 0.25, maxiter = 10, num_projection_ccd_iter = 5)
)
# 
# df1 <- data.frame(
#   celltype = pbmc_facs$samples$celltype,
#   PC1 = fit1$V[,1],
#   PC2 = fit1$V[,2]
# )
# 
# library(ggplot2)
# 
# ggplot(data = df1) +
#   geom_point(aes(x = PC1, y = PC2, color = celltype))
# 
# df2 <- data.frame(
#   celltype = pbmc_facs$samples$celltype,
#   PC1 = fit2$V[,1],
#   PC2 = fit2$V[,2]
# )
# 
# library(ggplot2)
# 
# ggplot(data = df2) +
#   geom_point(aes(x = PC1, y = PC2, color = celltype))
