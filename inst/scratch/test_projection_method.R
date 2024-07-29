library(fastglmpca)

set.seed(1)
fit1 <- fit_glmpca_pois(
  Y = pbmc_facs$counts, 
  K = 2,
  control = list(training_frac = 1, maxiter = 10)
)

# for some reason the calculated log-likelihood and the expected
# are not matching up
set.seed(1)
fit2 <- fit_glmpca_pois(
  Y = pbmc_facs$counts, 
  K = 2,
  control = list(training_frac = 0.99, maxiter = 10)
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
