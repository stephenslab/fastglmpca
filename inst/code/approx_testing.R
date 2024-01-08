library(fastglmpca)

set.seed(1)

library(RhpcBLASctl)
omp_set_num_threads(10)
blas_set_num_threads(1)

fit0 <- init_glmpca_pois(
  Y = pbmc_facs$counts,
  K = 2
)

fit_exact <- fit_glmpca_pois(
  Y = pbmc_facs$counts,
  fit0 = fit0,
  control = list(maxiter=100)
)

# want to see where the mean stands
m <- fit_exact$U %*% diag(fit_exact$d) %*% t(fit_exact$V) +
  fit_exact$X %*% t(fit_exact$B) + fit_exact$W %*% t(fit_exact$Z)

m2 <- fit_exact$U %*% diag(fit_exact$d) %*% t(fit_exact$V)

Y_0_idx <- which(as.matrix(pbmc_facs$counts) == 0)

ans <- pracma::polyApprox(exp,a = -7,b = 0,n = 2)
a1  <- ans$p[2]
a2  <- ans$p[1]

set.seed(1)

fit0_approx <- init_glmpca_pois(
  Y = as.matrix(pbmc_facs$counts),
  K = 2
)

fit_approx <- fit_glmpca_pois_approx(
  Y = as.matrix(pbmc_facs$counts),
  fit0 = fit0_approx,
  a1 = a1,
  a2 = a2,
  max_iter = 100
)

df <- data.frame(celltype = pbmc_facs$samples$subpop,
           PC1 = fit_exact$V[,1],
           PC2 = fit_exact$V[,2])

library(ggplot2)

ggplot(df,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  cowplot::theme_cowplot(font_size = 10)

out <- fastglmpca:::orthonormalize(
  t(fit_approx$LL[c(1,2),]),
  t(fit_approx$FF[c(1,2),])
  )

df2 <- data.frame(celltype = pbmc_facs$samples$subpop,
                 PC1 = out$V[,1],
                 PC2 = out$V[,2])

library(ggplot2)

ggplot(df2,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  cowplot::theme_cowplot(font_size = 10)


