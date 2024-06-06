library(fastglmpca)

set.seed(1)
counts <- pbmc_facs$counts
counts <- counts[sample(1:nrow(counts), 250), ]

counts <- counts[, sample(1:ncol(counts), 125)]

counts <- counts[Matrix::rowSums(counts) > 0, ]
counts <- counts[,Matrix::colSums(counts) > 0]

fit0 <- fit_glmpca_pois(counts, K = 2, control = list(maxiter = 25000), verbose = F)

original_mean <- exp(fit0$U %*% diag(fit0$d) %*% t(fit0$V) +
  fit0$X %*% t(fit0$B) + fit0$W %*% t(fit0$Z))

# now, I want to see if I can project a subset of cells onto the fit
# hold out 5 cells
held_out <- counts[, 1:5]

# remove first 5 cells from fit
fit0$V <- fit0$V[6:125,,drop=FALSE]
fit0$B <- fit0$B[6:125,,drop=FALSE]
fit0$Z <- fit0$Z[6:125,,drop=FALSE]

new_Z <- fit0$Z[1:5,,drop=FALSE]

projected_fit <- fastglmpca:::project_onto_U(
  fit0,
  new_Y = held_out,
  new_Z = new_Z
)


projected_mean <- exp(projected_fit$U %*% diag(projected_fit$d) %*% t(projected_fit$V) +
  projected_fit$X %*% t(projected_fit$B) + projected_fit$W %*% t(projected_fit$Z))

original_mean <- original_mean[,1:5,drop=FALSE]
projected_mean <- projected_mean[,1:5,drop=FALSE]

m_df <- data.frame(
  original = as.vector(original_mean),
  projected = as.vector(projected_mean)
)

plot(m_df$original, m_df$projected)
