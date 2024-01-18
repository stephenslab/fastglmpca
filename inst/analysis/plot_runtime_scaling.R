df <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/num_cells_scaling_df.rds"
)

df_long <- data.frame(
  ncells = rep(df$ncells, 3),
  time = c(df$fastglmpca_time, df$glmpca_time, df$scGBM_time),
  algorithm = c(
    rep("fastglmpca", nrow(df)),
    rep("glmpca-sgd", nrow(df)),
    rep("scGBM", nrow(df))
  )
)

library(ggplot2)

g1 <- ggplot(data = df_long) +
  geom_point(aes(x = ncells, y = time, color = algorithm)) +
  geom_line(aes(x = ncells, y = time, color = algorithm)) +
  cowplot::theme_cowplot() +
  xlab("Number of Cells") +
  ylab("Time (s)")


df <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/num_factors_scaling_df.rds"
)

df_long <- data.frame(
  K = rep(df$K, 3),
  time = c(df$fastglmpca_time, df$glmpca_time, df$scGBM_time),
  algorithm = c(
    rep("fastglmpca", nrow(df)),
    rep("glmpca-sgd", nrow(df)),
    rep("scGBM", nrow(df))
  )
)

g2 <- ggplot(data = df_long) +
  geom_point(aes(x = K, y = time, color = algorithm)) +
  geom_line(aes(x = K, y = time, color = algorithm)) +
  cowplot::theme_cowplot() +
  xlab("Number of Factors") +
  ylab("Time (s)")

