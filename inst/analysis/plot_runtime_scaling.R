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

ggplot(data = df_long) +
  geom_line(aes(x = ncells, y = time, color = algorithm)) +
  cowplot::theme_cowplot()
