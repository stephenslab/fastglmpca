load("~/Documents/fastglmpca/inst/analysis/results.RData")

df <- runtime_scaling_res$ncells

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


df <- runtime_scaling_res$nfactors

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

library(ggpubr)
ggarrange(g1, g2, common.legend = TRUE, legend = "right", labels = "AUTO")
ggsave(
  "~/Documents/fastglmpca/inst/scratch/runtime_scaling.pdf",
  device = "pdf",
  width = 15,
  height = 8
)

