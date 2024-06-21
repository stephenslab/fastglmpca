load("~/Documents/fastglmpca/inst/analysis/results.RData")
res_1core_fastglmpca <- pbmc_purified_results$fastglmpca_1_core$cluster_metrics_by_iter
res_28core_fastglmpca <- pbmc_purified_results$fastglmpca_28_cores$cluster_metrics_by_iter
res_glmpca <- pbmc_purified_results$glmpca$cluster_metrics_by_iter
res_scGBM <- pbmc_purified_results$scGBM$cluster_metrics_by_iter

res_1core_fastglmpca$hour <- seq(0, 10, length.out = nrow(res_1core_fastglmpca))
res_28core_fastglmpca$hour <- seq(0, 10, length.out = nrow(res_28core_fastglmpca))
res_glmpca$hour <- seq(0, 10, length.out = nrow(res_glmpca))
res_scGBM$hour <- seq(0, 10, length.out = nrow(res_scGBM))

res_df <- data.frame(
  hour = c(
    res_1core_fastglmpca$hour,
    res_28core_fastglmpca$hour,
    res_glmpca$hour,
    res_scGBM$hour
  ),
  nmi = c(
    res_1core_fastglmpca$nmi,
    res_28core_fastglmpca$nmi,
    res_glmpca$nmi,
    res_scGBM$nmi
  ),
  ari = c(
    res_1core_fastglmpca$ari,
    res_28core_fastglmpca$ari,
    res_glmpca$ari,
    res_scGBM$ari
  ),
  method = c(
    rep("fastglmpca-1core", nrow(res_1core_fastglmpca)),
    rep("fastglmpca-28core", nrow(res_28core_fastglmpca)),
    rep("glmpca-avagrad", nrow(res_glmpca)),
    rep("scGBM", nrow(res_scGBM))
  )
)

library(ggplot2)

g1 <- ggplot(data = res_df, aes(x = hour, y = nmi)) +
  geom_line(aes(color = method)) +
  xlab("Time (Hours)") +
  ylab("NMI") +
  cowplot::theme_cowplot()

g2 <- ggplot(data = res_df, aes(x = hour, y = ari)) +
  geom_line(aes(color = method)) +
  xlab("Time (Hours)") +
  ylab("ARI") +
  cowplot::theme_cowplot()

library(ggpubr)
ggarrange(g1, g2, nrow = 1, labels = "AUTO", common.legend = TRUE, legend = "right")
