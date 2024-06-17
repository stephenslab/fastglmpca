load("~/Documents/fastglmpca/inst/analysis/results.RData")
res_1core_fastglmpca <- pbmc_purified_results$fastglmpca_1_core$nmi_res_by_iter
res_28core_fastglmpca <- pbmc_purified_results$fastglmpca_28_cores$nmi_res_by_iter
res_glmpca <- pbmc_purified_results$glmpca$nmi_res_by_iter
res_scGBM <- pbmc_purified_results$scGBM$nmi_res_by_iter

res_1core_fastglmpca$hour <- seq(0, 10, length.out = 13)
res_28core_fastglmpca$hour <- seq(0, 10, length.out = 220)
res_glmpca$hour <- seq(0, 10, length.out = 247)
res_scGBM$hour <- seq(0, 10, length.out = 82)

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
  method = c(
    rep("fastglmpca-1core", 13),
    rep("fastglmpca-28core", 220),
    rep("glmpca-avagrad", 247),
    rep("scGBM", 82)
  )
)

library(ggplot2)

ggplot(data = res_df, aes(x = hour, y = nmi)) +
  geom_line(aes(color = method)) +
  xlab("Time (Hours)") +
  ylab("NMI") +
  cowplot::theme_cowplot()
