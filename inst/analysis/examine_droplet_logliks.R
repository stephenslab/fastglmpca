library(Matrix)
library(fastglmpca)
library(ggplot2)
library(cowplot)
source("../code/loglik.R")
load("~/git/fastTopics-experiments/data/droplet.RData")
counts <- t(counts)
fit1 <- readRDS(paste("droplets_fastglmpca_fit_10_factors_5105_iter_28",
                      "cores_dec_23.rds",sep = "_"))
fit2 <- readRDS(paste("droplets_glmpca_fit_10_factors_10_hrs_avagrad",
                      "optimizer_minibatch_stochastic_dec_23.rds",sep = "_"))
fit3 <- readRDS("droplets_scGBM_fit_10_factors_no_beta_infer_10_hrs.rds")
glmpca_pois_loglik <- function (Y, H)
  Y*H - exp(H) - lfactorial(as.matrix(Y))
H1 <- glmpca_pois_loglik(counts,logrates_fastglmpca(fit1))
H2 <- glmpca_pois_loglik(counts,logrates_glmpca(fit2))
H3 <- glmpca_pois_loglik(counts,logrates_scgbm(fit3))

tissue_colors <- c("royalblue",   # basal
                   "firebrick",   # ciliated
                   "forestgreen", # club
                   "gold",        # goblet
                   "darkmagenta", # ionocyte
                   "darkorange",  # neuroendocrine
                   "skyblue")     # tuft

# Compare cell-wise logliks.
pdat <- cbind(samples,
              data.frame(fastglmpca = colSums(H1),
                         glmpca     = colSums(H2),
                         scgbm      = colSums(H3)))
p1 <- ggplot(pdat,aes(x = glmpca,y = fastglmpca,color = tissue)) +
  geom_point(show.legend = FALSE) +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  scale_color_manual(values = tissue_colors) +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = scgbm,y = fastglmpca,color = tissue)) +
  geom_point(show.legend = FALSE) +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  scale_color_manual(values = tissue_colors) +
  theme_cowplot(font_size = 12)

# Compare gene-wise logliks.
pdat <- data.frame(fastglmpca = rowSums(H1),
                   glmpca     = rowSums(H2),
                   scgbm      = rowSums(H3))
pdat <- subset(pdat,scgbm > -1e5)
p3 <- ggplot(pdat,aes(x = glmpca,y = fastglmpca)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  theme_cowplot(font_size = 12)
p4 <- ggplot(pdat,aes(x = scgbm,y = fastglmpca)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",linetype = "dotted") +
  theme_cowplot(font_size = 12)

plot_grid(p1,p2,p3,p4)
