# sinteractive -p mstephens -c 4 --mem=16G --time=2:00:00
# module load R/3.5.1
library(Matrix)
library(ggplot2)
library(cowplot)
source("../code/loglik.R")
datadir <- "/project2/mstephens/eweine/fastglmpca_experiments"
load(file.path(datadir,"droplet.RData"))
counts <- t(counts)
fit1 <- readRDS(file.path(datadir,
                          paste("droplets_fastglmpca_fit_10_factors_5105",
                                "iter_28_cores_dec_23.rds",sep = "_")))
fit2 <- readRDS(file.path(datadir,
                          paste("droplets_glmpca_fit_10_factors_10_hrs",
                                "avagrad_optimizer_minibatch_stochastic",
                                "dec_23.rds",sep = "_")))
fit3 <- readRDS(file.path(datadir,
                          paste("droplets_scGBM_fit_10_factors_no_beta",
                                "infer_10_hrs.rds",sep = "_")))
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

# Compare the cell-wise logliks.
pdat <- cbind(samples,
              data.frame(fastglmpca = colSums(H1),
                         glmpca     = colSums(H2),
                         scgbm      = colSums(H3)))
x <- with(pdat,c(scgbm,glmpca,fastglmpca))
p1 <- ggplot(pdat,aes(x = glmpca,y = fastglmpca,color = tissue)) +
  geom_point(show.legend = FALSE) +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  scale_color_manual(values = tissue_colors) +
  xlim(min(x),max(x)) +
  ylim(min(x),max(x)) +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = scgbm,y = fastglmpca,color = tissue)) +
  geom_point(show.legend = FALSE) +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  scale_color_manual(values = tissue_colors) +
  xlim(min(x),max(x)) +
  ylim(min(x),max(x)) +
  theme_cowplot(font_size = 12)
ggsave("droplet_loglik_cellwise1.eps",p1,height = 3,width = 3.25)
ggsave("droplet_loglik_cellwise2.eps",p2,height = 3,width = 3.25)
