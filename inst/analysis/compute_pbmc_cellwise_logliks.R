# sinteractive -p mstephens -c 4 --mem=64G --time=2:00:00
# module load R/3.5.1
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
source("../code/loglik.R")
datadir <- "/project2/mstephens/eweine/fastglmpca_experiments"
load(file.path(datadir,"pbmc_68k.RData"))
T_celltypes <-
  c("CD4+ T Helper2","CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T",
    "CD4+/CD45RO+ Memory","CD8+/CD45RA+ Naive Cytotoxic")
samples$celltype <- as.character(samples$celltype)
samples$celltype <- with(samples,
                         if_else(celltype %in% T_celltypes,"T Cell",celltype))
samples$celltype <- as.factor(samples$celltype)
counts <- t(counts)
fit1 <- readRDS(file.path(datadir,
                          paste("pbmc_fastglmpca_fit_10_factors_450_iter",
                                "28_cores_dec_23.rds",sep = "_")))
fit2 <- readRDS(file.path(datadir,
                          paste("pbmc_glmpca_fit_10_factors_10_hrs_avagrad",
                                "optimizer_minibatch_stochastic_dec_23.rds",
                                sep = "_")))
fit3 <- readRDS(file.path(datadir,
                          paste("pbmc_scGBM_fit_10_factors_no_beta_infer",
                                "10_hrs.rds",sep = "_")))
glmpca_pois_loglik <- function (Y, H)
  Y*H - exp(H) - lfactorial(as.matrix(Y))
H1 <- glmpca_pois_loglik(counts,logrates_fastglmpca(fit1))
H2 <- glmpca_pois_loglik(counts,logrates_glmpca(fit2))
H3 <- glmpca_pois_loglik(counts,logrates_scgbm(fit3))

pbmc_colors <- c("forestgreen", # CD14+
                 "dodgerblue",  # B cells
                 "limegreen",   # CD34+
                 "gray",        # NK cells
                 "tomato",      # cytotoxic T cells
                 "darkmagenta", # dendritic
                 "gold")        # T cells

# Compare the cell-wise logliks.
pdat <- cbind(samples,
              data.frame(fastglmpca = colSums(H1),
                         glmpca     = colSums(H2),
                         scgbm      = colSums(H3)))
x <- with(pdat,c(scgbm,glmpca,fastglmpca))
p1 <- ggplot(pdat,aes(x = glmpca,y = fastglmpca,color = celltype)) +
  geom_point(show.legend = FALSE) +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  scale_color_manual(values = pbmc_colors) +
  xlim(min(x),max(x)) +
  ylim(min(x),max(x)) +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = scgbm,y = fastglmpca,color = celltype)) +
  geom_point(show.legend = FALSE) +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dotted") +
  scale_color_manual(values = pbmc_colors) +
  xlim(min(x),max(x)) +
  ylim(min(x),max(x)) +
  theme_cowplot(font_size = 12)
ggsave("pbmc_loglik_cellwise1.eps",p1,height = 3,width = 3.25)
ggsave("pbmc_loglik_cellwise2.eps",p2,height = 3,width = 3.25)
