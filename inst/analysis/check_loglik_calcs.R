# A script to check our log-likelihood calculations for different
# packages (glmpca, scGBM, fastglmpca).
library(Matrix)
library(fastTopics)
source("../code/loglik.R")
load("~/git/fastTopics-experiments/data/droplet.RData")
counts <- t(counts)
fit1 <- readRDS(paste("droplets_fastglmpca_fit_10_factors_5105_iter_28",
                      "cores_dec_23.rds",sep = "_"))
fit2 <- readRDS(paste("droplets_glmpca_fit_10_factors_10_hrs_avagrad",
                      "optimizer_minibatch_stochastic_dec_23.rds",sep = "_"))
fit3 <- readRDS("droplets_scGBM_fit_10_factors_no_beta_infer_10_hrs.rds")
ll1 <- glmpca_pois_loglik(counts,logrates_fastglmpca(fit1))
ll2 <- glmpca_pois_loglik(counts,logrates_glmpca(fit2))
ll3 <- glmpca_pois_loglik(counts,logrates_scgbm(fit3),const = 0)
print((fit1$loglik - ll1)/ll1)
print((max(fit2$lik) - ll2)/ll2)
print((max(fit3$loglik) - ll3)/ll3)
