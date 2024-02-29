# TO DO: Explain what this script is for.
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
ll1 <- loglik_fastglmpca(fit1,counts)
print(fit1$loglik - ll1)
