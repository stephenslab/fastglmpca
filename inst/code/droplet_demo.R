# sinteractive -p mstephens --account=pi-mstephens --time=24:00:00 \
#   -c 10 --mem=16G
# use_daarem     <- FALSE
# orthonormalize <- TRUE
# outfile        <- "droplet_k=5_daarem=no.RData"
use_daarem     <- TRUE
orthonormalize <- FALSE
outfile        <- "droplet_k=5_daarem=yes.RData"
library(Matrix)
library(RhpcBLASctl)
library(fastglmpca)
load(file.path("/project2/mstephens/pcarbo/git/fastTopics-experiments/data",
               "droplet.RData"))
omp_set_num_threads(10)
blas_set_num_threads(1)
set.seed(1)
counts <- t(counts)
cat("use_daarem     =",use_daarem,"\n")
cat("orthonormalize =",orthonormalize,"\n")
cat("outfile        =",outfile,"\n")
fit0 <- fit_glmpca_pois(counts,K = 5,
                        control = list(maxiter = 100,tol = 0,
                                       use_daarem = FALSE,
                                       orthonormalize = TRUE))
fit <- fit_glmpca_pois(counts,fit0 = fit0,
                       control = list(maxiter = 100,tol = 0,
                                      use_daarem = use_daarem,
                                      orthonormalize = orthonormalize))
save(list = "fit",file = outfile)

# Create the plot comparing the two runs:
#
# library(ggplot2)
# library(cowplot)
# load("droplet_k=5_daarem=no.RData")
# fit1 <- fit
# load("droplet_k=5_daarem=yes.RData")
# fit2 <- fit
# pdat <- rbind(data.frame(method = "fpiter",
#                          iter   = fit1$progress$iter,
#                          loglik = fit1$progress$loglik),
#               data.frame(method = "daarem",
#                          iter   = fit2$progress$iter,
#                          loglik = fit2$progress$loglik))
# best_loglik <- max(pdat$loglik)
# pdat <- transform(pdat,loglik = best_loglik - loglik + 1)
# p <- ggplot(pdat,aes(x = iter,y = loglik,color = method)) +
#   geom_line(size = 0.5) +
#   scale_color_manual(values = c("darkblue","darkorange"))+
#   scale_y_continuous(trans = "log10",
#                      breaks = c(1,10,100,1e3,1e4,1e5,1e6,1e7)) +
#   theme_cowplot(font_size = 10)
# print(p)
#
# Create plots comparing the estimates:
#
# fit2$U[,c(3,5)] <- -fit2$U[,c(3,5)]
# plot(fit1$U,fit2$U,pch = 20)
# abline(a = 0,b = 1,col = "magenta",lty = "dotted")
# 
# fit2$V[,c(3,5)] <- -fit2$V[,c(3,5)]
# plot(fit1$V,fit2$V,pch = 20)
# abline(a = 0,b = 1,col = "magenta",lty = "dotted")
# 
