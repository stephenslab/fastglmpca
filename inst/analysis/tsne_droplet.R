library(Matrix)
library(uwot)
library(Rtsne)
library(ggplot2)
library(cowplot)
set.seed(1)
load("results.RData")
load("~/git/fastTopics-experiments/data/droplet.RData")
# V1 <- droplets_res_list$fastglmpca_28_cores$PCs_10$V
V1 <- droplets_res_list$scGBM[["10_factors"]]$V
# V1 <- droplets_res_list$glmpca[["10_factors"]]$V

out <- Rtsne(V1,dims = 2,pca = FALSE,normalize = FALSE,eeta = 200,
             perplexity = 100,theta = 0.1,max_iter = 1000,
             check_duplicates = FALSE,verbose = TRUE,
             num_threads = 4)
colnames(out$Y) <- c("d1","d2")

# out <- umap(V1,n_components = 2,n_neighbors = 30,
#             metric = "euclidean", scale = FALSE,
#             pca = NULL,n_threads = 4,verbose = TRUE)
# colnames(out) <- c("d1","d2")

tissue_colors <- c("royalblue",   # basal
                   "firebrick",   # ciliated
                   "forestgreen", # club
                   "gold",        # goblet
                   "darkmagenta", # ionocyte
                   "darkorange",  # neuroendocrine
                   "skyblue")     # tuft

pdat <- cbind(samples,out$Y)
p1 <- ggplot(pdat,aes(x = d1,y = d2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  theme_cowplot(font_size = 12)
print(p1)
